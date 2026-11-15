// vres: "part of the key is known" -- the rest of the key is solved from the key-schedule outputs (S, b0b1 range, k780, k7b0).
//
// A hypothesis is the list of the 4-byte chunks of the key, in key order:
//    K hex len        known chunk (len = 4, or 1..3 for the last, partial chunk)
//    U                unknown full chunk (any bytes; with 2 or 3 U chunks their 4th bytes must be >= 0x33 for L = 16, >= 0x20 for L = 17)
//    M hex mask len   chunk with some unknown bytes: bit i of mask = byte i unknown; those bytes run over [blo, bhi] (outer loop)
// For every assignment of the M bytes the known contributions are divided out of the targets and the U chunks (0..3) are solved:
//    0 U: direct check          1 U: the candidate values of sigma are the chunk
//    2/3 U: the two sum-product systems on the 2^-L lattice (f3 side from k7b0 first, then the f1 side from k780 for the sigma values that
//           survive), then the chunks are rebuilt from the pairs of lattice values (points of a 3-D lattice in a box) and their sum is
//           compared with sigma
// Every candidate is checked with the complete key schedule (map of the last chunk, S, b0b1 range, k780 (any branch), k7b0, and
// k788/k7b8 when given).  Output: one line "KEY id assignment hexkey S b0b1 k780 k788 k7b0 k7b8" per ordering of the U chunks that
// gives the wanted map, and one "STAT ..." line per hypothesis.
//
// input file:
//    T S Tlo Thi map c3 h3 nbr c1a h1a [c1b h1b] has2 [k788 h2 k7b8 h4]     (hex S/Tlo/Thi; c/h = centre and half-width of the cells)
//    O L blo bhi sjmax                                                 (lattice bits 16|17, byte range of the M bytes, max sum of J)
//    H id nchunks        followed by nchunks chunk lines
// build: gcc -O2 -march=native -ffp-contract=off -fno-fast-math -fopenmp -o vres vres.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>
typedef unsigned __int128 u128;
#define M54 ((1ULL << 54) - 1)
#define TWO54 18014398509481984.0

static double K_PI, K_E, C01; static uint64_t MP, ME;
static double h2d(uint64_t h) { double d; memcpy(&d, &h, 8); return d; }
static inline double clr(double a) { uint64_t q; memcpy(&q, &a, 8); q &= ~0xFFULL; memcpy(&a, &q, 8); return a; }
static inline uint32_t ftol(double x) { if (x != x || isinf(x) || fabs(x) >= 9.223372036854775808e18) return 0; return (uint32_t)(int64_t)x; }
static inline int64_t isqrt64(int64_t D) { int64_t q = (int64_t)sqrt((double)D); while (q * q > D) q--; while ((q + 1) * (q + 1) <= D) q++; return q; }
// squares mod 64: 0 1 4 9 16 17 25 33 36 41 49 57
static const uint64_t SQ64 = (1ULL<<0)|(1ULL<<1)|(1ULL<<4)|(1ULL<<9)|(1ULL<<16)|(1ULL<<17)|(1ULL<<25)|(1ULL<<33)|(1ULL<<36)|(1ULL<<41)|(1ULL<<49)|(1ULL<<57);

// ---- target and options
static uint32_t S_t, Tlo, Thi; static int wantmap, nbr, has2; static double c3, h3, c1[2], h1[2], k788t, h2t, k7b8t, h4t;
static int LAT, blo, bhi, nb; static int64_t sjmax; static int64_t LO, HI; static uint32_t CMIN;
static double WIDEN = 1.25, WADD = 5e-13; static int XFIRST = 0;      // the cells are widened a little for the search; the caller decides with the keystream

typedef struct { double f[4]; uint32_t s, t; int map; } CV;
static inline CV cvals(uint32_t c) {
    CV v; double a1 = clr((double)c * K_PI), i1 = floor(a1), a2 = clr(i1 * K_PI), i2 = floor(a2);
    double a3 = clr((double)c * K_E), i3 = floor(a3), a4 = clr(i3 * K_E), i4 = floor(a4);
    v.f[0] = a1 - i1; v.f[1] = a2 - i2; v.f[2] = a3 - i3; v.f[3] = a4 - i4;
    v.s = ftol(i3 / i2) + ftol(i1); v.t = ftol(i1 / i4) + ftol(i3);
    double xm = i2; while (xm > 32767.0) xm *= C01;
    int mv = (int16_t)(ftol(xm) & 0xFFFF); v.map = mv % 3;
    return v;
}
// value on the 2^-LAT lattice (2^LAT .. 2^(LAT+1)-1) and class; 0 if the class is finer than LAT
static inline uint32_t latv(uint32_t c, double K, int *ell) {
    double a = clr((double)c * K); uint64_t q; memcpy(&q, &a, 8);
    int e = (int)((q >> 52) & 0x7FF) - 1023; *ell = 44 - e;
    if (*ell > LAT || *ell < 1) return 0;
    double f = a - floor(a);
    return ((1u << *ell) + (uint32_t)(f * (double)(1u << *ell))) << (LAT - *ell);
}

// ---- hypotheses
typedef struct { char type; uint32_t base, mask; int len; } CH;
typedef struct { char id[64]; int nch; CH *ch; int nU, nM, E; int64_t nassign; int lastU; } HYP;
static HYP *hyps; static int nhyp;

// ---- three-gap stepping tables for the classes 14 .. LAT
static const uint64_t CLB[5] = { 0x100000000ULL, 0xC90FDAA2ULL, 0x6487ED51ULL, 0x3243F6A8ULL, 0x1921FB54ULL };   // lower bounds of classes 13(none),14,15,16,17
static uint64_t g_up[4], d_up[4], g_dn[4], d_dn[4], WW[4], F1LO[4], F1HI[4]; static int NCLS;
static void init_steps(void) {
    NCLS = LAT - 13;
    for (int k = 0; k < NCLS; k++) {
        int ell = 14 + k; uint64_t w = 1ULL << (54 - ell); WW[k] = w + (w >> 8);
        F1HI[k] = CLB[k]; F1LO[k] = CLB[k + 1]; if (F1LO[k] < CMIN) F1LO[k] = CMIN;
        g_up[k] = g_dn[k] = 0; uint64_t t = 0;
        for (uint64_t g = 1; !(g_up[k] && g_dn[k]); g++) {
            t = (t + MP) & M54;
            if (!g_up[k] && t < WW[k]) { g_up[k] = g; d_up[k] = t; }
            if (!g_dn[k] && t > (M54 + 1) - WW[k]) { g_dn[k] = g; d_dn[k] = (M54 + 1) - t; }
        }
    }
}

// ---- chunks with a given pair of lattice values (f1 value, f3 value): points of a 3-D lattice in a box
// c in [clo, chi), (c*MP mod 2^54) in [A1, A1+W1), (c*ME mod 2^54) in [A3, A3+W3)  <=>  the lattice generated by (1, MP, ME), (0, 2^54, 0),
// (0, 0, 2^54) has a point in the box.  The basis is reduced once per pair of classes (the coordinates are scaled so that the box is roughly a cube);
// a query lists the integer combinations of the reduced basis around the centre of the box.
typedef __int128 i128;
static const uint64_t CLB3[5] = { 0x100000000ULL, 0xADF85458ULL, 0x56FC2A2CULL, 0x2B7E1516ULL, 0x15BF0A8BULL };   // f3: lower bounds of classes 13(none),14,15,16,17
typedef struct { int ell1, ell3, sh0, sh1, sh3; int64_t clo, chi; int64_t B[3][3]; long double Bi[3][3]; } CP;
static CP cps[16]; static int ncp;
static void lll3(i128 b[3][3]) {
    int k = 1, guard = 0;
    while (k < 3 && guard++ < 100000) {
        long double bs[3][3], mu[3][3], Bn[3];
        for (int i = 0; i <= k; i++) {
            for (int c = 0; c < 3; c++) bs[i][c] = (long double)b[i][c];
            for (int j = 0; j < i; j++) {
                long double d = 0; for (int c = 0; c < 3; c++) d += (long double)b[i][c] * bs[j][c];
                mu[i][j] = d / Bn[j]; for (int c = 0; c < 3; c++) bs[i][c] -= mu[i][j] * bs[j][c];
            }
            Bn[i] = 0; for (int c = 0; c < 3; c++) Bn[i] += bs[i][c] * bs[i][c];
        }
        int changed = 0;
        for (int j = k - 1; j >= 0; j--) {
            long double d = 0, n = 0; for (int c = 0; c < 3; c++) { d += (long double)b[k][c] * (long double)b[j][c]; }
            // coefficient of b[j] in b[k] with respect to the Gram-Schmidt vectors: recompute mu[k][j] from the current b[k]
            long double m = 0; { long double t[3]; for (int c = 0; c < 3; c++) t[c] = (long double)b[k][c]; for (int c = 0; c < 3; c++) m += t[c] * bs[j][c]; m /= Bn[j]; }
            (void)d; (void)n;
            long double q = roundl(m);
            if (q != 0) { i128 qi = (i128)q; for (int c = 0; c < 3; c++) b[k][c] -= qi * b[j][c]; changed = 1; }
        }
        if (changed) continue;                      // recompute the Gram-Schmidt data with the reduced b[k]
        if (Bn[k] >= (0.99L - mu[k][k - 1] * mu[k][k - 1]) * Bn[k - 1]) k++;
        else { for (int c = 0; c < 3; c++) { i128 t = b[k][c]; b[k][c] = b[k - 1][c]; b[k - 1][c] = t; } k = k > 1 ? k - 1 : 1; }
    }
}
static void init_cp(void) {
    ncp = 0;
    for (int k1 = 0; k1 < NCLS; k1++) for (int k3 = 0; k3 < NCLS; k3++) {
        int64_t lo1 = (int64_t)CLB[k1 + 1], hi1 = (int64_t)CLB[k1], lo3 = (int64_t)CLB3[k3 + 1], hi3 = (int64_t)CLB3[k3];
        int64_t lo = (lo1 > lo3 ? lo1 : lo3) - 8, hi = (hi1 < hi3 ? hi1 : hi3) + 8;
        if (lo < (int64_t)CMIN) lo = CMIN;
        if (hi > 0x100000000LL) hi = 0x100000000LL;
        if (hi - lo < 64) continue;
        CP *cp = &cps[ncp++]; cp->ell1 = 14 + k1; cp->ell3 = 14 + k3; cp->clo = lo; cp->chi = hi;
        int mn = cp->ell1 < cp->ell3 ? cp->ell1 : cp->ell3; cp->sh1 = cp->ell1 - mn; cp->sh3 = cp->ell3 - mn;
        double side = ldexp(1.0, 54 - mn); cp->sh0 = (int)floor(log2(side / (double)(hi - lo)) + 0.5); if (cp->sh0 < 0) cp->sh0 = 0;
        i128 b[3][3] = { { (i128)1 << cp->sh0, (i128)MP << cp->sh1, (i128)ME << cp->sh3 }, { 0, (i128)1 << (54 + cp->sh1), 0 }, { 0, 0, (i128)1 << (54 + cp->sh3) } };
        lll3(b);
        for (int i = 0; i < 3; i++) for (int c = 0; c < 3; c++) { if (b[i][c] > ((i128)1 << 62) || b[i][c] < -((i128)1 << 62)) { fprintf(stderr, "lattice basis too large\n"); exit(3); } cp->B[i][c] = (int64_t)b[i][c]; }
        long double m[3][3]; for (int i = 0; i < 3; i++) for (int c = 0; c < 3; c++) m[i][c] = (long double)cp->B[i][c];
        long double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        for (int i = 0; i < 3; i++) for (int c = 0; c < 3; c++) {
            int i1 = (i + 1) % 3, i2 = (i + 2) % 3, c1_ = (c + 1) % 3, c2 = (c + 2) % 3;
            cp->Bi[c][i] = (m[i1][c1_] * m[i2][c2] - m[i1][c2] * m[i2][c1_]) / det;       // inverse = adjugate / det (cyclic cofactors carry the sign)
        }
    }
}
// chunks of one class pair with the residues in [lo1, hi1] and [lo3, hi3] (inclusive); returns the number found (at most maxout)
static inline int cp_enum(const CP *cp, int64_t lo1, int64_t hi1, int64_t lo3, int64_t hi3, uint32_t *out, int maxout, uint64_t *ntest) {
    long double ctr[3] = { ldexpl(((long double)cp->clo + (long double)(cp->chi - 1)) * 0.5L, cp->sh0), ldexpl(((long double)lo1 + (long double)hi1) * 0.5L, cp->sh1), ldexpl(((long double)lo3 + (long double)hi3) * 0.5L, cp->sh3) };
    long double hw[3] = { ldexpl((long double)(cp->chi - 1 - cp->clo) * 0.5L, cp->sh0), ldexpl((long double)(hi1 - lo1) * 0.5L, cp->sh1), ldexpl((long double)(hi3 - lo3) * 0.5L, cp->sh3) };
    int64_t zl[3], zh[3];
    for (int i = 0; i < 3; i++) {
        long double t = 0, d = 0; for (int j = 0; j < 3; j++) { t += ctr[j] * cp->Bi[j][i]; d += hw[j] * fabsl(cp->Bi[j][i]); }
        d += 1e-6L + fabsl(t) * 1e-15L; zl[i] = (int64_t)ceill(t - d); zh[i] = (int64_t)floorl(t + d);
        if (zh[i] < zl[i]) return 0;
    }
    i128 X0lo = (i128)cp->clo << cp->sh0, X0hi = (i128)(cp->chi - 1) << cp->sh0, X1lo = (i128)lo1 * ((i128)1 << cp->sh1), X1hi = (i128)hi1 * ((i128)1 << cp->sh1), X3lo = (i128)lo3 * ((i128)1 << cp->sh3), X3hi = (i128)hi3 * ((i128)1 << cp->sh3);
    int n = 0;
    for (int64_t z0 = zl[0]; z0 <= zh[0]; z0++) for (int64_t z1 = zl[1]; z1 <= zh[1]; z1++) {
        i128 p0 = (i128)z0 * cp->B[0][0] + (i128)z1 * cp->B[1][0], p1 = (i128)z0 * cp->B[0][1] + (i128)z1 * cp->B[1][1], p3 = (i128)z0 * cp->B[0][2] + (i128)z1 * cp->B[1][2];
        for (int64_t z2 = zl[2]; z2 <= zh[2]; z2++) {
            (*ntest)++;
            i128 x0 = p0 + (i128)z2 * cp->B[2][0]; if (x0 < X0lo || x0 > X0hi) continue;
            i128 x1 = p1 + (i128)z2 * cp->B[2][1]; if (x1 < X1lo || x1 > X1hi) continue;
            i128 x3 = p3 + (i128)z2 * cp->B[2][2]; if (x3 < X3lo || x3 > X3hi) continue;
            if (n < maxout) out[n++] = (uint32_t)(x0 >> cp->sh0);
        }
    }
    return n;
}

// ---- per-thread work space
typedef struct { int64_t A, sig; int J, SJ, br; } Cand;
typedef struct { int64_t sig; int J, SJ, br; uint32_t N[3]; } Sol;
typedef struct { int64_t sig; int j3, w; uint64_t fr3; } SigC;
typedef struct { uint32_t c; int next; } PE;
typedef struct {
    Cand *cx, *cy; size_t ncx, ncy, capx, capy;
    Sol *sx, *sy; size_t nsx, nsy, capsx, capsy;
    SigC *sg, *ds; size_t nsg, capsg, capds;
    uint32_t *xv, *yv; int *head; size_t caphead; PE *pool; size_t cappool;
    char *out; size_t lout, capout; void *y2; size_t capy2;
    uint64_t n_assign, n_feas, n_xsol, n_ysol, n_grp, n_cand, n_keys, n_iter, n_fib;
} WS;
#define PUSH(arr, n, cap, T) do { if ((n) == (cap)) { (cap) = (cap) ? (cap) * 2 : 1024; (arr) = realloc((arr), (cap) * sizeof(T)); } } while (0)
static void outs(WS *w, const char *s, int l) { if (w->lout + l + 1 > w->capout) { w->capout = (w->capout + l) * 2 + 4096; w->out = realloc(w->out, w->capout); } memcpy(w->out + w->lout, s, l); w->lout += l; }

static int cmp_cand(const void *a, const void *b) { const Cand *x = a, *y = b;
    if (x->br != y->br) return x->br - y->br; if (x->J != y->J) return x->J - y->J;
    return x->A < y->A ? -1 : x->A > y->A ? 1 : (x->sig < y->sig ? -1 : x->sig > y->sig); }
static int cmp_sigc(const void *a, const void *b) { const SigC *x = a, *y = b; return x->sig < y->sig ? -1 : x->sig > y->sig; }
static int cmp_sol(const void *a, const void *b) { const Sol *x = a, *y = b; return x->sig < y->sig ? -1 : x->sig > y->sig; }
static int cmp_u32(const void *a, const void *b) { uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b; return x < y ? -1 : x > y; }
static inline void bounds(double F, int n, double *lo, double *hi) {
    double j = floor(F), b = 1.0 + (F - j); *lo = j < 0.5 ? b : j < 1.5 ? 2.0 * b : j < 2.5 ? 4.0 * b : ldexp(b, (int)j);
    double m = 1.0 + F / n; *hi = n == 2 ? m * m : n == 3 ? m * m * m : pow(m, n);
}

// candidates of sigma_U (sum of the unknown full chunks) allowed by S and the b0b1 range
static void sigma_cands(WS *w, uint32_t Sp, uint32_t TK, int nU) {
    w->nsg = 0;
    for (int wr = 0; wr <= nU; wr++) {
        int64_t Sw = (int64_t)Sp + ((int64_t)wr << 32);
        int64_t lo = (int64_t)floor((double)(Sw - 3 * nU - sjmax) / K_PI) - 3, hi = (int64_t)floor((double)(Sw - 3 * nU - 4 + nU) / K_PI) + 3;
        { int64_t cm = nU == 1 ? 0x01000000LL : (int64_t)nU * CMIN; if (lo < cm) lo = cm; }     // 1 U: no lattice is used, any 4th byte >= 1
        if (hi > (int64_t)nU * 0xFFFFFFFFLL) hi = (int64_t)nU * 0xFFFFFFFFLL;
        if (hi < lo) continue;
        int64_t vlo = (int64_t)(((u128)ME * (u128)lo) >> 54), vhi = (int64_t)(((u128)ME * (u128)hi) >> 54);
        for (int j3 = 0; j3 < nU; j3++) {
            int64_t off = (2 * nU + (int64_t)TK - j3) & 0xFFFF;
            int64_t x = (vlo + off) & 0xFFFF;
            int64_t v = vlo + (x <= (int64_t)Thi ? (x < (int64_t)Tlo ? (int64_t)Tlo - x : 0) : 65536 - x + (int64_t)Tlo);
            while (v <= vhi) {
                int64_t t16 = (v + off) & 0xFFFF;
                if (t16 < (int64_t)Tlo || t16 > (int64_t)Thi) { v += (t16 > (int64_t)Thi) ? 65536 - t16 + (int64_t)Tlo : (int64_t)Tlo - t16; continue; }
                int64_t vend = v + ((int64_t)Thi - t16); if (vend > vhi) vend = vhi;          // the run of v inside this occurrence of the window
                int64_t sg = (int64_t)((((u128)v << 54) + ME - 1) / ME); u128 pp = (u128)ME * (u128)sg;
                for (; (int64_t)(pp >> 54) <= vend; sg++, pp += ME) {
                    if (sg < lo || sg > hi) continue;
                    PUSH(w->sg, w->nsg, w->capsg, SigC);
                    w->sg[w->nsg].sig = sg; w->sg[w->nsg].j3 = j3; w->sg[w->nsg].w = wr; w->sg[w->nsg].fr3 = (uint64_t)(pp & M54); w->nsg++;
                }
                v = vend + 1;
            }
        }
    }
}

// one run of candidates with the same (br, J): n = 3 -> the smallest element is enumerated; n = 2 -> closed directly
static void enum_run(WS *w, char sys, int n, const Cand *l, int nl, long double B, double tol) {
    if (n == 2) {
        long double Bf = ldexpl(B, 2 * LAT);
        double preal = (double)Bf, wd = preal * tol + 2e-5;
        int64_t rlo = (int64_t)ceil(preal - wd), rhi = (int64_t)floor(preal + wd);
        for (int a = 0; a < nl; a++) for (int64_t r = rlo; r <= rhi; r++) {
            int64_t s = l[a].A, D = s * s - 4 * r; w->n_iter++;
            if (D < 0 || !((SQ64 >> (D & 63)) & 1)) continue;
            int64_t q = isqrt64(D); if (q * q != D) continue;
            int64_t x = (s - q) / 2, y = s - x; if (x < LO || y >= HI || x * y != r) continue;
            Sol so; so.sig = l[a].sig; so.J = l[a].J; so.SJ = l[a].SJ; so.br = l[a].br; so.N[0] = (uint32_t)x; so.N[1] = (uint32_t)y; so.N[2] = 0;
            if (sys == 'X') { PUSH(w->sx, w->nsx, w->capsx, Sol); w->sx[w->nsx++] = so; } else { PUSH(w->sy, w->nsy, w->capsy, Sol); w->sy[w->nsy++] = so; }
        }
        return;
    }
    long double Bf = ldexpl(B, 3 * LAT);                  // Na * x * y
    double Bd = (double)Bf;                                // < 2^54: the quotient below is good to about 3e-5, covered by the slack in wd
    int64_t Amin = l[0].A, Amax = l[nl - 1].A;
    int64_t nlo = LO, nhi = HI - 1;
    if (Amax / 3 < nhi) nhi = Amax / 3;
    { int64_t t = (int64_t)floor(cbrt((double)Bf) * (1.0 + 1e-9)) + 1; if (t < nhi) nhi = t; }
    { int64_t t = Amin - 2 * (HI - 1); if (t > nlo) nlo = t; }
    { int64_t t = (int64_t)floor((double)(Bf / ((long double)(HI - 1) * (long double)(HI - 1))) * (1.0 - 1e-9)); if (t > nlo) nlo = t; }
    for (int64_t Na = nlo; Na <= nhi; Na++) {
        w->n_iter++;
        double preal = Bd / (double)Na, wd = preal * tol + 1e-4;
        if (wd < 0.49 && fabs(rint(preal) - preal) > wd) continue;
        int64_t rlo = (int64_t)ceil(preal - wd), rhi = (int64_t)floor(preal + wd);
        for (int64_t r = rlo; r <= rhi; r++) {
            if (r < Na * Na) continue;
            int64_t smin = (int64_t)ceil(2.0 * sqrt((double)r) - 1e-6);
            int64_t xmin = (r + HI - 2) / (HI - 1); if (xmin < Na) xmin = Na;         // y = r / x < HI
            int64_t smax = xmin + r / xmin + 1;
            if (smax + Na < Amin || smin + Na > Amax) continue;
            int a = 0, b = nl; while (a < b) { int m = (a + b) >> 1; if (l[m].A < smin + Na) a = m + 1; else b = m; }
            for (; a < nl && l[a].A <= smax + Na; a++) {
                int64_t s = l[a].A - Na, D = s * s - 4 * r;
                if (D < 0 || !((SQ64 >> (D & 63)) & 1)) continue;
                int64_t q = isqrt64(D); if (q * q != D) continue;
                int64_t x = (s - q) / 2, y = s - x; if (x < Na || y >= HI || x * y != r) continue;
                Sol so; so.sig = l[a].sig; so.J = l[a].J; so.SJ = l[a].SJ; so.br = l[a].br; so.N[0] = (uint32_t)Na; so.N[1] = (uint32_t)x; so.N[2] = (uint32_t)y;
                if (sys == 'X') { PUSH(w->sx, w->nsx, w->capsx, Sol); w->sx[w->nsx++] = so; } else { PUSH(w->sy, w->nsy, w->capsy, Sol); w->sy[w->nsy++] = so; }
            }
        }
    }
}

// complete key schedule of the chunks c[0..nc-1] (len of the last chunk irrelevant here: values already zero-padded)
static int full_check(const uint32_t *c, int nc, double *k, uint32_t *So, uint32_t *To, int *mapo) {
    double P[4] = { 1, 1, 1, 1 }; uint32_t S = 0, T = 0; int map = -9;
    for (int i = 0; i < nc; i++) { CV v = cvals(c[i]); for (int j = 0; j < 4; j++) P[j] *= 1.0 + v.f[j]; S += v.s; T += v.t; map = v.map; }
    for (int j = 0; j < 4; j++) { if (P[j] > 1.0) { double ip = floor(P[j]); k[j] = P[j] - ip; S += (uint32_t)ip; } else k[j] = P[j]; }
    *So = S; *To = T; *mapo = map;
    if (S != S_t) return 0;
    uint32_t t16 = T & 0xFFFF; if (t16 < Tlo || t16 > Thi) return 0;
    int ok = 0; for (int b = 0; b < nbr; b++) if (fabs(k[0] - c1[b]) <= h1[b] * WIDEN + WADD) ok = 1;
    if (!ok) return 0;
    if (fabs(k[2] - c3) > h3 * WIDEN + WADD) return 0;
    if (has2 && (fabs(k[1] - k788t) > h2t * WIDEN + WADD || fabs(k[3] - k7b8t) > h4t * WIDEN + WADD)) return 0;
    return 1;
}

static void emit_key(WS *w, const HYP *h, int64_t asg, const uint32_t *c) {
    double k[4]; uint32_t S, T; int map;
    full_check(c, h->nch, k, &S, &T, &map);
    char line[512]; int o = snprintf(line, sizeof line, "KEY %s %lld ", h->id, (long long)asg);
    for (int i = 0; i < h->nch; i++) for (int b = 0; b < h->ch[i].len; b++) o += snprintf(line + o, sizeof line - o, "%02x", (c[i] >> (8 * b)) & 0xFF);
    o += snprintf(line + o, sizeof line - o, " %08x %04x %.17g %.17g %.17g %.17g\n", S, T & 0xFFFF, k[0], k[1], k[2], k[3]);
    outs(w, line, o); w->n_keys++;
}

// a set of values for the U chunks: check it and print every ordering that gives the wanted map
static void try_set(WS *w, const HYP *h, int64_t asg, uint32_t *cfull, const uint32_t *u, int nU) {
    int slot[3], ns = 0; for (int i = 0; i < h->nch; i++) if (h->ch[i].type == 'U') slot[ns++] = i;
    for (int i = 0; i < nU; i++) cfull[slot[i]] = u[i];
    double k[4]; uint32_t S, T; int map; w->n_cand++;
    if (!full_check(cfull, h->nch, k, &S, &T, &map)) return;
    static const int P3[6][3] = { {0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0} };
    static const int P2[2][3] = { {0,1,0},{1,0,0} }, P1[1][3] = { {0,0,0} };
    const int (*PM)[3] = nU == 3 ? P3 : nU == 2 ? P2 : P1;
    int np = nU == 3 ? 6 : nU == 2 ? 2 : 1;
    for (int p = 0; p < np; p++) {
        uint32_t v[3]; for (int i = 0; i < nU; i++) v[i] = u[PM[p][i]];
        int dup = 0; for (int q = 0; q < p && !dup; q++) { int same = 1; for (int i = 0; i < nU; i++) if (u[PM[q][i]] != v[i]) same = 0; dup = same; }
        if (dup) continue;
        for (int i = 0; i < nU; i++) cfull[slot[i]] = v[i];
        if (cvals(cfull[h->nch - 1]).map != wantmap) continue;
        emit_key(w, h, asg, cfull);
    }
}

// join of the X and Y solutions of one sigma: chunks from pairs of lattice values, sum of the chunks = sigma
static void join_sigma(WS *w, const HYP *h, int64_t asg, uint32_t *cfull, int nU, int64_t sigma, const Sol *X, int nx, const Sol *Y, int ny) {
    int nxv = 0, nyv = 0;
    if ((size_t)(nx + ny) * 3 > w->caphead) { w->caphead = (size_t)(nx + ny) * 6 + 64; w->xv = realloc(w->xv, w->caphead * sizeof(uint32_t)); w->yv = realloc(w->yv, w->caphead * sizeof(uint32_t)); }
    for (int i = 0; i < nx; i++) for (int k = 0; k < nU; k++) w->xv[nxv++] = X[i].N[k];
    for (int i = 0; i < ny; i++) for (int k = 0; k < nU; k++) w->yv[nyv++] = Y[i].N[k];
    qsort(w->xv, nxv, sizeof(uint32_t), cmp_u32); qsort(w->yv, nyv, sizeof(uint32_t), cmp_u32);
    { int m = 0; for (int i = 0; i < nxv; i++) if (!m || w->xv[i] != w->xv[m - 1]) w->xv[m++] = w->xv[i]; nxv = m; }
    { int m = 0; for (int i = 0; i < nyv; i++) if (!m || w->yv[i] != w->yv[m - 1]) w->yv[m++] = w->yv[i]; nyv = m; }
    int *head = calloc((size_t)nxv * nyv, sizeof(int)); size_t pn = 1;
    if (!w->cappool) { w->cappool = 1 << 16; w->pool = malloc(w->cappool * sizeof(PE)); }
    for (int iy = 0; iy < nyv; iy++) {
        uint32_t yv = w->yv[iy];
        for (int ix = 0; ix < nxv; ix++) {
            uint32_t xv = w->xv[ix];
            for (int q = 0; q < ncp; q++) {
                const CP *cp = &cps[q];
                if ((yv & ((1u << (LAT - cp->ell1)) - 1)) || (xv & ((1u << (LAT - cp->ell3)) - 1))) continue;
                int64_t wd1 = 1LL << (54 - cp->ell1), a1 = (int64_t)((uint64_t)(yv - (1u << LAT)) >> (LAT - cp->ell1)) << (54 - cp->ell1);
                int64_t wd3 = 1LL << (54 - cp->ell3), a3 = (int64_t)((uint64_t)(xv - (1u << LAT)) >> (LAT - cp->ell3)) << (54 - cp->ell3);
                uint32_t found[32];
                int nf = cp_enum(cp, a1 - (wd1 >> 8), a1 + wd1 - 1, a3 - (wd3 >> 8), a3 + wd3 - 1, found, 32, &w->n_fib);
                for (int e = 0; e < nf; e++) {
                    uint32_t c = found[e]; int l1, l3;
                    if (latv(c, K_PI, &l1) != yv || l1 != cp->ell1 || latv(c, K_E, &l3) != xv || l3 != cp->ell3) continue;
                    if (pn == w->cappool) { w->cappool *= 2; w->pool = realloc(w->pool, w->cappool * sizeof(PE)); }
                    size_t hh = (size_t)iy * nxv + ix; w->pool[pn].c = c; w->pool[pn].next = head[hh]; head[hh] = (int)pn; pn++;
                }
            }
        }
    }
    static const int P3[6][3] = { {0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0} };
    static const int P2[2][3] = { {0,1,0},{1,0,0} };
    int np = nU == 3 ? 6 : 2;
    for (int i = 0; i < nx; i++) {
        int ix[3]; for (int k = 0; k < nU; k++) { int lo2 = 0, hi2 = nxv; while (lo2 < hi2) { int md = (lo2 + hi2) >> 1; if (w->xv[md] < X[i].N[k]) lo2 = md + 1; else hi2 = md; } ix[k] = lo2; }
        for (int j = 0; j < ny; j++) {
            if (Y[j].SJ - Y[j].J - X[i].J < 2) continue;
            int iyy[3]; for (int k = 0; k < nU; k++) { int lo2 = 0, hi2 = nyv; while (lo2 < hi2) { int md = (lo2 + hi2) >> 1; if (w->yv[md] < Y[j].N[k]) lo2 = md + 1; else hi2 = md; } iyy[k] = lo2; }
            int hd[3][3];
            for (int a = 0; a < nU; a++) for (int b = 0; b < nU; b++) hd[a][b] = head[(size_t)iyy[b] * nxv + ix[a]];
            for (int p = 0; p < np; p++) {
                const int *pm = nU == 3 ? P3[p] : P2[p];
                int ok = 1; for (int a = 0; a < nU; a++) if (!hd[a][pm[a]]) { ok = 0; break; }
                if (!ok) continue;
                uint32_t cc[3];
                for (int e0 = hd[0][pm[0]]; e0; e0 = w->pool[e0].next) for (int e1 = hd[1][pm[1]]; e1; e1 = w->pool[e1].next) {
                    cc[0] = w->pool[e0].c; cc[1] = w->pool[e1].c;
                    if (nU == 2) { if ((int64_t)cc[0] + cc[1] == sigma) try_set(w, h, asg, cfull, cc, 2); }
                    else for (int e2 = hd[2][pm[2]]; e2; e2 = w->pool[e2].next) { cc[2] = w->pool[e2].c; if ((int64_t)cc[0] + cc[1] + cc[2] == sigma) try_set(w, h, asg, cfull, cc, 3); }
                }
            }
        }
    }
    free(head);
}

// f1 side for two unknown chunks, without candidate lists: N1a * N1b = r must be an integer within the tolerance of the target, which happens
// for only a few (branch, J); only for those the sums (sigma, j1, d) are tried.
typedef struct { int64_t sig, SJ[2], b1[2]; int Ya[2][2], Yb[2][2], dmax, ok[2]; } Y2;
typedef struct { int br, J; int64_t rlo, rhi; } H2;
// the (branch, J) for which N1a * N1b can be an integer within the tolerance; returns the number found (at most maxh; more than that is not expected)
static int hits2(WS *w, const double *D, H2 *hh, int maxh) {
    int nh = 0;
    for (int br = 0; br < nbr; br++) {
        int Jlo = (int)ceil(D[0] * (1.0 - 2e-3) - c1[br] - 1e-6), Jhi = (int)floor(D[0] * (4.0 + 2e-3) - c1[br] + 1e-6); if (Jlo < 1) Jlo = 1;
        for (int J = Jlo; J <= Jhi; J++) {
            double preal = ldexp(((double)J + c1[br]) / D[0], 2 * LAT), tol = (h1[br] * WIDEN + WADD) / (J + c1[br]) + 5e-15, wd = preal * tol + 1e-4;
            double rr = rint(preal); w->n_iter++;
            if (wd < 0.49 && fabs(rr - preal) > wd) continue;
            int64_t rlo = (int64_t)ceil(preal - wd), rhi = (int64_t)floor(preal + wd);
            if (rhi < rlo) continue;
            if (nh == maxh) { fprintf(stderr, "hits2: too many\n"); exit(4); }
            hh[nh].br = br; hh[nh].J = J; hh[nh].rlo = rlo; hh[nh].rhi = rhi; nh++;
        }
    }
    return nh;
}
static void yside2(WS *w, const double *D, uint32_t Sp, uint32_t TK, const H2 *hh, int nh) {
    const int nU = 2; size_t nd = 0;
    if (w->nsg > w->capy2) { w->capy2 = w->nsg * 2 + 64; w->y2 = realloc(w->y2, w->capy2 * sizeof(Y2)); }
    Y2 *y2 = (Y2 *)w->y2;
    for (size_t i = 0; i < w->nsg; i++) {
        const SigC *s = &w->sg[i];
        if (s->j3 > 0) {                                  // the same sigma listed with a smaller j3: already taken
            int64_t v = (int64_t)(((u128)ME * (u128)s->sig) >> 54); int dup = 0;
            for (int jj = 0; jj < s->j3 && !dup; jj++) { int64_t t16 = (v + 2 * nU + (int64_t)TK - jj) & 0xFFFF; if (t16 >= (int64_t)Tlo && t16 <= (int64_t)Thi) dup = 1; }
            if (dup) continue;
        }
        Y2 *y = &y2[nd]; y->sig = s->sig; double lo, hi;
        int64_t Sw = (int64_t)Sp + ((int64_t)s->w << 32);
        u128 p1 = (u128)MP * (u128)s->sig; int64_t fl1 = (int64_t)(p1 >> 54); uint64_t fr1 = (uint64_t)(p1 & M54);
        y->dmax = nU << (LAT - 14); { double t = K_PI * (double)s->sig / ldexp(1.0, 44 - LAT); if (t < y->dmax) y->dmax = (int)t; }
        int any = 0;
        for (int j1 = 0; j1 < nU; j1++) {
            int64_t I1 = fl1 - j1, SJ = Sw - 3 * nU - I1; y->ok[j1] = !(SJ < 4 || SJ > sjmax); if (!y->ok[j1]) continue;
            y->SJ[j1] = SJ; double F1 = j1 + (double)fr1 / TWO54; bounds(F1, nU, &lo, &hi);
            y->b1[j1] = ((int64_t)nU << LAT) + ((int64_t)j1 << LAT) + (int64_t)(fr1 >> (54 - LAT));
            for (int br = 0; br < nbr; br++) {
                int Ya = (int)ceil((lo - 2e-3) * D[0] - c1[br] - 1e-6), Yb = (int)floor((hi + 2e-3) * D[0] - c1[br] + 1e-6); if (Ya < 1) Ya = 1; if (Yb > SJ - 3) Yb = (int)(SJ - 3);
                y->Ya[j1][br] = Ya; y->Yb[j1][br] = Yb; if (Yb >= Ya) any = 1;
            }
        }
        if (any) nd++;
    }
    if (!nd) return;
    for (int k = 0; k < nh; k++) {
        int br = hh[k].br, J = hh[k].J;
        for (size_t i = 0; i < nd; i++) for (int j1 = 0; j1 < nU; j1++) {
            const Y2 *y = &y2[i]; if (!y->ok[j1] || J < y->Ya[j1][br] || J > y->Yb[j1][br]) continue;
            for (int d = 0; d <= y->dmax; d++) for (int64_t r = hh[k].rlo; r <= hh[k].rhi; r++) {
                int64_t sm = y->b1[j1] - d, DD = sm * sm - 4 * r;
                if (DD < 0 || !((SQ64 >> (DD & 63)) & 1)) continue;
                int64_t q = isqrt64(DD); if (q * q != DD) continue;
                int64_t x = (sm - q) / 2, yy = sm - x; if (x < LO || yy >= HI || x * yy != r) continue;
                Sol so; so.sig = y->sig; so.J = J; so.SJ = (int)y->SJ[j1]; so.br = br; so.N[0] = (uint32_t)x; so.N[1] = (uint32_t)yy; so.N[2] = 0;
                PUSH(w->sy, w->nsy, w->capsy, Sol); w->sy[w->nsy++] = so;
            }
        }
    }
}

// one assignment of the M bytes
static void solve_assignment(WS *w, const HYP *h, int64_t asg) {
    uint32_t cfull[32]; double D[4] = { 1, 1, 1, 1 }; uint32_t SK = 0, TK = 0; int lastmap = -9;
    int64_t v = asg;
    for (int i = 0; i < h->nch; i++) {
        const CH *c = &h->ch[i];
        if (c->type == 'U') { cfull[i] = 0; continue; }
        uint32_t val = c->base;
        if (c->type == 'M') for (int b = 0; b < c->len; b++) if (c->mask & (1u << b)) { val = (val & ~(0xFFu << (8 * b))) | ((uint32_t)(blo + v % nb) << (8 * b)); v /= nb; }
        cfull[i] = val; CV cv = cvals(val);
        for (int j = 0; j < 4; j++) D[j] *= 1.0 + cv.f[j];
        SK += cv.s; TK += cv.t; if (i == h->nch - 1) lastmap = cv.map;
    }
    w->n_assign++;
    if (!h->lastU && lastmap != wantmap) return;
    int nU = h->nU;
    if (nU == 0) {
        double k[4]; uint32_t S, T; int map; w->n_feas++; w->n_cand++;
        if (full_check(cfull, h->nch, k, &S, &T, &map)) emit_key(w, h, asg, cfull);
        return;
    }
    uint32_t Sp = S_t - SK;
    H2 hh[256]; int nh = 0;
    if (nU == 2 && !XFIRST) { nh = hits2(w, D, hh, 256); if (!nh) { w->n_feas++; return; } }     // no (branch, J) can give an integer product: nothing to solve
    sigma_cands(w, Sp, TK, nU);
    if (!w->nsg) return;
    w->n_feas++;
    if (nU == 1) {
        int64_t last = -1;
        for (size_t i = 0; i < w->nsg; i++) {
            int64_t sg = w->sg[i].sig; if (sg == last || sg > 0xFFFFFFFFLL) continue; last = sg;
            double a3 = clr((double)(uint32_t)sg * K_E), p3 = D[2] * (1.0 + (a3 - floor(a3))), f = p3 - floor(p3);
            if (fabs(f - c3) > h3 * WIDEN + WADD + 1e-12) continue;
            uint32_t u[1] = { (uint32_t)sg }; try_set(w, h, asg, cfull, u, 1);
        }
        return;
    }
    if (XFIRST) {
    // f3 side (k7b0) first: one branch only, and each sigma comes with its own j3
    w->ncx = w->ncy = 0; w->nsx = w->nsy = 0;
    for (size_t i = 0; i < w->nsg; i++) {
        const SigC *s = &w->sg[i]; double lo, hi;
        double F = s->j3 + (double)s->fr3 / TWO54; bounds(F, nU, &lo, &hi);
        int64_t base = ((int64_t)nU << LAT) + ((int64_t)s->j3 << LAT) + (int64_t)(s->fr3 >> (54 - LAT));
        int dmax = nU << (LAT - 14); { double t = K_E * (double)s->sig / ldexp(1.0, 44 - LAT); if (t < dmax) dmax = (int)t; }
        int Ja = (int)ceil((lo - 2e-3) * D[2] - c3 - 1e-6), Jb = (int)floor((hi + 2e-3) * D[2] + 1e-6 - c3); if (Ja < 1) Ja = 1;
        for (int J = Ja; J <= Jb; J++) for (int d = 0; d <= dmax; d++) { PUSH(w->cx, w->ncx, w->capx, Cand); Cand *c = &w->cx[w->ncx++]; c->A = base - d; c->sig = s->sig; c->J = J; c->SJ = s->w; c->br = 0; }
    }
    if (!w->ncx) return;
    qsort(w->cx, w->ncx, sizeof(Cand), cmp_cand);
    for (size_t a = 0; a < w->ncx; ) {
        size_t b = a; while (b < w->ncx && w->cx[b].J == w->cx[a].J) b++;
        int J = w->cx[a].J;
        enum_run(w, 'X', nU, w->cx + a, (int)(b - a), ((long double)J + (long double)c3) / (long double)D[2], (h3 * WIDEN + WADD) / (J + c3) + 5e-15);
        a = b;
    }
    w->n_xsol += w->nsx;
    if (!w->nsx) return;
    qsort(w->sx, w->nsx, sizeof(Sol), cmp_sol);
    // f1 side (k780, both branches), only for the sigma values that have a solution on the f3 side
    {
        int64_t last = -1;
        for (size_t i = 0; i < w->nsx; i++) {
            int64_t sig = w->sx[i].sig; if (sig == last) continue; last = sig;
            double lo, hi; int wr = w->sx[i].SJ;
            int64_t Sw = (int64_t)Sp + ((int64_t)wr << 32);
            u128 p1 = (u128)MP * (u128)sig; int64_t fl1 = (int64_t)(p1 >> 54); uint64_t fr1 = (uint64_t)(p1 & M54);
            int dmax1 = nU << (LAT - 14); { double t = K_PI * (double)sig / ldexp(1.0, 44 - LAT); if (t < dmax1) dmax1 = (int)t; }
            for (int j1 = 0; j1 < nU; j1++) {
                int64_t I1 = fl1 - j1, SJ = Sw - 3 * nU - I1; if (SJ < 4 || SJ > sjmax) continue;
                double F1 = j1 + (double)fr1 / TWO54; bounds(F1, nU, &lo, &hi);
                int64_t b1 = ((int64_t)nU << LAT) + ((int64_t)j1 << LAT) + (int64_t)(fr1 >> (54 - LAT));
                for (int br = 0; br < nbr; br++) {
                    int Ya = (int)ceil((lo - 2e-3) * D[0] - c1[br] - 1e-6), Yb = (int)floor((hi + 2e-3) * D[0] - c1[br] + 1e-6); if (Ya < 1) Ya = 1; if (Yb > SJ - 3) Yb = (int)(SJ - 3);
                    for (int J = Ya; J <= Yb; J++) for (int d = 0; d <= dmax1; d++) { PUSH(w->cy, w->ncy, w->capy, Cand); Cand *c = &w->cy[w->ncy++]; c->A = b1 - d; c->sig = sig; c->J = J; c->SJ = (int)SJ; c->br = br; }
                }
            }
        }
    }
    if (!w->ncy) return;
    qsort(w->cy, w->ncy, sizeof(Cand), cmp_cand);
    for (size_t a = 0; a < w->ncy; ) {
        size_t b = a; while (b < w->ncy && w->cy[b].J == w->cy[a].J && w->cy[b].br == w->cy[a].br) b++;
        int J = w->cy[a].J, br = w->cy[a].br;
        enum_run(w, 'Y', nU, w->cy + a, (int)(b - a), ((long double)J + (long double)c1[br]) / (long double)D[0], (h1[br] * WIDEN + WADD) / (J + c1[br]) + 5e-15);
        a = b;
    }
    w->n_ysol += w->nsy;
    if (!w->nsy) return;
    qsort(w->sy, w->nsy, sizeof(Sol), cmp_sol);
    } else {
    // f1 side (k780) first: it does not depend on j3
    w->ncx = w->ncy = 0; w->nsx = w->nsy = 0;
    if (nU == 2) yside2(w, D, Sp, TK, hh, nh);
    else {
    {
        if (w->nsg > w->capds) { w->capds = w->nsg * 2 + 64; w->ds = realloc(w->ds, w->capds * sizeof(SigC)); }
        memcpy(w->ds, w->sg, w->nsg * sizeof(SigC)); qsort(w->ds, w->nsg, sizeof(SigC), cmp_sigc);
        size_t nd = 0; for (size_t i = 0; i < w->nsg; i++) if (!nd || w->ds[i].sig != w->ds[nd - 1].sig) w->ds[nd++] = w->ds[i];
        for (size_t i = 0; i < nd; i++) {
            const SigC *s = &w->ds[i]; double lo, hi;
            int64_t Sw = (int64_t)Sp + ((int64_t)s->w << 32);
            u128 p1 = (u128)MP * (u128)s->sig; int64_t fl1 = (int64_t)(p1 >> 54); uint64_t fr1 = (uint64_t)(p1 & M54);
            int dmax1 = nU << (LAT - 14); { double t = K_PI * (double)s->sig / ldexp(1.0, 44 - LAT); if (t < dmax1) dmax1 = (int)t; }
            for (int j1 = 0; j1 < nU; j1++) {
                int64_t I1 = fl1 - j1, SJ = Sw - 3 * nU - I1; if (SJ < 4 || SJ > sjmax) continue;
                double F1 = j1 + (double)fr1 / TWO54; bounds(F1, nU, &lo, &hi);
                int64_t b1 = ((int64_t)nU << LAT) + ((int64_t)j1 << LAT) + (int64_t)(fr1 >> (54 - LAT));
                for (int br = 0; br < nbr; br++) {
                    int Ya = (int)ceil((lo - 2e-3) * D[0] - c1[br] - 1e-6), Yb = (int)floor((hi + 2e-3) * D[0] - c1[br] + 1e-6); if (Ya < 1) Ya = 1; if (Yb > SJ - 3) Yb = (int)(SJ - 3);
                    for (int J = Ya; J <= Yb; J++) for (int d = 0; d <= dmax1; d++) { PUSH(w->cy, w->ncy, w->capy, Cand); Cand *c = &w->cy[w->ncy++]; c->A = b1 - d; c->sig = s->sig; c->J = J; c->SJ = (int)SJ; c->br = br; }
                }
            }
        }
    }
    if (!w->ncy) return;
    qsort(w->cy, w->ncy, sizeof(Cand), cmp_cand);
    for (size_t a = 0; a < w->ncy; ) {
        size_t b = a; while (b < w->ncy && w->cy[b].J == w->cy[a].J && w->cy[b].br == w->cy[a].br) b++;
        int J = w->cy[a].J, br = w->cy[a].br;
        enum_run(w, 'Y', nU, w->cy + a, (int)(b - a), ((long double)J + (long double)c1[br]) / (long double)D[0], (h1[br] * WIDEN + WADD) / (J + c1[br]) + 5e-15);
        a = b;
    }
    }
    w->n_ysol += w->nsy;
    if (!w->nsy) return;
    qsort(w->sy, w->nsy, sizeof(Sol), cmp_sol);
    // f3 side (k7b0), only for the sigma values that have a solution on the f1 side
    for (size_t i = 0; i < w->nsg; i++) {
        const SigC *s = &w->sg[i]; double lo, hi;
        { size_t lo2 = 0, hi2 = w->nsy; while (lo2 < hi2) { size_t md = (lo2 + hi2) >> 1; if (w->sy[md].sig < s->sig) lo2 = md + 1; else hi2 = md; }
          if (lo2 >= w->nsy || w->sy[lo2].sig != s->sig) continue; }
        double F = s->j3 + (double)s->fr3 / TWO54; bounds(F, nU, &lo, &hi);
        int64_t base = ((int64_t)nU << LAT) + ((int64_t)s->j3 << LAT) + (int64_t)(s->fr3 >> (54 - LAT));
        int dmax = nU << (LAT - 14); { double t = K_E * (double)s->sig / ldexp(1.0, 44 - LAT); if (t < dmax) dmax = (int)t; }
        int Ja = (int)ceil((lo - 2e-3) * D[2] - c3 - 1e-6), Jb = (int)floor((hi + 2e-3) * D[2] - c3 + 1e-6); if (Ja < 1) Ja = 1;
        for (int J = Ja; J <= Jb; J++) for (int d = 0; d <= dmax; d++) { PUSH(w->cx, w->ncx, w->capx, Cand); Cand *c = &w->cx[w->ncx++]; c->A = base - d; c->sig = s->sig; c->J = J; c->SJ = 0; c->br = 0; }
    }
    if (!w->ncx) return;
    qsort(w->cx, w->ncx, sizeof(Cand), cmp_cand);
    for (size_t a = 0; a < w->ncx; ) {
        size_t b = a; while (b < w->ncx && w->cx[b].J == w->cx[a].J) b++;
        int J = w->cx[a].J;
        enum_run(w, 'X', nU, w->cx + a, (int)(b - a), ((long double)J + (long double)c3) / (long double)D[2], (h3 * WIDEN + WADD) / (J + c3) + 5e-15);
        a = b;
    }
    w->n_xsol += w->nsx;
    if (!w->nsx) return;
    qsort(w->sx, w->nsx, sizeof(Sol), cmp_sol);
    }
    size_t i = 0, j = 0;
    while (i < w->nsx && j < w->nsy) {
        if (w->sx[i].sig < w->sy[j].sig) { i++; continue; }
        if (w->sx[i].sig > w->sy[j].sig) { j++; continue; }
        int64_t sg = w->sx[i].sig; size_t i2 = i, j2 = j;
        while (i2 < w->nsx && w->sx[i2].sig == sg) i2++;
        while (j2 < w->nsy && w->sy[j2].sig == sg) j2++;
        w->n_grp++;
        join_sigma(w, h, asg, cfull, nU, sg, w->sx + i, (int)(i2 - i), w->sy + j, (int)(j2 - j));
        i = i2; j = j2;
    }
}

int main(int argc, char **argv) {
    K_PI = h2d(0x3fd45f306dda35a5ULL); K_E = h2d(0x3fd78b563618cf81ULL); C01 = h2d(0x3fb999999999999aULL);
    MP = (0x3fd45f306dda35a5ULL & ((1ULL << 52) - 1)) | (1ULL << 52); ME = (0x3fd78b563618cf81ULL & ((1ULL << 52) - 1)) | (1ULL << 52);
    if (argc < 2) { fprintf(stderr, "usage: vres hypotheses.txt\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    char tok[64]; int capH = 0; LAT = 16; blo = 0x20; bhi = 0xFC; sjmax = 4092;
    while (fscanf(f, "%63s", tok) == 1) {
        if (!strcmp(tok, "T")) {
            unsigned s_, tl, th; if (fscanf(f, "%x %x %x %d %lf %lf %d", &s_, &tl, &th, &wantmap, &c3, &h3, &nbr) != 7) return 2;
            S_t = s_; Tlo = tl; Thi = th;
            for (int b = 0; b < nbr; b++) if (fscanf(f, "%lf %lf", &c1[b], &h1[b]) != 2) return 2;
            if (fscanf(f, "%d", &has2) != 1) return 2;
            if (has2 && fscanf(f, "%lf %lf %lf %lf", &k788t, &h2t, &k7b8t, &h4t) != 4) return 2;
        } else if (!strcmp(tok, "O")) {
            long long sj; if (fscanf(f, "%d %i %i %lld", &LAT, &blo, &bhi, &sj) != 4) return 2; sjmax = sj;
        } else if (!strcmp(tok, "H")) {
            if (nhyp == capH) { capH = capH ? capH * 2 : 256; hyps = realloc(hyps, capH * sizeof(HYP)); }
            HYP *h = &hyps[nhyp]; memset(h, 0, sizeof *h);
            if (fscanf(f, "%63s %d", h->id, &h->nch) != 2 || h->nch < 1 || h->nch > 32) return 2;
            h->ch = calloc(h->nch, sizeof(CH));
            for (int i = 0; i < h->nch; i++) {
                char t[8]; if (fscanf(f, "%7s", t) != 1) return 2; CH *c = &h->ch[i]; c->type = t[0]; c->len = 4;
                if (t[0] == 'K') { if (fscanf(f, "%x %d", &c->base, &c->len) != 2) return 2; }
                else if (t[0] == 'M') { if (fscanf(f, "%x %x %d", &c->base, &c->mask, &c->len) != 3) return 2; h->nM++; h->E += __builtin_popcount(c->mask); }
                else if (t[0] == 'U') h->nU++;
                else return 2;
            }
            h->lastU = h->ch[h->nch - 1].type == 'U';
            nhyp++;
        } else { fprintf(stderr, "bad token %s\n", tok); return 2; }
    }
    fclose(f);
    nb = bhi - blo + 1; LO = 1LL << LAT; HI = 2LL << LAT; CMIN = LAT == 16 ? 0x33000000u : 0x20000000u;
    if (LAT != 16 && LAT != 17) { fprintf(stderr, "L must be 16 or 17\n"); return 2; }
    init_steps(); init_cp();
    if (getenv("VRES_XFIRST")) XFIRST = atoi(getenv("VRES_XFIRST"));
    // tasks: (hypothesis, block of assignments)
    typedef struct { int h; int64_t a0, a1; } TASK; TASK *tasks = NULL; size_t nt = 0, capt = 0;
    for (int i = 0; i < nhyp; i++) {
        HYP *h = &hyps[i]; if (h->nU > 3) { fprintf(stderr, "%s: more than 3 unknown full chunks\n", h->id); return 2; }
        h->nassign = 1; for (int e = 0; e < h->E; e++) h->nassign *= nb;
        int64_t blk = h->nU >= 2 ? 4 : 4096;
        for (int64_t a = 0; a < h->nassign; a += blk) { PUSH(tasks, nt, capt, TASK); tasks[nt].h = i; tasks[nt].a0 = a; tasks[nt].a1 = a + blk < h->nassign ? a + blk : h->nassign; nt++; }
    }
    uint64_t (*st)[9] = calloc(nhyp, sizeof *st); double t0 = omp_get_wtime();
    #pragma omp parallel
    {
        WS w; memset(&w, 0, sizeof w);
        #pragma omp for schedule(dynamic, 1)
        for (size_t t = 0; t < nt; t++) {
            const HYP *h = &hyps[tasks[t].h];
            uint64_t b[9] = { w.n_assign, w.n_feas, w.n_xsol, w.n_ysol, w.n_grp, w.n_cand, w.n_keys, w.n_iter, w.n_fib };
            for (int64_t a = tasks[t].a0; a < tasks[t].a1; a++) solve_assignment(&w, h, a);
            uint64_t e[9] = { w.n_assign, w.n_feas, w.n_xsol, w.n_ysol, w.n_grp, w.n_cand, w.n_keys, w.n_iter, w.n_fib };
            #pragma omp critical
            { for (int k = 0; k < 9; k++) st[tasks[t].h][k] += e[k] - b[k]; if (w.lout) { fwrite(w.out, 1, w.lout, stdout); fflush(stdout); w.lout = 0; } }
        }
    }
    for (int i = 0; i < nhyp; i++)
        printf("STAT %s nU %d E %d assignments %llu solved %llu ysol %llu xsol %llu groups %llu candidates %llu keys %llu iter %llu fibre %llu\n", hyps[i].id, hyps[i].nU, hyps[i].E,
               (unsigned long long)st[i][0], (unsigned long long)st[i][1], (unsigned long long)st[i][3], (unsigned long long)st[i][2], (unsigned long long)st[i][4],
               (unsigned long long)st[i][5], (unsigned long long)st[i][6], (unsigned long long)st[i][7], (unsigned long long)st[i][8]);
    fprintf(stderr, "vres: %d hypotheses, %zu tasks, %.1f s\n", nhyp, nt, omp_get_wtime() - t0);
    return 0;
}
