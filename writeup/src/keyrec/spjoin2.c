// spjoin2: join of the f3 system (X) and the f1 system (Y) for keys of n full chunks with arbitrary bytes
// (every chunk >= CMIN, so that both lattice classes are <= 16; all values are given on the 2^-16 lattice).
// For each sigma candidate: the chunks having one of the Y values as f1 value are enumerated (three-gap stepping along the
// rotation by K_PI), kept when their f3 value occurs in an X solution, and every (X solution, Y solution, pairing) whose
// chunks sum to sigma is checked against the key-schedule targets.
// build: gcc -O2 -march=native -ffp-contract=off -fno-fast-math -fopenmp -o spjoin2 spjoin2.c -lm
// usage: spjoin2 n params.txt X.txt Y.txt   (X: tag J N..; Y: tag SJ J br N..; both sorted by tag)
// params.txt: S Tlo Thi k780 k788 k7b0 k7b8 tol nbr(1|2) ; ntags ; then "tag sigma" lines.  tol < 0: print every sum-consistent candidate
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define CMIN 0x33000000u
#define M54 ((1ULL << 54) - 1)
static double K_PI, K_E; static uint64_t MP, ME;
static double h2d(uint64_t h) { double d; memcpy(&d, &h, 8); return d; }
static inline double clr(double a) { uint64_t q; memcpy(&q, &a, 8); q &= ~0xFFULL; memcpy(&a, &q, 8); return a; }
// value on the 2^-16 lattice (65536..131071) and class; 0 if the class is finer than 16
static inline uint32_t lat16(uint32_t c, double K, int *ell) {
    double a = clr((double)c * K); uint64_t q; memcpy(&q, &a, 8);
    int e = (int)((q >> 52) & 0x7FF) - 1023; *ell = 44 - e;
    if (*ell > 16 || *ell < 1) return 0;
    double f = a - floor(a);
    return ((1u << *ell) + (uint32_t)(f * (double)(1u << *ell))) << (16 - *ell);
}
typedef struct { int tag, J; uint32_t N[4]; } XS;
typedef struct { int tag, SJ, J, br; uint32_t N[4]; } YS;
typedef struct { uint32_t c; int next; } PE;
static int n; static uint32_t S_t, Tlo, Thi; static double kt[4], tolk; static int nbr;
static int64_t *sigma_of; static int ntags;
static const uint64_t F1LO[3] = { 0xC90FDAA2ULL, 0x6487ED51ULL, CMIN }, F1HI[3] = { 0x100000000ULL, 0xC90FDAA2ULL, 0x6487ED51ULL };   // class 14, 15, 16
static uint64_t g_up[3], d_up[3], g_dn[3], d_dn[3], WW[3];

static void init_steps(void) {
    for (int k = 0; k < 3; k++) {
        int ell = 14 + k; uint64_t w = 1ULL << (54 - ell); WW[k] = w + (w >> 8);
        g_up[k] = g_dn[k] = 0; uint64_t t = 0;
        for (uint64_t g = 1; !(g_up[k] && g_dn[k]); g++) {
            t = (t + MP) & M54;
            if (!g_up[k] && t < WW[k]) { g_up[k] = g; d_up[k] = t; }
            if (!g_dn[k] && t > (M54 + 1) - WW[k]) { g_dn[k] = g; d_dn[k] = (M54 + 1) - t; }
        }
    }
}

static int verify(const uint32_t *c, FILE *out, int tag) {
    double P[4] = { 1, 1, 1, 1 }; uint32_t S = 0, T = 0;
    for (int i = 0; i < n; i++) {
        double a1 = clr((double)c[i] * K_PI), i1 = floor(a1), a2 = clr(i1 * K_PI), i2 = floor(a2);
        double a3 = clr((double)c[i] * K_E), i3 = floor(a3), a4 = clr(i3 * K_E), i4 = floor(a4);
        P[0] *= 1.0 + (a1 - i1); P[1] *= 1.0 + (a2 - i2); P[2] *= 1.0 + (a3 - i3); P[3] *= 1.0 + (a4 - i4);
        S += (uint32_t)(int64_t)(i3 / i2) + (uint32_t)(int64_t)i1; T += (uint32_t)(int64_t)(i1 / i4) + (uint32_t)(int64_t)i3;
    }
    double k[4];
    for (int j = 0; j < 4; j++) { double ip = floor(P[j]); if (P[j] > 1.0) { k[j] = P[j] - ip; S += (uint32_t)ip; } else k[j] = P[j]; }
    if (tolk >= 0) {
        if (S != S_t) return 0;
        uint32_t t16 = T & 0xFFFF; if (t16 < Tlo || t16 > Thi) return 0;
        int ok0 = fabs(k[0] - kt[0]) <= tolk || (nbr == 2 && fabs(k[0] - (1.0 - kt[0])) <= tolk);
        if (!ok0 || fabs(k[2] - kt[2]) > tolk) return 0;
        if (kt[1] >= 0 && (fabs(k[1] - kt[1]) > tolk || fabs(k[3] - kt[3]) > tolk)) return 0;
    }
    #pragma omp critical
    { fprintf(out, "%d", tag); for (int i = 0; i < n; i++) fprintf(out, " %u", c[i]); fprintf(out, " %08x %04x %.17g %.17g %.17g %.17g\n", S, T & 0xFFFF, k[0], k[1], k[2], k[3]); }
    return 1;
}


static int selftest(void) {
    // the three-gap stepping must list exactly the chunks found by a plain scan
    uint64_t seed = 88172645463325252ULL; int bad = 0; uint64_t tot = 0;
    for (int k = 0; k < 3; k++) for (int rep = 0; rep < 40; rep++) {
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        int ell = 14 + k; uint64_t idx = seed % (1ULL << ell), w = 1ULL << (54 - ell), a = idx << (54 - ell), m = w >> 8, W = WW[k];
        uint64_t span = F1HI[k] - F1LO[k]; uint64_t lo = F1LO[k] + (seed >> 20) % (span - (1ULL << 25)), hi = lo + (1ULL << 25);
        uint64_t cnt1 = 0, sum1 = 0, cnt2 = 0, sum2 = 0;
        uint64_t c = lo, o = ((c * MP) - a + m) & M54;
        while (o >= W && c < hi) { c++; o = (o + MP) & M54; }
        while (c < hi) { cnt1++; sum1 += c * 2654435761ULL;
            int up = (o + d_up[k] < W), dn = (o >= d_dn[k]);
            if (up && dn) { if (g_up[k] < g_dn[k]) { c += g_up[k]; o += d_up[k]; } else { c += g_dn[k]; o -= d_dn[k]; } }
            else if (up) { c += g_up[k]; o += d_up[k]; } else if (dn) { c += g_dn[k]; o -= d_dn[k]; } else { c += g_up[k] + g_dn[k]; o = o + d_up[k] - d_dn[k]; } }
        uint64_t th = ((lo * MP) - a + m) & M54;
        for (uint64_t x = lo; x < hi; x++) { if (th < W) { cnt2++; sum2 += x * 2654435761ULL; } th = (th + MP) & M54; }
        tot += cnt2; if (cnt1 != cnt2 || sum1 != sum2) { bad++; fprintf(stderr, "MISMATCH class %d idx %llu: %llu vs %llu\n", ell, (unsigned long long)idx, (unsigned long long)cnt1, (unsigned long long)cnt2); }
    }
    // the window must contain every chunk whose true lattice value is idx (check on random chunks)
    int miss = 0;
    for (int rep = 0; rep < 2000000; rep++) {
        seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17;
        uint32_t c = CMIN + (uint32_t)(seed % (0x100000000ULL - CMIN)); int l1; uint32_t v = lat16(c, K_PI, &l1); if (!v) continue;
        int k = l1 - 14; uint64_t idx = (uint64_t)(v - 65536) >> (16 - l1), w = 1ULL << (54 - l1), a = idx << (54 - l1), m = w >> 8;
        uint64_t o = (((uint64_t)c * MP) - a + m) & M54; if (o >= WW[k]) miss++;
        if (c < F1LO[k] - 8 || c >= F1HI[k] + 8) miss++;
    }
    fprintf(stderr, "selftest: 120 fibres, %llu chunks, stepping mismatches %d; window/class misses in 2,000,000 random chunks: %d\n", (unsigned long long)tot, bad, miss);
    fprintf(stderr, "steps: "); for (int k = 0; k < 3; k++) fprintf(stderr, "class %d g_up %llu g_dn %llu  ", 14 + k, (unsigned long long)g_up[k], (unsigned long long)g_dn[k]); fprintf(stderr, "\n");
    return bad || miss;
}

static const int PERM3[6][4] = { {0,1,2,0},{0,2,1,0},{1,0,2,0},{1,2,0,0},{2,0,1,0},{2,1,0,0} };
static int PERM4[24][4];

int main(int argc, char **argv) {
    K_PI = h2d(0x3fd45f306dda35a5ULL); K_E = h2d(0x3fd78b563618cf81ULL);
    MP = (0x3fd45f306dda35a5ULL & ((1ULL << 52) - 1)) | (1ULL << 52); ME = (0x3fd78b563618cf81ULL & ((1ULL << 52) - 1)) | (1ULL << 52);   // K = M / 2^54
    if (argc == 2 && !strcmp(argv[1], "selftest")) { init_steps(); return selftest(); }
    if (argc < 5) return 1;
    n = atoi(argv[1]);
    { int k = 0; for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++) for (int c = 0; c < 4; c++) for (int d = 0; d < 4; d++)
        if (a != b && a != c && a != d && b != c && b != d && c != d) { PERM4[k][0] = a; PERM4[k][1] = b; PERM4[k][2] = c; PERM4[k][3] = d; k++; } }
    FILE *f = fopen(argv[2], "r"); unsigned s_, tl, th;
    if (fscanf(f, "%x %x %x %lf %lf %lf %lf %lf %d %d", &s_, &tl, &th, &kt[0], &kt[1], &kt[2], &kt[3], &tolk, &nbr, &ntags) != 10) return 1;
    S_t = s_; Tlo = tl; Thi = th; sigma_of = calloc(ntags, sizeof(int64_t));
    for (int i = 0; i < ntags; i++) { int t; long long sg; if (fscanf(f, "%d %lld", &t, &sg) != 2) return 1; sigma_of[t] = sg; }
    fclose(f);
    size_t nx = 0, ny = 0, cx = 1 << 20, cy = 1 << 20; XS *X = malloc(cx * sizeof(XS)); YS *Y = malloc(cy * sizeof(YS));
    f = fopen(argv[3], "r");
    for (;;) { XS s; memset(&s, 0, sizeof s); if (fscanf(f, "%d %d", &s.tag, &s.J) != 2) break; for (int k = 0; k < n; k++) if (fscanf(f, "%u", &s.N[k]) != 1) return 1;
        if (nx == cx) { cx *= 2; X = realloc(X, cx * sizeof(XS)); } X[nx++] = s; }
    fclose(f); f = fopen(argv[4], "r");
    for (;;) { YS s; memset(&s, 0, sizeof s); if (fscanf(f, "%d %d %d %d", &s.tag, &s.SJ, &s.J, &s.br) != 4) break; for (int k = 0; k < n; k++) if (fscanf(f, "%u", &s.N[k]) != 1) return 1;
        if (ny == cy) { cy *= 2; Y = realloc(Y, cy * sizeof(YS)); } Y[ny++] = s; }
    fclose(f);
    init_steps();
    // tag ranges
    size_t *xs0 = calloc(ntags + 1, sizeof(size_t)), *ys0 = calloc(ntags + 1, sizeof(size_t));
    for (size_t i = 0; i < nx; i++) xs0[X[i].tag + 1]++;
    for (size_t i = 0; i < ny; i++) ys0[Y[i].tag + 1]++;
    for (int t = 0; t < ntags; t++) { xs0[t + 1] += xs0[t]; ys0[t + 1] += ys0[t]; }
    uint64_t tot_fib = 0, tot_keep = 0, tot_pairs = 0, tot_sum = 0, tot_out = 0;
    #pragma omp parallel for schedule(dynamic, 1) reduction(+:tot_fib,tot_keep,tot_pairs,tot_sum,tot_out)
    for (int tag = 0; tag < ntags; tag++) {
        size_t x0 = xs0[tag], x1 = xs0[tag + 1], y0 = ys0[tag], y1 = ys0[tag + 1];
        if (x0 == x1 || y0 == y1) continue;
        int64_t sigma = sigma_of[tag];
        int *idxX = malloc(65536 * sizeof(int)), *idxY = malloc(65536 * sizeof(int)); memset(idxX, -1, 65536 * sizeof(int)); memset(idxY, -1, 65536 * sizeof(int));
        int nxv = 0, nyv = 0; uint32_t *yv = malloc(4 * (y1 - y0) * sizeof(uint32_t));
        for (size_t i = x0; i < x1; i++) for (int k = 0; k < n; k++) { uint32_t v = X[i].N[k] - 65536; if (idxX[v] < 0) idxX[v] = nxv++; }
        for (size_t i = y0; i < y1; i++) for (int k = 0; k < n; k++) { uint32_t v = Y[i].N[k] - 65536; if (idxY[v] < 0) { idxY[v] = nyv; yv[nyv++] = v + 65536; } }
        int *head = calloc((size_t)nxv * nyv, sizeof(int)); size_t pcap = 1 << 20, pn = 1; PE *pool = malloc(pcap * sizeof(PE));
        // fibres of the Y values
        for (int iy = 0; iy < nyv; iy++) {
            uint32_t y16 = yv[iy];
            for (int k = 0; k < 3; k++) {
                int ell = 14 + k; if (y16 & ((1u << (16 - ell)) - 1)) continue;
                uint64_t idx = (uint64_t)(y16 - 65536) >> (16 - ell), w = 1ULL << (54 - ell), a = (idx << (54 - ell)), m = w >> 8, W = WW[k];
                uint64_t lo = F1LO[k] > 8 + CMIN ? F1LO[k] - 8 : CMIN, hi = F1HI[k] + 8 > 0x100000000ULL ? 0x100000000ULL : F1HI[k] + 8;
                uint64_t c = lo, o = ((c * MP) - a + m) & M54;
                while (o >= W && c < hi) { c++; o = (o + MP) & M54; }
                while (c < hi) {
                    tot_fib++;
                    int l1, l3; uint32_t v1 = lat16((uint32_t)c, K_PI, &l1);
                    if (v1 == y16 && l1 == ell) {
                        uint32_t v3 = lat16((uint32_t)c, K_E, &l3);
                        if (v3 && idxX[v3 - 65536] >= 0) {
                            if (pn == pcap) { pcap *= 2; pool = realloc(pool, pcap * sizeof(PE)); }
                            size_t h = (size_t)iy * nxv + idxX[v3 - 65536]; pool[pn].c = (uint32_t)c; pool[pn].next = head[h]; head[h] = (int)pn; pn++; tot_keep++;
                        }
                    }
                    int up = (o + d_up[k] < W), dn = (o >= d_dn[k]);
                    if (up && dn) { if (g_up[k] < g_dn[k]) { c += g_up[k]; o += d_up[k]; } else { c += g_dn[k]; o -= d_dn[k]; } }
                    else if (up) { c += g_up[k]; o += d_up[k]; }
                    else if (dn) { c += g_dn[k]; o -= d_dn[k]; }
                    else { c += g_up[k] + g_dn[k]; o = o + d_up[k] - d_dn[k]; }
                }
            }
        }
        // pairs
        const int np = (n == 3) ? 6 : 24;
        for (size_t i = x0; i < x1; i++) {
            int ix[4]; for (int k = 0; k < n; k++) ix[k] = idxX[X[i].N[k] - 65536];
            for (size_t j = y0; j < y1; j++) {
                if (Y[j].SJ - Y[j].J - X[i].J < 2) continue;
                tot_pairs++;
                int h[4][4]; int any = 1;
                for (int a = 0; a < n; a++) { int rowany = 0; for (int b = 0; b < n; b++) { h[a][b] = head[(size_t)idxY[Y[j].N[b] - 65536] * nxv + ix[a]]; rowany |= h[a][b]; } if (!rowany) { any = 0; break; } }
                if (!any) continue;
                for (int p = 0; p < np; p++) {
                    const int *pm = (n == 3) ? PERM3[p] : PERM4[p];
                    int ok = 1; for (int a = 0; a < n; a++) if (!h[a][pm[a]]) { ok = 0; break; }
                    if (!ok) continue;
                    // skip duplicate pairings when Y has equal values
                    int dup = 0; for (int a = 0; a < n && !dup; a++) for (int b = a + 1; b < n; b++) if (X[i].N[a] == X[i].N[b] && pm[a] > pm[b]) { dup = 1; break; }
                    if (dup) continue;
                    uint32_t cc[4];
                    for (int e0 = h[0][pm[0]]; e0; e0 = pool[e0].next) for (int e1 = h[1][pm[1]]; e1; e1 = pool[e1].next) for (int e2 = h[2][pm[2]]; e2; e2 = pool[e2].next) {
                        cc[0] = pool[e0].c; cc[1] = pool[e1].c; cc[2] = pool[e2].c;
                        if (n == 3) { if ((int64_t)cc[0] + cc[1] + cc[2] == sigma) { tot_sum++; tot_out += verify(cc, stdout, tag); } }
                        else for (int e3 = h[3][pm[3]]; e3; e3 = pool[e3].next) { cc[3] = pool[e3].c;
                            if ((int64_t)cc[0] + cc[1] + cc[2] + cc[3] == sigma) { tot_sum++; tot_out += verify(cc, stdout, tag); } }
                    }
                }
            }
        }
        free(idxX); free(idxY); free(yv); free(head); free(pool);
    }
    fprintf(stderr, "X %zu Y %zu tags %d fibre %llu kept %llu pairs %llu sum_ok %llu out %llu\n", nx, ny, ntags, (unsigned long long)tot_fib,
            (unsigned long long)tot_keep, (unsigned long long)tot_pairs, (unsigned long long)tot_sum, (unsigned long long)tot_out);
    return 0;
}
