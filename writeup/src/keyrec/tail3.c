// tail3: keys of 3 full chunks followed by a tail of 1..3 bytes (13..15 bytes).
// For every tail value (bytes in [blo, bhi]) whose i2 selects the wanted map, the contribution of the tail is divided out and
// the two sum-product systems of the 3 full chunks (2^-16 lattice) are enumerated (one element + quadratic closing).
// output lines:  X tail sigma J A Na x y      Y tail sigma SJ J A Na x y
// build: gcc -O2 -march=native -ffp-contract=off -fno-fast-math -fopenmp -o tail3 tail3.c -lm
// usage: tail3 S Tlo Thi map k7b0 half3 k780 half1 taillen blo bhi [part nparts]
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>
typedef unsigned __int128 u128;
static double K_PI, K_E, C01; static uint64_t MP, ME;
static double h2d(uint64_t h) { double d; memcpy(&d, &h, 8); return d; }
static inline double clr(double a) { uint64_t q; memcpy(&q, &a, 8); q &= ~0xFFULL; memcpy(&a, &q, 8); return a; }
static inline uint32_t ftol(double x) { if (x != x || isinf(x) || fabs(x) >= 9.223372036854775808e18) return 0; return (uint32_t)(int64_t)x; }
typedef struct { int64_t A; int64_t sig; int SJ; } Cand;
static int cmpc(const void *a, const void *b) { const Cand *x = a, *y = b; return x->A < y->A ? -1 : x->A > y->A ? 1 : (x->sig < y->sig ? -1 : x->sig > y->sig); }
static inline int64_t isqrt64(int64_t D) { int64_t q = (int64_t)sqrt((double)D); while (q * q > D) q--; while ((q + 1) * (q + 1) <= D) q++; return q; }
#define L 16
#define LO (1LL << L)
#define HI (2LL << L)
#define DMAX 11
#define NF 3

typedef struct { Cand *c; int n, cap; } CL;
static void push(CL *l, int64_t A, int64_t sig, int SJ) { if (l->n == l->cap) { l->cap = l->cap ? l->cap * 2 : 256; l->c = realloc(l->c, l->cap * sizeof(Cand)); } l->c[l->n].A = A; l->c[l->n].sig = sig; l->c[l->n].SJ = SJ; l->n++; }
static void bounds(double F, double *lo, double *hi) { double j = floor(F); *lo = ldexp(1.0 + (F - j), (int)j); *hi = pow(1.0 + F / NF, NF); }

static void enumerate(char sys, uint32_t tail, int J, long double Bp, double tol, CL *l, char **buf, size_t *len, size_t *cap, uint64_t *iters) {
    if (!l->n) return;
    qsort(l->c, l->n, sizeof(Cand), cmpc);
    int64_t Amax = l->c[l->n - 1].A, Amin = l->c[0].A;
    long double Bfull = Bp * 281474976710656.0L;                 // * 2^48 : Na * x * y
    for (int64_t Na = LO; Na <= Amax / 3 && Na < HI; Na++) {
        (*iters)++;
        double preal = (double)(Bfull / (long double)Na), w = preal * tol + 2e-5;
        int64_t rlo = (int64_t)ceil(preal - w), rhi = (int64_t)floor(preal + w);
        for (int64_t r = rlo; r <= rhi; r++) {
            if (r < Na * Na) continue;
            int64_t smin = (int64_t)ceil(2.0 * sqrt((double)r) - 1e-6), smax = Na + r / Na;
            if (smax + Na < Amin || smin + Na > Amax) continue;
            int a = 0, b = l->n; while (a < b) { int m = (a + b) >> 1; if (l->c[m].A < smin + Na) a = m + 1; else b = m; }
            for (; a < l->n && l->c[a].A <= smax + Na; a++) {
                int64_t s = l->c[a].A - Na, D = s * s - 4 * r; if (D < 0) continue;
                int64_t q = isqrt64(D); if (q * q != D) continue;
                int64_t x = (s - q) / 2, y = s - x; if (x < Na || y >= HI || x * y != r) continue;
                char line[200]; int o;
                if (sys == 'X') o = snprintf(line, sizeof line, "X %u %lld %d %lld %lld %lld %lld\n", tail, (long long)l->c[a].sig, J, (long long)l->c[a].A, (long long)Na, (long long)x, (long long)y);
                else o = snprintf(line, sizeof line, "Y %u %lld %d %d %lld %lld %lld %lld\n", tail, (long long)l->c[a].sig, l->c[a].SJ, J, (long long)l->c[a].A, (long long)Na, (long long)x, (long long)y);
                if (*len + o + 1 > *cap) { *cap = (*cap + o) * 2 + 4096; *buf = realloc(*buf, *cap); }
                memcpy(*buf + *len, line, o); *len += o;
            }
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 12) { fprintf(stderr, "usage: tail3 S Tlo Thi map k7b0 half3 k780 half1 taillen blo bhi [part nparts]\n"); return 1; }
    K_PI = h2d(0x3fd45f306dda35a5ULL); K_E = h2d(0x3fd78b563618cf81ULL); C01 = h2d(0x3fb999999999999aULL);
    MP = (0x3fd45f306dda35a5ULL & ((1ULL << 52) - 1)) | (1ULL << 52); ME = (0x3fd78b563618cf81ULL & ((1ULL << 52) - 1)) | (1ULL << 52);
    uint32_t S = strtoul(argv[1], 0, 16), Tlo = strtoul(argv[2], 0, 16), Thi = strtoul(argv[3], 0, 16); int wantmap = atoi(argv[4]);
    double k7b0 = atof(argv[5]), h3 = atof(argv[6]), k780 = atof(argv[7]), h1 = atof(argv[8]);
    int tl = atoi(argv[9]), blo = strtol(argv[10], 0, 0), bhi = strtol(argv[11], 0, 0);
    int part = argc > 13 ? atoi(argv[12]) : 0, nparts = argc > 13 ? atoi(argv[13]) : 1;
    int nb = bhi - blo + 1; int64_t ntail = 1; for (int i = 0; i < tl; i++) ntail *= nb;
    uint64_t tot_tails = 0, tot_iter = 0, tot_c = 0;
    #pragma omp parallel reduction(+:tot_tails,tot_iter,tot_c)
    {
        CL lx[16], ly[16]; memset(lx, 0, sizeof lx); memset(ly, 0, sizeof ly);
        char *buf = NULL; size_t len = 0, cap = 0;
        #pragma omp for schedule(dynamic, 64)
        for (int64_t ti = 0; ti < ntail; ti++) {
            if (ti % nparts != part) continue;
            uint32_t t = 0; int64_t v = ti; for (int i = 0; i < tl; i++) { t |= (uint32_t)(blo + v % nb) << (8 * i); v /= nb; }
            double a1 = clr((double)t * K_PI), i1 = floor(a1), f1 = a1 - i1, a2 = clr(i1 * K_PI), i2 = floor(a2);
            double a3 = clr((double)t * K_E), i3 = floor(a3), f3 = a3 - i3, a4 = clr(i3 * K_E), i4 = floor(a4);
            double xm = i2; while (xm > 32767.0) xm *= C01;
            int mv = (int16_t)(ftol(xm) & 0xFFFF); if (mv % 3 != wantmap) continue;
            tot_tails++;
            uint32_t r1 = ftol(i3 / i2), r3 = ftol(i1 / i4);
            uint32_t Sp = S - (uint32_t)i1 - r1; uint32_t Toff = (uint32_t)i3 + r3;
            for (int J = 0; J < 16; J++) lx[J].n = ly[J].n = 0;
            for (int w = 0; w < 2; w++) {
                int64_t Sw = (int64_t)Sp + ((int64_t)w << 32);
                int64_t slo = (int64_t)((Sw - 9 - 60) / K_PI) - 3, shi = (int64_t)((Sw - 9 - 4 + NF) / K_PI) + 3;
                if (shi < 0) continue; if (slo < 0) slo = 0; if (slo > 3LL * 0xFFFFFFFFLL) continue;
                for (int64_t sg = slo; sg <= shi; sg++) {
                    u128 p3 = (u128)ME * (u128)sg; int64_t fl3 = (int64_t)(p3 >> 54); uint64_t fr3 = (uint64_t)(p3 & (((u128)1 << 54) - 1));
                    int used = 0;
                    for (int j3 = 0; j3 < NF; j3++) {
                        int64_t I3 = fl3 - j3; uint32_t T16 = (uint32_t)((I3 + 2 * NF + Toff) & 0xFFFF);
                        if (T16 < Tlo || T16 > Thi) continue;
                        double F = j3 + (double)fr3 / 18014398509481984.0, lo, hi; bounds(F, &lo, &hi);
                        int64_t base = ((int64_t)NF << L) + ((int64_t)j3 << L) + (int64_t)(fr3 >> (54 - L));
                        int Ja = (int)ceil((lo - 2e-3) * (1.0 + f3) - k7b0), Jb = (int)floor((hi + 2e-3) * (1.0 + f3) - k7b0);
                        if (Ja < 1) Ja = 1; if (Jb > 15) Jb = 15;
                        for (int J = Ja; J <= Jb; J++) { used = 1; for (int d = 0; d < DMAX; d++) push(&lx[J], base - d, sg, 0); }
                    }
                    if (!used) continue;
                    u128 p1 = (u128)MP * (u128)sg; int64_t fl1 = (int64_t)(p1 >> 54); uint64_t fr1 = (uint64_t)(p1 & (((u128)1 << 54) - 1));
                    for (int j1 = 0; j1 < NF; j1++) {
                        int64_t I1 = fl1 - j1, SJ = Sw - 3 * NF - I1; if (SJ < 4 || SJ > 60) continue;
                        double F = j1 + (double)fr1 / 18014398509481984.0, lo, hi; bounds(F, &lo, &hi);
                        int64_t base = ((int64_t)NF << L) + ((int64_t)j1 << L) + (int64_t)(fr1 >> (54 - L));
                        int Ja = (int)ceil((lo - 2e-3) * (1.0 + f1) - k780), Jb = (int)floor((hi + 2e-3) * (1.0 + f1) - k780);
                        if (Ja < 1) Ja = 1; if (Jb > 15) Jb = 15;
                        for (int J = Ja; J <= Jb; J++) for (int d = 0; d < DMAX; d++) push(&ly[J], base - d, sg, (int)SJ);
                    }
                }
            }
            for (int J = 1; J < 16; J++) {
                tot_c += lx[J].n + ly[J].n;
                enumerate('X', t, J, ((long double)J + (long double)k7b0) / (1.0L + (long double)f3), h3 / (J + k7b0) + 3.6e-15, &lx[J], &buf, &len, &cap, &tot_iter);
                enumerate('Y', t, J, ((long double)J + (long double)k780) / (1.0L + (long double)f1), h1 / (J + k780) + 3.6e-15, &ly[J], &buf, &len, &cap, &tot_iter);
            }
            if (len > (1 << 20)) {
                #pragma omp critical
                { fwrite(buf, 1, len, stdout); }
                len = 0;
            }
        }
        #pragma omp critical
        { if (len) fwrite(buf, 1, len, stdout); }
    }
    fprintf(stderr, "tails with map %d: %llu of %lld; candidates %llu; iterations %llu\n", wantmap, (unsigned long long)tot_tails, (long long)ntail, (unsigned long long)tot_c, (unsigned long long)tot_iter);
    return 0;
}
