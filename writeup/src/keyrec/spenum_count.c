// spenum: enumerate the solutions of a "sum-product system"
//   N_1 <= ... <= N_n integers in [2^L, 2^(L+1)),  sum N_i = A (A in a candidate list),  prod N_i = (J+k) * 2^(nL) within a relative tolerance.
// The n-2 smallest elements are enumerated; the two largest are obtained from x+y = s, x*y = p (p must be an integer, s^2-4p a perfect square).
// Uniform lattice (all n elements on the same 2^-L lattice).  n = 3, 4, 5.
// build: gcc -O2 -march=native -ffp-contract=off -fno-fast-math -fopenmp -o spenum_count spenum_count.c -lm
// usage: spenum n L infile [lo_frac hi_frac]   (optional: restrict the smallest element to a fraction of its range, for scaled-down tests)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

typedef struct { int J; long double B; double tol; int nA; int64_t *A; int64_t Amin, Amax; } JB;
static int n, L; static JB *jb; static int nJ;
static double *inv;                 // inv[N] = 2^L / N
static int64_t LO, HI;

static inline int64_t isqrt64(int64_t D) { int64_t q = (int64_t)sqrt((double)D); while (q * q > D) q--; while ((q + 1) * (q + 1) <= D) q++; return q; }

typedef struct { char *buf; size_t len, cap; uint64_t nsurv, ntest, nsol, niter; } TB;
static void emit(TB *t, const char *s) { size_t l = strlen(s); if (t->len + l + 1 > t->cap) { t->cap = (t->cap + l) * 2 + 4096; t->buf = realloc(t->buf, t->cap); } memcpy(t->buf + t->len, s, l); t->len += l; }

// closing step: given the enumerated elements (sum Ns, last = Nl) and p_real (expected x*y), tolerance tol
static inline void closing(TB *t, const JB *j, const int64_t *N, int ne, int64_t Ns, int64_t Nl, double preal) {
    double w = preal * j->tol + 2e-5;
    int64_t rlo = (int64_t)ceil(preal - w), rhi = (int64_t)floor(preal + w);
    for (int64_t r = rlo; r <= rhi; r++) {
        if (r < Nl * Nl) continue;                       // need x >= Nl, y >= x
        t->nsurv++;
        int64_t smin = (int64_t)ceil(2.0 * sqrt((double)r) - 1e-6), smax = Nl + r / Nl;
        int64_t Alo = smin + Ns, Ahi = smax + Ns;
        if (Ahi < j->Amin || Alo > j->Amax) continue;
        int a = 0, b = j->nA;                            // lower_bound
        while (a < b) { int m = (a + b) >> 1; if (j->A[m] < Alo) a = m + 1; else b = m; }
        for (; a < j->nA && j->A[a] <= Ahi; a++) {
            int64_t s = j->A[a] - Ns, D = s * s - 4 * r; t->ntest++;
            if (D < 0) continue;
            int64_t q = isqrt64(D); if (q * q != D) continue;
            int64_t x = (s - q) / 2, y = s - x;
            if (x < Nl || y >= HI || x * y != r) continue;
            char line[256]; int o = snprintf(line, sizeof line, "%d %lld", j->J, (long long)j->A[a]);
            for (int i = 0; i < ne; i++) o += snprintf(line + o, sizeof line - o, " %lld", (long long)N[i]);
            snprintf(line + o, sizeof line - o, " %lld %lld\n", (long long)x, (long long)y);
            emit(t, line); t->nsol++;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: spenum n L infile [lo_frac hi_frac]\n"); return 1; }
    n = atoi(argv[1]); L = atoi(argv[2]); LO = 1LL << L; HI = 2LL << L;
    double f0 = argc > 5 ? atof(argv[4]) : 0.0, f1 = argc > 5 ? atof(argv[5]) : 1.0;
    FILE *fp = fopen(argv[3], "r"); if (!fp) { perror("in"); return 1; }
    if (fscanf(fp, "%d", &nJ) != 1) return 1;
    jb = calloc(nJ, sizeof(JB));
    for (int i = 0; i < nJ; i++) {
        if (fscanf(fp, "%d %Lf %lf %d", &jb[i].J, &jb[i].B, &jb[i].tol, &jb[i].nA) != 4) return 1;
        jb[i].A = malloc(sizeof(int64_t) * (jb[i].nA + 1));
        for (int a = 0; a < jb[i].nA; a++) { long long v; if (fscanf(fp, "%lld", &v) != 1) return 1; jb[i].A[a] = v; }
        jb[i].Amin = jb[i].A[0]; jb[i].Amax = jb[i].A[jb[i].nA - 1];
    }
    fclose(fp);
    inv = malloc(sizeof(double) * HI);
    for (int64_t v = LO; v < HI; v++) inv[v] = (double)((long double)LO / (long double)v);
    uint64_t tot_surv = 0, tot_test = 0, tot_sol = 0, tot_iter = 0; double t0 = omp_get_wtime();
    for (int ji = 0; ji < nJ; ji++) {
        const JB *j = &jb[ji];
        long double Bfull = j->B * powl(2.0L, 2.0L * L);            // (J+k) * 2^(2L): x*y = Bfull * prod_enum (2^L / N_i)
        int64_t amax = j->Amax / n, a_lo = LO + (int64_t)((amax + 1 - LO) * f0), a_hi = LO + (int64_t)((amax + 1 - LO) * f1);
        #pragma omp parallel
        {
            TB t = {0};
            #pragma omp for schedule(dynamic, 16)
            for (int64_t Na = a_lo; Na < a_hi; Na++) {
                int64_t N[3]; N[0] = Na;
                if (n == 3) { double preal = (double)(Bfull * (long double)inv[Na]); t.niter++; closing(&t, j, N, 1, Na, Na, preal); continue; }
                int64_t bmax = (j->Amax - Na) / (n - 1);
                for (int64_t Nb = Na; Nb <= bmax && Nb < HI; Nb++) {
                    N[1] = Nb;
                    long double Bq_l = Bfull * (long double)inv[Na] * (long double)inv[Nb];
                    if (n == 4) { t.niter++; closing(&t, j, N, 2, Na + Nb, Nb, (double)Bq_l); continue; }
                    double Bq = (double)Bq_l;                      // x*y*Nc/2^L = Bq
                    int64_t cmax = (j->Amax - Na - Nb) / 3; if (cmax >= HI) cmax = HI - 1;
                    // x,y >= Nc  =>  Nc^3 <= Bq * 2^L ;  x,y < 2^(L+1)  =>  Nc > Bq*2^L / 2^(2L+2)
                    int64_t c2 = (int64_t)cbrt(Bq * (double)LO) + 1; if (c2 < cmax) cmax = c2;
                    int64_t cmin = (int64_t)(Bq * (double)LO / ((double)HI * (double)HI)); if (cmin < Nb) cmin = Nb;
                    int64_t c3 = j->Amin - Na - Nb - 2 * (HI - 1); if (c3 > cmin) cmin = c3;
                    double tol = j->tol; int64_t Ns2 = Na + Nb;
                    t.niter += (cmax >= cmin) ? (uint64_t)(cmax - cmin + 1) : 0;
                    if (getenv("COUNT_ONLY")) continue;
                    for (int64_t Nc = cmin; Nc <= cmax; Nc++) {
                        double preal = Bq * inv[Nc];
                        double r = nearbyint(preal);
                        if (fabs(preal - r) > preal * tol + 2e-5) continue;
                        N[2] = Nc; closing(&t, j, N, 3, Ns2 + Nc, Nc, preal);
                    }
                }
            }
            #pragma omp critical
            { if (t.len) fwrite(t.buf, 1, t.len, stdout); tot_surv += t.nsurv; tot_test += t.ntest; tot_sol += t.nsol; tot_iter += t.niter; free(t.buf); }
        }
        fprintf(stderr, "J=%d nA=%d done: iter %llu surv %llu tests %llu sol %llu  %.1fs\n", j->J, j->nA, (unsigned long long)tot_iter,
                (unsigned long long)tot_surv, (unsigned long long)tot_test, (unsigned long long)tot_sol, omp_get_wtime() - t0);
    }
    return 0;
}
