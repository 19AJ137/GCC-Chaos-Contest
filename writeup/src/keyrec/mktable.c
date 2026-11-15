// mktable: chunk table (n3, n1, c) for a byte set, sorted by (n3, n1).  All chunks must be in lattice class 15 for f3 and f1.
// usage: mktable out.bin lo4 hi4  (bytes 1-3: 0x21-0x7E and 0xA1-0xDF; byte 4: lo4..hi4)
// build: gcc -O2 -march=native -ffp-contract=off -fno-fast-math -o mktable mktable.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef struct { uint32_t n3, n1, c; } Ent;
static double h2d(uint64_t h) { double d; memcpy(&d, &h, 8); return d; }
static inline uint32_t latN(uint32_t c, double K, int *ell) {
    double a = (double)c * K; uint64_t q; memcpy(&q, &a, 8); q &= ~0xFFULL; memcpy(&a, &q, 8);
    int e; frexp(a, &e); e -= 1; *ell = 44 - e;
    double f = a - floor(a); return (1u << *ell) + (uint32_t)llrint(ldexp(f, *ell));
}
static int cmp(const void *a, const void *b) { const Ent *x = a, *y = b; return x->n1 < y->n1 ? -1 : x->n1 > y->n1 ? 1 : (x->c < y->c ? -1 : x->c > y->c); }
int main(int argc, char **argv) {
    double K_PI = h2d(0x3fd45f306dda35a5ULL), K_E = h2d(0x3fd78b563618cf81ULL);
    int lo4 = strtol(argv[2], 0, 0), hi4 = strtol(argv[3], 0, 0);
    int bs[256], nb = 0; for (int b = 0x21; b <= 0x7E; b++) bs[nb++] = b; for (int b = 0xA1; b <= 0xDF; b++) bs[nb++] = b;
    size_t N = (size_t)nb * nb * nb * (hi4 - lo4 + 1); Ent *t = malloc(N * sizeof(Ent)); size_t *cnt = calloc(32770, sizeof(size_t));
    fprintf(stderr, "%zu chunks\n", N);
    for (int pass = 0; pass < 2; pass++) {
        if (pass == 1) { size_t s = 0; for (int i = 0; i < 32769; i++) { size_t c = cnt[i]; cnt[i] = s; s += c; } }
        for (int d = lo4; d <= hi4; d++) for (int i = 0; i < nb; i++) for (int j = 0; j < nb; j++) for (int k = 0; k < nb; k++) {
            uint32_t c = (uint32_t)bs[i] | (uint32_t)bs[j] << 8 | (uint32_t)bs[k] << 16 | (uint32_t)d << 24; int l3, l1;
            uint32_t n3 = latN(c, K_E, &l3), n1 = latN(c, K_PI, &l1);
            if (l3 != 15 || l1 != 15) { fprintf(stderr, "class %d %d for %08x\n", l3, l1, c); return 1; }
            if (pass == 0) cnt[n3 - 32768]++; else { Ent e = { n3, n1, c }; t[cnt[n3 - 32768]++] = e; }
        }
    }
    size_t s = 0; for (int i = 0; i < 32768; i++) { qsort(t + s, cnt[i] - s, sizeof(Ent), cmp); s = cnt[i]; }
    FILE *f = fopen(argv[1], "wb"); fwrite(t, sizeof(Ent), N, f); fclose(f); return 0;
}
