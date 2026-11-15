// hiscan: list the 4-byte chunk values whose four fractional parts (f1..f4) are all >= thr.
//   ./hiscan thr [lo hi]      output: one line per chunk  "hex f1 f2 f3 f4 map"
// build: gcc -O2 -march=native -ffp-contract=off -fno-fast-math -fopenmp -o short/hiscan short/hiscan.c -lm
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>
static double h2d(uint64_t h) { double d; memcpy(&d, &h, 8); return d; }
static inline double clr(double a) { uint64_t q; memcpy(&q, &a, 8); q &= ~0xFFULL; memcpy(&a, &q, 8); return a; }
static inline uint32_t ftol(double x) { if (x != x || isinf(x) || fabs(x) >= 9.223372036854775808e18) return 0; return (uint32_t)(int64_t)x; }
int main(int argc, char **argv) {
    double K_PI = h2d(0x3FD45F306DDA35A5ULL), K_E = h2d(0x3FD78B563618CF81ULL), C01 = h2d(0x3FB999999999999AULL);
    double thr = atof(argv[1]); uint64_t lo = argc > 2 ? strtoull(argv[2], 0, 16) : 0, hi = argc > 3 ? strtoull(argv[3], 0, 16) : 0x100000000ULL;
    #pragma omp parallel for schedule(static)
    for (int64_t blk = (int64_t)(lo >> 20); blk < (int64_t)((hi + 0xFFFFF) >> 20); blk++) {
        char buf[1 << 16]; int n = 0;
        for (uint64_t c = (uint64_t)blk << 20; c < (((uint64_t)blk + 1) << 20); c++) {
            if (c < lo || c >= hi) continue;
            double a1 = clr((double)c * K_PI), i1 = floor(a1); if (a1 - i1 < thr) continue;
            double a3 = clr((double)c * K_E), i3 = floor(a3); if (a3 - i3 < thr) continue;
            double a2 = clr(i1 * K_PI), i2 = floor(a2); if (a2 - i2 < thr) continue;
            double a4 = clr(i3 * K_E), i4 = floor(a4); if (a4 - i4 < thr) continue;
            double xm = i2; while (xm > 32767.0) xm *= C01;
            int mv = (int16_t)(ftol(xm) & 0xFFFF);
            n += snprintf(buf + n, sizeof(buf) - n, "%08llx %.17g %.17g %.17g %.17g %d\n", (unsigned long long)c, a1 - i1, a2 - i2, a3 - i3, a4 - i4, mv % 3);
            if (n > (int)sizeof(buf) - 256) { 
                #pragma omp critical
                fwrite(buf, 1, n, stdout);
                n = 0; }
        }
        if (n) {
            #pragma omp critical
            fwrite(buf, 1, n, stdout);
        }
    }
    return 0;
}
