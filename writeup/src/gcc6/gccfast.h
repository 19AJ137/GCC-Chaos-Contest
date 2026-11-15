/* GCC chaos cipher (Gcc6_32.dll 1997-01-23) — fast C port. Compile with
 * -O2 -ffp-contract=off -fno-fast-math (SSE2 doubles == x87 with PC=53). */
#include <math.h>
#include <stdint.h>
#include <string.h>

typedef union { double d; uint8_t b[8]; uint64_t u; } dbl;

static const dbl K_INV_PI = {.u = 0x3fd45f306dda35a5ULL};
static const dbl K_INV_E  = {.u = 0x3fd78b563618cf81ULL};
static const dbl K_TENTH  = {.u = 0x3fb999999999999aULL};
static const dbl K_TENTH2 = {.u = 0x3fb99999999999a0ULL};
static const dbl K_X0DEF  = {.u = 0x3fe98f7c84e0cfe4ULL};
static const dbl K_3_9    = {.u = 0x400f333333333333ULL};
static const dbl K_1_9    = {.u = 0x3ffe666666666666ULL};
static const dbl K_0_4    = {.u = 0x3fd999999999999aULL};
static const dbl K_0_9    = {.u = 0x3feccccccccccccdULL};
static const dbl K_M014   = {.u = 0xbfc1eb851eb851ecULL};
static const dbl K_0_8    = {.u = 0x3fe999999999999aULL};
static const dbl K_1_72   = {.u = 0x3ffb851eb851eb85ULL};

typedef struct {
    int mode, warm;
    uint32_t S, T;
    uint8_t b0, b1;
    double k780, k788, k7b0, k7b8, k7b0_saved;
    double x, r, y, u, v, Q, A, C;
} gcc_state;

static inline double setlow(double x, uint8_t b0, uint8_t b1) {
    dbl t; t.d = x; t.b[0] = b0; t.b[1] = b1; return t.d;
}
static inline double clearlow(double x) { dbl t; t.d = x; t.b[0] = 0; return t.d; }

static inline uint32_t ftol32(double x) {
    if (!isfinite(x) || fabs(x) >= 9.223372036854775808e18) return 0;
    return (uint32_t)(uint64_t)(int64_t)x;
}

static inline void gcc_setup(gcc_state *s) {
    if (s->mode == 0) {
        double x = s->k780;
        if (x == 0.5 || !(0.0 < x)) { x = K_X0DEF.d; s->k780 = x; }
        s->x = x;
        s->r = s->k7b0 * K_TENTH2.d + K_3_9.d;
    } else if (s->mode == 1) {
        s->y = s->k780 * 2.0 - 1.0;
        s->r = s->k7b0 * K_TENTH2.d + K_1_9.d;
        s->y = setlow(s->y, s->b0, s->b1);
    } else {
        s->u = s->k780; s->v = s->k788;
        s->r = s->k7b0 * K_0_4.d + K_0_9.d;
        double Q = (s->k7b8 * K_M014.d - s->r * K_0_8.d) + K_1_72.d;
        s->Q = setlow(Q, s->b0, s->b1);
        s->A = setlow(Q * s->r + 1.0, s->b0, s->b1);
        s->C = setlow(1.0 - Q, s->b0, s->b1);
    }
    s->r = setlow(s->r, s->b0, s->b1);
}

static inline double gcc_step(gcc_state *s) {
    if (s->mode == 0) {
        s->x = setlow(((1.0 - s->x) * s->x) * s->r, s->b0, s->b1);
        return s->x;
    } else if (s->mode == 1) {
        double t = setlow(s->y * s->y - s->r, s->b0, s->b1);
        s->y = t; s->x = (1.0 + t) * 0.5; return s->x;
    } else {
        double a = s->Q * s->u, b = a * s->u;
        double nu = s->A * s->u - a * s->v;
        double nv = s->C * s->v + b;
        s->u = setlow(nu, s->b0, s->b1);
        s->v = setlow(nv, s->b0, s->b1);
        s->x = s->v; return s->x;
    }
}

static inline double fdiv_x87(double a, double b) { return a / b; } /* IEEE: inf/nan like x87 masked */

/* returns 0 ok, -1 if the DLL would leave the map uninitialised (empty key / negative mode) */
static int gcc_key_schedule(const uint8_t *key, int len, gcc_state *s) {
    double k780 = 1.0, k788 = 1.0, k7b0 = 1.0, k7b8 = 1.0, ip1, ip2 = 0, ip3, ip4, f1, f2, f3, f4;
    uint32_t S = 0, T = 0;
    int n = len & 0xffff; if (n >= 0x8000 || n == 0) return -1;
    int off = 0, rem = n;
    while (rem > 0) {
        int m = rem < 4 ? rem : 4;
        uint32_t chunk = 0;
        for (int i = 0; i < m; i++) chunk |= (uint32_t)key[off + i] << (8 * i);
        double dc = (double)chunk;
        f1 = modf(clearlow(dc * K_INV_PI.d), &ip1);
        f2 = modf(clearlow(ip1 * K_INV_PI.d), &ip2);
        k780 = (1.0 + f1) * k780;
        k788 = (f2 + 1.0) * k788;
        f3 = modf(clearlow(dc * K_INV_E.d), &ip3);
        f4 = modf(clearlow(ip3 * K_INV_E.d), &ip4);
        k7b0 = (1.0 + f3) * k7b0;
        k7b8 = (f4 + 1.0) * k7b8;
        S += ftol32(fdiv_x87(ip3, ip2)) + ftol32(ip1);
        T += ftol32(fdiv_x87(ip1, ip4)) + ftol32(ip3);
        off += m; rem -= 4;
    }
    double *ks[4] = {&k780, &k788, &k7b0, &k7b8};
    for (int i = 0; i < 4; i++) if (1.0 < *ks[i]) { double ip; *ks[i] = modf(*ks[i], &ip); S += ftol32(ip); }
    int warm = (int)(S % 79u); if (warm < 11) warm = 11;
    s->b0 = T & 0xff; s->b1 = (T >> 8) & 0xff;
    double v = ip2;
    if (!(v <= 32767.0)) { while (!(v <= 32767.0)) v = K_TENTH.d * v; }
    int16_t ax = (int16_t)(ftol32(v) & 0xffff);
    int mode = ax % 3;   /* C '%' truncates toward zero like idiv */
    if (mode < 0) return -1;
    s->mode = mode; s->warm = warm; s->S = S; s->T = T;
    s->k780 = k780; s->k788 = k788; s->k7b0 = k7b0; s->k7b8 = k7b8; s->k7b0_saved = k7b0;
    gcc_setup(s);
    for (int i = 0; i < warm; i++) gcc_step(s);
    return 0;
}

/* keystream word from output double: bytes (4,2,3,5); returns 1 if NOT chained (byte6 > 0xE6) */
static inline int gcc_ks(double x, uint32_t *ks) {
    dbl t; t.d = x;
    *ks = (uint32_t)t.b[4] | ((uint32_t)t.b[2] << 8) | ((uint32_t)t.b[3] << 16) | ((uint32_t)t.b[5] << 24);
    return t.b[6] > 0xE6;
}

/* decrypt nwords dwords of one 1024-byte block (File_cryption mode 2, first block) */
static inline void gcc_decrypt_words(gcc_state *s, const uint32_t *c, uint32_t *p, int nwords) {
    uint32_t chain = s->S;
    for (int i = 0; i < nwords; i++) {
        uint32_t ks; int nc = gcc_ks(gcc_step(s), &ks);
        if (nc) p[i] = c[i] ^ ks;
        else { p[i] = c[i] ^ ks ^ chain; chain = c[i]; }
    }
}
