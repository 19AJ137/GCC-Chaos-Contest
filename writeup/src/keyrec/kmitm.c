// kmitm.c — 「語を辞書から取る」形の鍵 v の候補を、中間一致で調べる（11-14、3 回目）。
//   鍵 = 固定のバイト列（IISI・区切り）と、語を入れる枠（スロット）の並び。枠ごとに語の一覧を渡す。
//   枠を 左の組（1〜2 個）・右の組（1〜2 個）・外側（0〜1 個）に分ける。1 つの塊が 2 つの組にまたがる分け方は使えない。
//   条件: S = (塊ごとの s の和) + (積の整数部 4 つの和 SJ)。左右の組の s の和を A_L・A_R とすると、
//         A_L + A_R + A_fix + A_outer ≡ S − SJ (mod 2^32)、SJ は [jlo, jhi]。
//         左の表・右の表を A の順に並べ、外側の語ごとに 2 本の表を 1 回ずつなぞって、A が合う組だけ積 P1 = Π(1+f1) を見る。
//         P1 の小数部が k780 のセル（2 通り）に入った組を、候補として鍵の 16 進で出す。最終の判定は呼び出し側（鍵スケジュール全体と鍵ストリーム）。
//   v2（11-14）: 大きいほうの表に、幅 256 の区画ごとの「要素があるか」のビット列と先頭位置の表を付ける。小さいほうの表を 1 回なぞり、
//         相手の区画が空なら、次に要素のある区画に届くところまで読み飛ばす（大きいほうの表は、要素のある区画だけを見る）。結果は v1（kmitm_v1.c）と同じ。
//   gcc -O2 -march=native -ffp-contract=off -fno-fast-math -fopenmp -o kmitm kmitm.c -lm
// 入力（テキスト）:
//   T S(16進) c1 h1 c1b h1b pad      k780 のセル 2 つ（中心と半幅）と、積の順序の違いに対する余裕
//   J jlo jhi
//   K L 16進                          鍵のひな形（枠の中のバイトは何でもよい）
//   W id len count path               語の一覧（count*len バイトのバイナリ）
//   G g ns off list [off list]        g = 0 左、1 右、2 外側。ns = 枠の数（外側は 0 か 1）
//   R                                 実行
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <immintrin.h>

static double K_PI, K_E;
static inline double clr(double x) { uint64_t u; memcpy(&u, &x, 8); u &= ~0xFFULL; memcpy(&x, &u, 8); return x; }
static inline uint32_t ftol32(double x) { if (x != x || isinf(x) || fabs(x) >= 9223372036854775808.0) return 0; return (uint32_t)(int64_t)x; }
// 塊 c の s（S への寄与）と 1 + f1
static inline void chunk_sp(uint32_t c, uint32_t *s, double *p) {
    double x1 = clr((double)c * K_PI), i1 = floor(x1), f1 = x1 - i1;
    double x2 = clr(i1 * K_PI), i2 = floor(x2);
    double x3 = clr((double)c * K_E), i3 = floor(x3);
    *s = ftol32(i3 / i2) + ftol32(i1); *p = 1.0 + f1;
}

typedef struct { int len, count; uint8_t *bytes; } WList;
typedef struct { int off, list; } Slot;
typedef struct { uint32_t A, idx; double P; } Ent;
typedef struct { int ns; Slot s[2]; uint32_t *a[2]; double *p[2]; int nmulti; int mch[16]; } Group;

static int L, NCH; static uint8_t TPL[64];
static WList WL[16]; static Group G[3];
static uint32_t S_T, JLO, JHI; static double C1, H1, C1B, H1B, PAD;
static uint32_t A_FIX; static double P_FIX;

static inline uint32_t chunk_val(const uint8_t *key, int i) { uint32_t c = 0; for (int j = 0; j < 4 && 4 * i + j < L; j++) c |= (uint32_t)key[4 * i + j] << (8 * j); return c; }
static int slot_touch(const Slot *s, int ch) { int len = WL[s->list].len; return s->off < 4 * ch + 4 && s->off + len > 4 * ch; }

static void die(const char *m) { fprintf(stderr, "kmitm: %s\n", m); exit(2); }

static void setup(void) {
    NCH = (L + 3) / 4; A_FIX = 0; P_FIX = 1.0;
    int owner[16], mask[16];
    for (int i = 0; i < NCH; i++) {
        owner[i] = -1; mask[i] = 0;
        for (int g = 0; g < 3; g++) for (int k = 0; k < G[g].ns; k++) if (slot_touch(&G[g].s[k], i)) {
            if (owner[i] >= 0 && owner[i] != g) die("a chunk touches two groups");
            owner[i] = g; mask[i] |= 1 << k;
        }
        if (owner[i] < 0) { uint32_t s; double p; chunk_sp(chunk_val(TPL, i), &s, &p); A_FIX += s; P_FIX *= p; }
    }
    for (int g = 0; g < 3; g++) {
        Group *gr = &G[g]; gr->nmulti = 0;
        for (int i = 0; i < NCH; i++) if (owner[i] == g && mask[i] == 3) gr->mch[gr->nmulti++] = i;
        if (g == 2 && gr->nmulti) die("outer slot shares a chunk");
        for (int k = 0; k < gr->ns; k++) {
            WList *w = &WL[gr->s[k].list]; int n = w->count;
            gr->a[k] = malloc(sizeof(uint32_t) * n); gr->p[k] = malloc(sizeof(double) * n);
            #pragma omp parallel for schedule(static)
            for (int x = 0; x < n; x++) {
                uint8_t key[64]; memcpy(key, TPL, L); memcpy(key + gr->s[k].off, w->bytes + (size_t)x * w->len, w->len);
                uint32_t a = 0; double p = 1.0;
                for (int i = 0; i < NCH; i++) if (owner[i] == g && mask[i] == (1 << k)) { uint32_t s; double q; chunk_sp(chunk_val(key, i), &s, &q); a += s; p *= q; }
                gr->a[k][x] = a; gr->p[k][x] = p;
            }
        }
    }
}

static int cmp_ent(const void *x, const void *y) { uint32_t a = ((const Ent *)x)->A, b = ((const Ent *)y)->A; return a < b ? -1 : a > b; }

// 組の表（A の順）
static Ent *build(Group *gr, size_t *nout) {
    int n0 = WL[gr->s[0].list].count, n1 = gr->ns == 2 ? WL[gr->s[1].list].count : 1;
    size_t n = (size_t)n0 * n1; Ent *T = malloc(sizeof(Ent) * (n ? n : 1)); if (!T) die("out of memory");
    size_t *cnt = calloc(65538, sizeof(size_t));
    int len0 = WL[gr->s[0].list].len, len1 = gr->ns == 2 ? WL[gr->s[1].list].len : 0;
    for (int pass = 0; pass < 2; pass++) {
        #pragma omp parallel
        {
            size_t *lc = pass == 0 ? calloc(65538, sizeof(size_t)) : NULL; uint8_t key[64]; memcpy(key, TPL, L);
            #pragma omp for schedule(dynamic, 16)
            for (int w0 = 0; w0 < n0; w0++) {
                uint32_t a0 = gr->a[0][w0]; double p0 = gr->p[0][w0];
                if (gr->nmulti) memcpy(key + gr->s[0].off, WL[gr->s[0].list].bytes + (size_t)w0 * len0, len0);
                for (int w1 = 0; w1 < n1; w1++) {
                    uint32_t A = a0; double P = p0;
                    if (gr->ns == 2) { A += gr->a[1][w1]; P *= gr->p[1][w1]; }
                    if (gr->nmulti) {
                        memcpy(key + gr->s[1].off, WL[gr->s[1].list].bytes + (size_t)w1 * len1, len1);
                        for (int m = 0; m < gr->nmulti; m++) { uint32_t s; double q; chunk_sp(chunk_val(key, gr->mch[m]), &s, &q); A += s; P *= q; }
                    }
                    if (pass == 0) lc[(A >> 16) + 1]++;
                    else { size_t pos;
                        #pragma omp atomic capture
                        pos = cnt[A >> 16]++;
                        T[pos].A = A; T[pos].idx = (uint32_t)((size_t)w0 * n1 + w1); T[pos].P = P; }
                }
            }
            if (pass == 0) {
                #pragma omp critical
                for (int b = 0; b < 65538; b++) cnt[b] += lc[b];
                free(lc);
            }
        }
        if (pass == 0) { for (int b = 1; b < 65538; b++) cnt[b] += cnt[b - 1]; }      // cnt[b] = b 番目のバケツの先頭
    }
    // 2 回目のあと cnt[b] は b 番目のバケツの末尾。バケツごとに A の順に並べる
    #pragma omp parallel for schedule(dynamic, 64)
    for (int b = 0; b < 65536; b++) { size_t e = cnt[b], s = b ? cnt[b - 1] : 0; if (e > s + 1) qsort(T + s, e - s, sizeof(Ent), cmp_ent); }
    free(cnt); *nout = n; return T;
}

typedef struct { int o; uint32_t il, ir; } Cand;
static Cand *CANDS; static size_t NCAND, CAPC; static unsigned long long NCHK;

static size_t lower_bound(const Ent *T, size_t n, uint64_t v) { size_t lo = 0, hi = n; while (lo < hi) { size_t m = (lo + hi) / 2; if (T[m].A < v) lo = m + 1; else hi = m; } return lo; }

static uint64_t *BM64; static uint32_t *START, *NXT; static int SH; static size_t NCELL;
// 表 T（A の順）に、幅 2^SH の区画ごとのビット列（要素があるか）と、区画の先頭位置を付ける。2^SH >= W
static void build_index(const Ent *T, size_t n, uint32_t W) {
    SH = 8; while ((1u << SH) < W) SH++;
    size_t ncell = (size_t)1 << (32 - SH), j = 0, ng = ncell / 64; NCELL = ncell;
    BM64 = calloc(ng + 2, 8); START = malloc(sizeof(uint32_t) * (ncell + 2)); NXT = malloc(sizeof(uint32_t) * (ng + 2)); if (!BM64 || !START || !NXT) die("out of memory");
    for (size_t c = 0; c < ncell; c++) {
        START[c] = (uint32_t)j; uint64_t lim = ((uint64_t)(c + 1)) << SH; size_t j0 = j;
        while (j < n && T[j].A < lim) j++;
        if (j > j0) BM64[c >> 6] |= 1ULL << (c & 63);
    }
    START[ncell] = START[ncell + 1] = (uint32_t)n;
    NXT[ng] = NXT[ng + 1] = (uint32_t)ncell;                      // NXT[g] = 64g 番目以降で最初に要素のある区画（なければ ncell）
    for (size_t g = ng; g-- > 0;) NXT[g] = BM64[g] ? (uint32_t)(g * 64 + __builtin_ctzll(BM64[g])) : NXT[g + 1];
}
// c 番目以降で最初に要素のある区画（なければ NCELL）
static inline size_t next_occ(size_t c) { uint64_t w = BM64[c >> 6] >> (c & 63); return w ? c + __builtin_ctzll(w) : NXT[(c >> 6) + 1]; }
static size_t lower_bound32(const uint32_t *A, size_t n, uint64_t v) { size_t lo = 0, hi = n; while (lo < hi) { size_t m = (lo + hi) / 2; if (A[m] < v) lo = m + 1; else hi = m; } return lo; }
#define CHECK(j) do { double x = T[j].P * pr, f = x - floor(x); chk++; \
    if (fabs(f - C1) <= H1 + PAD || fabs(f - C1B) <= H1B + PAD) { if (no < cap) { out[no].o = o; out[no].il = T[j].idx; out[no].ir = S[r].idx; } no++; } } while (0)
// r を、(uint32_t)(Cp − SA[r−1]) < tu のあいだ下げる（rlo まで）。この範囲では u は r が小さいほど大きい
static inline size_t skip_down(const uint32_t *SA, size_t r, size_t rlo, uint32_t Cp, uint32_t tu) {
#if defined(__AVX512F__)
    const __m512i cp = _mm512_set1_epi32((int)Cp), tv = _mm512_set1_epi32((int)tu);
    while (r >= rlo + 16) {
        __mmask16 m = _mm512_cmplt_epu32_mask(_mm512_sub_epi32(cp, _mm512_loadu_si512((const void *)(SA + r - 16))), tv);
        if (m == 0xFFFF) { r -= 16; continue; }
        return r - (size_t)__builtin_popcount((unsigned)m);       // u < tu のものは番号の大きい側にまとまっている
    }
#elif defined(__AVX2__)
    const __m256i cp = _mm256_set1_epi32((int)Cp), tv = _mm256_set1_epi32((int)(tu ^ 0x80000000u)), sg = _mm256_set1_epi32((int)0x80000000u);
    while (r >= rlo + 8) {
        __m256i u = _mm256_xor_si256(_mm256_sub_epi32(cp, _mm256_loadu_si256((const __m256i *)(SA + r - 8))), sg);
        unsigned m = (unsigned)_mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(tv, u)));
        if (m == 0xFF) { r -= 8; continue; }
        return r - (size_t)__builtin_popcount(m);
    }
#endif
    while (r > rlo && (uint32_t)(Cp - SA[r - 1]) < tu) r--;
    return r;
}
// 小さいほうの表 S（キーは SA）の rhi 番目から rlo 番目までを、番号の大きいほうからなぞる（この範囲では u = Cp − A_s (mod 2^32) が折り返さずに増える）。
// 大きいほうの表 T は A が [u, u+W-1] のもの。相手の区画が空のあいだは、次に要素のある区画に届くまで S を読み飛ばす
static void sweep(const Ent *T, size_t n, const uint32_t *SA, const Ent *S, size_t rhi, size_t rlo, uint32_t Cp, uint32_t W, double pc, int o, Cand *out, size_t *nout, size_t cap, unsigned long long *nchk) {
    unsigned long long chk = 0; size_t no = *nout, lo = 0; const uint32_t ulim = 0xFFFFFFFFu - (W - 1);
    size_t r = rhi + 1;                                          // これから見るのは r − 1
    while (r > rlo) {
        r--;
        uint32_t u = Cp - SA[r];
        if (__builtin_expect(u > ulim, 0)) {                     // 窓が 2^32 をまたぐ（まれ）
            double pr = S[r].P * pc; uint64_t tw = (uint64_t)u + W - 1 - 0x100000000ULL;
            for (size_t j = lower_bound(T, n, u); j < n; j++) CHECK(j);
            for (size_t j = 0; j < n && T[j].A <= tw; j++) CHECK(j);
            continue;
        }
        size_t c = u >> SH, cn = next_occ(c);
        if (cn > c + 1) {                                        // 区画 c・c+1 は空。u が区画 cn−1 に届くまで読み飛ばす
            uint32_t tu = (uint32_t)(((cn >= NCELL ? NCELL : cn) - 1) << SH);      // cn = NCELL（もう要素がない）なら、残るのは折り返しの窓だけ
            r = skip_down(SA, r, rlo, Cp, tu);
            continue;
        }
        size_t s = START[c]; if (s > lo) lo = s;
        while (lo < n && T[lo].A < u) lo++;
        uint32_t top = u + (W - 1);
        if (lo < n && T[lo].A <= top) { double pr = S[r].P * pc; for (size_t j = lo; j < n && T[j].A <= top; j++) CHECK(j); }
    }
    *nout = no; *nchk += chk;
}
// S を「u が小さい順」に並べたときの k0〜k1 番目（k = 0 は A_s <= Cp のうち最大のもの）
static void merge(const Ent *T, size_t n, const uint32_t *SA, const Ent *S, size_t ns, uint32_t Cp, uint32_t W, double pc, int o, size_t k0, size_t k1, Cand *out, size_t *nout, size_t cap, unsigned long long *nchk) {
    if (!n || !ns || k1 <= k0) return;
    size_t i0 = lower_bound32(SA, ns, (uint64_t)Cp + 1);         // A_s <= Cp の個数。k < i0 は r = i0−1−k、k >= i0 は r = ns−1−(k−i0)
    if (k0 < i0) { size_t ke = k1 < i0 ? k1 : i0; sweep(T, n, SA, S, i0 - 1 - k0, i0 - ke, Cp, W, pc, o, out, nout, cap, nchk); }
    if (k1 > i0) { size_t kb = k0 > i0 ? k0 : i0; sweep(T, n, SA, S, ns - 1 - (kb - i0), ns - (k1 - i0), Cp, W, pc, o, out, nout, cap, nchk); }
}

static void put_word(uint8_t *key, const Slot *s, uint32_t w) { WList *l = &WL[s->list]; memcpy(key + s->off, l->bytes + (size_t)w * l->len, l->len); }

int main(int argc, char **argv) {
    uint64_t a = 0x3FD45F306DDA35A5ULL, b = 0x3FD78B563618CF81ULL; memcpy(&K_PI, &a, 8); memcpy(&K_E, &b, 8);
    if (argc < 2) die("usage: kmitm spec");
    FILE *f = fopen(argv[1], "r"); if (!f) die("cannot open spec");
    char line[4096]; PAD = 2e-12; JLO = 0; JHI = 4092;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'T') { if (sscanf(line + 1, "%x %lf %lf %lf %lf %lf", &S_T, &C1, &H1, &C1B, &H1B, &PAD) < 5) die("bad T"); }
        else if (line[0] == 'J') sscanf(line + 1, "%u %u", &JLO, &JHI);
        else if (line[0] == 'K') { char hex[200]; sscanf(line + 1, "%d %199s", &L, hex); for (int i = 0; i < L; i++) { unsigned v; sscanf(hex + 2 * i, "%2x", &v); TPL[i] = v; } }
        else if (line[0] == 'W') { int id, len, cnt; char path[2048]; sscanf(line + 1, "%d %d %d %2047s", &id, &len, &cnt, path);
            WL[id].len = len; WL[id].count = cnt; WL[id].bytes = malloc((size_t)len * cnt + 1); FILE *w = fopen(path, "rb"); if (!w || fread(WL[id].bytes, len, cnt, w) != (size_t)cnt) die("bad word list"); fclose(w); }
        else if (line[0] == 'G') { int g, ns, o0 = 0, l0 = 0, o1 = 0, l1 = 0; int n = sscanf(line + 1, "%d %d %d %d %d %d", &g, &ns, &o0, &l0, &o1, &l1); if (n < 2) die("bad G");
            G[g].ns = ns; G[g].s[0].off = o0; G[g].s[0].list = l0; G[g].s[1].off = o1; G[g].s[1].list = l1; }
        else if (line[0] == 'R') break;
    }
    fclose(f);
    if (G[0].ns < 1 || G[1].ns < 1) die("need left and right groups");
    double t0 = omp_get_wtime();
    setup();
    size_t nL, nR; Ent *TL = build(&G[0], &nL), *TR = build(&G[1], &nR);
    uint32_t W = JHI - JLO + 1;
    int swapped = nL < nR;                                        // 大きいほうを T（区画の表つき）、小さいほうを S（なぞる側）にする
    const Ent *T = swapped ? TR : TL, *S = swapped ? TL : TR; size_t nT = swapped ? nR : nL, nS = swapped ? nL : nR;
    uint32_t *SA = malloc(sizeof(uint32_t) * (nS ? nS : 1)); if (!SA) die("out of memory");
    for (size_t i = 0; i < nS; i++) SA[i] = S[i].A;
    build_index(T, nT, W);
    double t1 = omp_get_wtime();
    CAPC = 1 << 20; CANDS = malloc(sizeof(Cand) * CAPC);
    int no = G[2].ns ? WL[G[2].s[0].list].count : 1;
    if (G[2].ns) {
        #pragma omp parallel
        {
            Cand *loc = malloc(sizeof(Cand) * 4096); unsigned long long chk = 0;
            #pragma omp for schedule(dynamic, 1)
            for (int o = 0; o < no; o++) {
                size_t nl = 0; uint32_t Cp = S_T - A_FIX - G[2].a[0][o] - JHI;
                merge(T, nT, SA, S, nS, Cp, W, P_FIX * G[2].p[0][o], o, 0, nS, loc, &nl, 4096, &chk);
                if (nl) {
                    #pragma omp critical
                    { for (size_t i = 0; i < nl && i < 4096; i++) if (NCAND < CAPC) CANDS[NCAND++] = loc[i]; else NCAND++; }
                }
            }
            #pragma omp atomic
            NCHK += chk;
            free(loc);
        }
    } else {
        int nseg = 256; uint32_t Cp = S_T - A_FIX - JHI;
        #pragma omp parallel
        {
            Cand *loc = malloc(sizeof(Cand) * 4096); unsigned long long chk = 0;
            #pragma omp for schedule(dynamic, 1)
            for (int sgi = 0; sgi < nseg; sgi++) {
                size_t k0 = nS * (size_t)sgi / nseg, k1 = nS * (size_t)(sgi + 1) / nseg, nl = 0;
                if (k1 > k0) merge(T, nT, SA, S, nS, Cp, W, P_FIX, 0, k0, k1, loc, &nl, 4096, &chk);
                if (nl) {
                    #pragma omp critical
                    { for (size_t i = 0; i < nl && i < 4096; i++) if (NCAND < CAPC) CANDS[NCAND++] = loc[i]; else NCAND++; }
                }
            }
            #pragma omp atomic
            NCHK += chk;
            free(loc);
        }
    }
    double t2 = omp_get_wtime();
    size_t shown = NCAND < CAPC ? NCAND : CAPC;
    for (size_t i = 0; i < shown; i++) {
        uint8_t key[64]; memcpy(key, TPL, L);
        int n1l = G[0].ns == 2 ? WL[G[0].s[1].list].count : 1, n1r = G[1].ns == 2 ? WL[G[1].s[1].list].count : 1;
        uint32_t xl = swapped ? CANDS[i].ir : CANDS[i].il, xr = swapped ? CANDS[i].il : CANDS[i].ir;
        put_word(key, &G[0].s[0], xl / n1l); if (G[0].ns == 2) put_word(key, &G[0].s[1], xl % n1l);
        put_word(key, &G[1].s[0], xr / n1r); if (G[1].ns == 2) put_word(key, &G[1].s[1], xr % n1r);
        if (G[2].ns) put_word(key, &G[2].s[0], CANDS[i].o);
        printf("CAND "); for (int j = 0; j < L; j++) printf("%02x", key[j]); printf("\n");
    }
    printf("STAT nL %zu nR %zu outer %d checks %llu cands %zu build %.1f merge %.1f\n", nL, nR, no, NCHK, NCAND, t1 - t0, t2 - t1);
    return 0;
}
