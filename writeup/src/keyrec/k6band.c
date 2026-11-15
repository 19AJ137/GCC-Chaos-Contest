// k6band.c — 語 6 個の形（2 文字の語 6 個を「・」でつないだ 34 バイト）のための中間一致（3 個＋3 個）。帯ごとに貯めて引く版。docs/password.md §8.10（2026-11-15）
//   貯める側 A: 表 W の単位 3 個の組（順不同。i ≤ j ≤ k）。s の和（桁あふれなしの本当の和）の狭い帯ごとに数万組を、積 P3 = Π(1+f3) で引ける小さい表（2 次キャッシュに載る大きさ）にする。
//   探す側 B:   表 M の単位 2 個（順不同）と表 Z の単位 1 個。B は「相手の A の s の和の下限 u = S' − jmax − sB」がどの帯に入るかで 1 つの帯にだけ割り当てる。
//               帯 [u, u+β) の表には s の和が [u, u+β+(jmax−jmin)] の A を入れる（隣の帯と重なる）。B ごとに、積の整数部 J3 の候補すべてについて P3 の目標 (J3 + c3) / P3B を作って表を引く（8 個ずつまとめて引く）。
//   当たった組は、k7b0 のセル、k780 のセル（2 通り）、S（積の整数部 4 つを足して一致）、b0b1 の範囲を確かめて CAND として出す。最後の判定は呼ぶ側（k6.py）で行う。
//   単位 = 鍵の中の 1 語が占める塊（1〜2 個）。単位ごとに s（塊の s の和。桁あふれなし）・t・P1〜P4（Π(1+f)）を渡す。
//   取り出し方: 単位 2 個の和の表（s の順）を W と M について作り、3 個目の単位ごとに読む位置を覚えて、帯を s の小さい側から順に進める。
//   表: 1 マス 8 バイト = 4 バイトの枠 2 つ（札 = P3 を刻んだ値の下位ビット、組の番号）。同じマスに 3 組以上入ったら別の列につなぐ。組の中身（i, j, k）は別の配列。
//   コンパイル: gcc -O3 -march=native -ffp-contract=off -fno-fast-math -fopenmp -o k6band k6band.c -lm
//   使い方: k6band 指定ファイル
//     指定ファイルの行: T S(16進) c1 h1 c1b h1b c3 h3 Tlo Thi slack / J jmin jmax / W n path / M n path / Z n path / N 帯の組数 表のビット数 区間の数 / L 秒（時間の上限。0 = なし） / R a b（する区間。両端を含む。何行でも）
//     出力: NRANGE 区間の数 / CAND i j k a b z / DONE 区間 … / STAT …
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <omp.h>
#include <sys/mman.h>
#ifdef __AVX512DQ__
#include <immintrin.h>
#endif

typedef struct { uint64_t s; uint32_t t, pad; double P[4]; } unit_t;      // 48 バイト
typedef struct { uint64_t sk; double p3; } pair_t;                         // sk = s << 28 | i << 14 | j（i ≤ j）
#define IB 14
#define IM ((1u << IB) - 1)
#define FLAG 0x80000000u
#define END 0x7FFFFFFFu
#define BLK 128
typedef struct { uint64_t Sp, U0, U1; } range_t;
typedef struct {
    uint32_t *tab, *cs, *cn, *bz; uint64_t *tid, *ovf;
    uint64_t nA, nB, nprobe, nhit, nl3, nc1, nS, ncand, nband, nsplit, nslow; double tbuild, tprobe;
} thr_t;

static unit_t *W, *M, *Z; static int nW, nM, nZ;
static pair_t *pairsW, *pairsM; static uint32_t nPW, nPM;
static uint64_t *Ws, *Zs; static double *Wp3, *Zp3;
static uint32_t TS; static double c1, h1, c1b, h1b, c3, h3, slack, tol; static int Tlo, Thi, jmin, jmax, jw;
static int ntarget = 60000, dbits = 18, nrange_req = 2048, ib, tb, tg, xb; static uint32_t nb, cap, gmask, xmask; static uint64_t imask; static double scale, tlimit = 0;
static const double maxq = 4294967295.0, margin = 2e-3;

static unit_t *load(const char *path, int n) {
    FILE *f = fopen(path, "rb"); if (!f) { perror(path); exit(1); }
    unit_t *u = malloc(sizeof(unit_t) * (size_t)n); if (fread(u, sizeof(unit_t), n, f) != (size_t)n) { fprintf(stderr, "short read %s\n", path); exit(1); } fclose(f); return u;
}
static void *big_alloc(size_t bytes) {                 // 2 MB 境界にそろえて、大きいページを頼む（表を引くときの TLB の外れを減らす）
    size_t sz = (bytes + (2u << 20) - 1) & ~((size_t)(2u << 20) - 1); void *p = aligned_alloc(2u << 20, sz); if (!p) { fprintf(stderr, "out of memory\n"); exit(1); }
    madvise(p, sz, MADV_HUGEPAGE); return p;
}
static int cmp_pair(const void *a, const void *b) { uint64_t x = ((const pair_t *)a)->sk, y = ((const pair_t *)b)->sk; return x < y ? -1 : x > y; }
static int cmp_u64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }
static pair_t *make_pairs(const unit_t *U, int n, uint32_t *np) {
    uint32_t m = (uint32_t)((uint64_t)n * (n + 1) / 2), q = 0; pair_t *p = big_alloc(sizeof(pair_t) * (size_t)m);
    for (int i = 0; i < n; i++) for (int j = i; j < n; j++) { p[q].sk = ((U[i].s + U[j].s) << 28) | ((uint64_t)i << IB) | (uint64_t)j; p[q].p3 = U[i].P[2] * U[j].P[2]; q++; }
    qsort(p, m, sizeof(pair_t), cmp_pair); *np = m; return p;
}
static inline uint32_t lbp(const pair_t *pr, uint32_t n, int64_t sval) {       // s >= sval になる最初の位置
    if (sval <= 0) return 0;
    uint64_t key = (uint64_t)sval << 28; uint32_t lo = 0, hi = n;
    while (lo < hi) { uint32_t m = lo + ((hi - lo) >> 1); if (pr[m].sk < key) lo = m + 1; else hi = m; }
    return lo;
}

// 帯の A（s の和が [u, u+beta+jw)）を取り出して表に入れる。数が上限を超えたら -1。次の帯の読み始め（s の和 >= u+beta）を cn に入れる
static int64_t build_band(thr_t *T, uint64_t u, uint64_t beta, double *ppmin, double *ppmax) {
    uint32_t n = 0, ovn = 0, *tab = T->tab; uint64_t *ovf = T->ovf, *tid = T->tid; double pmin = 1e9, pmax = 0;
    memset(tab, 0, sizeof(uint32_t) * 2 * (size_t)nb);
    for (int k = 0; k < nW; k++) {
        uint32_t p = T->cs[k]; int64_t lim1 = (int64_t)(u + beta) - (int64_t)Ws[k], lim2 = lim1 + jw; double wk = Wp3[k];
        uint64_t key1 = lim1 > 0 ? (uint64_t)lim1 << 28 : 0, key2 = lim2 > 0 ? (uint64_t)lim2 << 28 : 0;
        for (int part = 0; part < 2; part++) {
            uint64_t key = part ? key2 : key1;
            while (p < nPW && pairsW[p].sk < key) {
                uint32_t j = (uint32_t)(pairsW[p].sk & IM);
                if (j <= (uint32_t)k) {
                    if (n >= cap) return -1;
                    double P = pairsW[p].p3 * wk; uint32_t i = (uint32_t)((pairsW[p].sk >> IB) & IM), q = (uint32_t)((P - 1.0) * scale), val = (((q >> 1) & gmask) << xb) | (n + 1), *bk = tab + 2 * (size_t)(q >> tb), v = bk[0];
                    if (v == 0) bk[0] = val;
                    else if (v & FLAG) { ovf[ovn] = ((uint64_t)(v & END) << 32) | val; bk[0] = FLAG | ovn; ovn++; }
                    else if (bk[1] == 0) bk[1] = val;
                    else { ovf[ovn] = ((uint64_t)END << 32) | v; ovf[ovn + 1] = ((uint64_t)ovn << 32) | bk[1]; ovf[ovn + 2] = ((uint64_t)(ovn + 1) << 32) | val; bk[0] = FLAG | (ovn + 2); bk[1] = 0; ovn += 3; }
                    tid[n] = (uint64_t)i | ((uint64_t)j << ib) | ((uint64_t)k << (2 * ib));
                    if (P < pmin) pmin = P; if (P > pmax) pmax = P; n++;
                }
                p++;
            }
            if (!part) T->cn[k] = p;
        }
    }
    *ppmin = pmin; *ppmax = pmax; return n;
}
static void verify(thr_t *T, uint64_t id, uint32_t p, int z, uint64_t Sp) {
    int i = (int)(id & imask), j = (int)((id >> ib) & imask), k = (int)((id >> (2 * ib)) & imask), a = (int)((pairsM[p].sk >> IB) & IM), b = (int)(pairsM[p].sk & IM);
    double P3 = (W[i].P[2] * W[j].P[2] * W[k].P[2]) * (M[a].P[2] * M[b].P[2] * Z[z].P[2]), f3 = P3 - floor(P3);
    if (fabs(f3 - c3) > h3 + slack) return;
    T->nl3++;
    double P1 = (W[i].P[0] * W[j].P[0] * W[k].P[0]) * (M[a].P[0] * M[b].P[0] * Z[z].P[0]), f1 = P1 - floor(P1);
    if (!(fabs(f1 - c1) <= h1 + slack || fabs(f1 - c1b) <= h1b + slack)) return;
    T->nc1++;
    double P2 = (W[i].P[1] * W[j].P[1] * W[k].P[1]) * (M[a].P[1] * M[b].P[1] * Z[z].P[1]), P4 = (W[i].P[3] * W[j].P[3] * W[k].P[3]) * (M[a].P[3] * M[b].P[3] * Z[z].P[3]);
    uint64_t Jsum = (uint64_t)floor(P1) + (uint64_t)floor(P2) + (uint64_t)floor(P3) + (uint64_t)floor(P4);
    if (W[i].s + W[j].s + W[k].s + M[a].s + M[b].s + Z[z].s + Jsum != Sp) return;
    T->nS++;
    int b0b1 = (int)((W[i].t + W[j].t + W[k].t + M[a].t + M[b].t + Z[z].t) & 0xFFFF);
    if (b0b1 < Tlo || b0b1 > Thi) return;
    T->ncand++;
    #pragma omp critical
    { printf("CAND %d %d %d %d %d %d\n", i, j, k, a, b, z); fflush(stdout); }
}
static void __attribute__((noinline)) slow_probe(thr_t *T, uint64_t ql, uint64_t qh, uint32_t p, int z, uint64_t Sp) {
    if (qh < ql) return;
    T->nslow++;
    uint64_t gl = ql >> 1, gh = qh >> 1;
    for (uint64_t bb = ql >> tb; bb <= (qh >> tb); bb++) {
        uint32_t v0 = T->tab[2 * bb], v1 = T->tab[2 * bb + 1];
        if (v0 & FLAG) {
            for (uint32_t c = v0 & END; c != END; c = (uint32_t)(T->ovf[c] >> 32)) {
                uint32_t val = (uint32_t)T->ovf[c]; uint64_t q = (bb << tg) | ((val >> xb) & gmask);
                if (q >= gl && q <= gh) { T->nhit++; verify(T, T->tid[(val & xmask) - 1], p, z, Sp); }
            }
        } else {
            if (v0) { uint64_t q = (bb << tg) | ((v0 >> xb) & gmask); if (q >= gl && q <= gh) { T->nhit++; verify(T, T->tid[(v0 & xmask) - 1], p, z, Sp); } }
            if (v1) { uint64_t q = (bb << tg) | ((v1 >> xb) & gmask); if (q >= gl && q <= gh) { T->nhit++; verify(T, T->tid[(v1 & xmask) - 1], p, z, Sp); } }
        }
    }
}
static inline void move_down(const pair_t *pr, uint32_t *pos, int64_t sval) {           // *pos を「s >= sval になる最初の位置」まで下げる（近くから順に。遠いときは幅を倍々にして探す）
    uint32_t x = *pos; if (sval <= 0) { *pos = 0; return; }
    uint64_t key = (uint64_t)sval << 28; int steps = 0;
    while (x > 0 && pr[x - 1].sk >= key) {
        x--;
        if (++steps > 256) {
            uint32_t hi = x, w = 256;
            while (hi >= w && pr[hi - w].sk >= key) { hi -= w; w <<= 1; }
            uint32_t lo = hi >= w ? hi - w : 0;                 // pr[lo] < key（lo = 0 のときは不明）、pr[hi] >= key
            while (lo < hi) { uint32_t m = lo + ((hi - lo) >> 1); if (pr[m].sk < key) lo = m + 1; else hi = m; }
            x = lo; break;
        }
    }
    *pos = x;
}
// 帯 [u, u+beta) に割り当たる B を順に取り出して表を引く
static void bloop(thr_t *T, uint64_t Sp, uint64_t u, uint64_t beta, double pmin, double pmax) {
    const uint32_t *tab = T->tab; uint64_t nB = 0, nprobe = 0;
#ifdef __AVX512DQ__
    const __m512d vmax = _mm512_set1_pd(maxq), vzero = _mm512_setzero_pd(), viota = _mm512_set_pd(7, 6, 5, 4, 3, 2, 1, 0);
    const __m512i vgmask = _mm512_set1_epi64(gmask), vflag = _mm512_set1_epi64(FLAG);
#endif
    double x0a[BLK], dxa[BLK], cna[BLK];
    for (int z = 0; z < nZ; z++) {
        uint32_t hi = T->bz[z]; move_down(pairsM, &T->bz[z], (int64_t)Sp - jmax - (int64_t)u - (int64_t)beta + 1 - (int64_t)Zs[z]);
        double zp = Zp3[z];
        for (uint32_t p0 = T->bz[z]; p0 < hi; p0 += BLK) {
            int m = hi - p0 < BLK ? (int)(hi - p0) : BLK; const pair_t *pm = pairsM + p0;
            for (int q = 0; q < m; q++) {                    // まず B ごとの値（J3 の範囲、目標の初めと刻み）をまとめて作る
                double P3B = pm[q].p3 * zp, inv = 1.0 / P3B, jl = ceil(P3B * pmin - c3 - 1e-7), jh = floor(P3B * pmax - c3 + 1e-7);
                if (jl < 1.0) jl = 1.0;
                cna[q] = jh - jl + 1.0; dxa[q] = inv * scale; x0a[q] = ((jl + c3) * inv - 1.0) * scale;
            }
            nB += (uint64_t)m;
            for (int q = 0; q < m; q++) {
                int cnt = (int)cna[q]; if (cnt <= 0) continue; nprobe += (uint64_t)cnt;
                double dx = dxa[q], x0 = x0a[q], tolq = tol * dx + margin; uint32_t p = p0 + (uint32_t)q;
#ifdef __AVX512DQ__
                const __m512d vtol = _mm512_set1_pd(tolq), vstep = _mm512_set1_pd(8.0 * dx); __m512d xv = _mm512_add_pd(_mm512_set1_pd(x0), _mm512_mul_pd(viota, _mm512_set1_pd(dx)));
                for (int n = 0; n < cnt; n += 8) {
                    __m512d xl = _mm512_min_pd(_mm512_max_pd(_mm512_sub_pd(xv, vtol), vzero), vmax), xh = _mm512_min_pd(_mm512_max_pd(_mm512_add_pd(xv, vtol), vzero), vmax);
                    __m512i ql = _mm512_cvttpd_epu64(xl), qh = _mm512_cvttpd_epu64(xh);
                    __m512i v = _mm512_i64gather_epi64(_mm512_srli_epi64(ql, tb), (const void *)tab, 8), gl = _mm512_srli_epi64(ql, 1);
                    __m512i tl = _mm512_and_epi64(gl, vgmask), span = _mm512_sub_epi64(_mm512_srli_epi64(qh, 1), gl);
                    __m512i tag0 = _mm512_and_epi64(_mm512_srli_epi64(v, xb), vgmask), tag1 = _mm512_and_epi64(_mm512_srli_epi64(v, 32 + xb), vgmask);
                    __mmask8 mk = _mm512_cmple_epu64_mask(_mm512_sub_epi64(tag0, tl), span) | _mm512_cmple_epu64_mask(_mm512_sub_epi64(tag1, tl), span) | _mm512_test_epi64_mask(v, vflag) | _mm512_cmpgt_epu64_mask(_mm512_add_epi64(tl, span), vgmask);
                    if (mk) {
                        uint64_t aql[8], aqh[8]; _mm512_storeu_si512((void *)aql, ql); _mm512_storeu_si512((void *)aqh, qh);
                        unsigned mm = mk; while (mm) { int l = __builtin_ctz(mm); mm &= mm - 1; if (n + l < cnt) slow_probe(T, aql[l], aqh[l], p, z, Sp); }
                    }
                    xv = _mm512_add_pd(xv, vstep);
                }
#else
                for (int n = 0; n < cnt; n++) {
                    double x = x0 + n * dx, xl = x - tolq, xh = x + tolq;
                    if (xl < 0) xl = 0; if (xl > maxq) xl = maxq; if (xh < 0) xh = 0; if (xh > maxq) xh = maxq;
                    uint64_t ql = (uint64_t)xl, qh = (uint64_t)xh; uint32_t v0 = tab[2 * (ql >> tb)], v1 = tab[2 * (ql >> tb) + 1], tl = (uint32_t)(ql >> 1) & gmask, span = (uint32_t)((qh >> 1) - (ql >> 1));
                    if (((((v0 >> xb) & gmask) - tl) <= span) | ((((v1 >> xb) & gmask) - tl) <= span) | (v0 >> 31) | ((tl + span) > gmask)) slow_probe(T, ql, qh, p, z, Sp);
                }
#endif
            }
        }
    }
    T->nB += nB; T->nprobe += nprobe;
}
static void do_range(thr_t *T, const range_t *R) {
    uint64_t Sp = R->Sp, U0 = R->U0, U1 = R->U1, est = 0;
    for (int k = 0; k < nW; k++) { T->cs[k] = lbp(pairsW, nPW, (int64_t)U0 - (int64_t)Ws[k]); est += lbp(pairsW, nPW, (int64_t)U1 + jw - (int64_t)Ws[k]) - T->cs[k]; }
    est /= 3;
    double rho = (double)est / (double)(U1 - U0 + jw);
    for (int z = 0; z < nZ; z++) T->bz[z] = lbp(pairsM, nPM, (int64_t)Sp - jmax - (int64_t)U0 - (int64_t)Zs[z] + 1);
    uint64_t u = U0;
    while (u < U1) {
        // 帯の幅: ふつうは表に ntarget 組が入る幅。混んだところでは重なり（jw）と同じ幅までは保ち（表の組数は増える）、表の上限の 6 割を超えるときだけ狭める
        double wd = rho > 0 ? ntarget / rho : 1e18; if (wd < 2.0 * jw + 2) wd = 2.0 * jw + 2; if (rho > 0 && wd * rho > 0.6 * cap) wd = 0.6 * cap / rho;
        uint64_t beta = wd > jw + 2 ? (wd < 1e18 ? (uint64_t)(wd - jw - 1) : U1 - u) : 1; if (beta > U1 - u) beta = U1 - u;
        double pmin, pmax, t0 = omp_get_wtime(); int64_t n;
        for (;;) {
            n = build_band(T, u, beta, &pmin, &pmax);
            if (n >= 0) break;
            if (beta <= 1) { fprintf(stderr, "band overflow at the smallest width (raise the table bits)\n"); exit(1); }
            beta = (beta + 1) / 2; T->nsplit++;
        }
        double t1 = omp_get_wtime();
        if (n > 0) { bloop(T, Sp, u, beta, pmin, pmax); T->nA += (uint64_t)n; T->nband++; }
        else for (int z = 0; z < nZ; z++) move_down(pairsM, &T->bz[z], (int64_t)Sp - jmax - (int64_t)u - (int64_t)beta + 1 - (int64_t)Zs[z]);
        T->tbuild += t1 - t0; T->tprobe += omp_get_wtime() - t1;
        uint32_t *tmp = T->cs; T->cs = T->cn; T->cn = tmp; u += beta;
        double r2 = (double)n / (double)(beta + jw); rho = n > 0 ? (r2 > 4 * rho ? 4 * rho : r2) : rho / 4;
    }
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: k6band spec\n"); return 1; }
    FILE *f = fopen(argv[1], "r"); if (!f) { perror(argv[1]); return 1; }
    char line[1024], path[900]; unsigned ts; int n; static int todo_a[4096], todo_b[4096]; int ntodo = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'T') { sscanf(line + 1, "%x %lf %lf %lf %lf %lf %lf %d %d %lf", &ts, &c1, &h1, &c1b, &h1b, &c3, &h3, &Tlo, &Thi, &slack); TS = ts; }
        else if (line[0] == 'J') sscanf(line + 1, "%d %d", &jmin, &jmax);
        else if (line[0] == 'W') { sscanf(line + 1, "%d %899s", &n, path); nW = n; W = load(path, n); }
        else if (line[0] == 'M') { sscanf(line + 1, "%d %899s", &n, path); nM = n; M = load(path, n); }
        else if (line[0] == 'Z') { sscanf(line + 1, "%d %899s", &n, path); nZ = n; Z = load(path, n); }
        else if (line[0] == 'N') sscanf(line + 1, "%d %d %d", &ntarget, &dbits, &nrange_req);
        else if (line[0] == 'L') sscanf(line + 1, "%lf", &tlimit);
        else if (line[0] == 'R' && ntodo < 4096) { sscanf(line + 1, "%d %d", &todo_a[ntodo], &todo_b[ntodo]); ntodo++; }
    }
    fclose(f);
    if (nW > (1 << IB) || nM > (1 << IB)) { fprintf(stderr, "too many units (max %d)\n", 1 << IB); return 1; }
    if (dbits < 10 || dbits > 24) { fprintf(stderr, "table bits must be 10..24\n"); return 1; }
    double T0 = omp_get_wtime();
    jw = jmax - jmin; ib = 1; while ((1 << ib) < nW) ib++; imask = (1ULL << ib) - 1;
    nb = 1u << dbits; tb = 32 - dbits; tg = tb - 1; gmask = (1u << tg) - 1; xb = dbits; xmask = (1u << xb) - 1; cap = xmask - 1; if ((uint32_t)ntarget > cap / 2) ntarget = (int)(cap / 2);
    scale = 4294967296.0 / 7.000001; tol = h3 + slack;
    #pragma omp parallel sections
    {
        #pragma omp section
        pairsW = make_pairs(W, nW, &nPW);
        #pragma omp section
        pairsM = make_pairs(M, nM, &nPM);
    }
    Ws = malloc(8 * (size_t)nW); Wp3 = malloc(8 * (size_t)nW); for (int i = 0; i < nW; i++) { Ws[i] = W[i].s; Wp3[i] = W[i].P[2]; }
    Zs = malloc(8 * (size_t)nZ); Zp3 = malloc(8 * (size_t)nZ); for (int i = 0; i < nZ; i++) { Zs[i] = Z[i].s; Zp3[i] = Z[i].P[2]; }
    // 区間: 桁あふれの回数 w ごとに、u（= A の s の和の下限）のうち相手がありうる範囲を、仕事がほぼ等しくなるように分ける（乱数の標本の分位点）
    uint64_t wmin = ~0ULL, wmax = 0, mmin = ~0ULL, mmax = 0, zmin = ~0ULL, zmax = 0;
    for (int i = 0; i < nW; i++) { if (W[i].s < wmin) wmin = W[i].s; if (W[i].s > wmax) wmax = W[i].s; }
    for (int i = 0; i < nM; i++) { if (M[i].s < mmin) mmin = M[i].s; if (M[i].s > mmax) mmax = M[i].s; }
    for (int i = 0; i < nZ; i++) { if (Z[i].s < zmin) zmin = Z[i].s; if (Z[i].s > zmax) zmax = Z[i].s; }
    int64_t Amin = (int64_t)(3 * wmin), Amax = (int64_t)(3 * wmax), Bmin = (int64_t)(2 * mmin + zmin), Bmax = (int64_t)(2 * mmax + zmax);
    range_t *rg = malloc(sizeof(range_t) * ((size_t)nrange_req + 64) * 4); int nrg = 0;
    int ns = 300000; uint64_t *smp = malloc(8 * ((size_t)ns * 4 + 16)); unsigned long long r = 88172645463325252ULL;
    for (int w = 0; w < 16; w++) {
        int64_t Sp = (int64_t)TS + ((int64_t)w << 32), lo = Sp - jmax - Bmax, hi = Sp - jmax - Bmin + 1;
        if (lo < Amin - jw) lo = Amin - jw; if (hi > Amax + 1) hi = Amax + 1; if (lo < 0) lo = 0;
        if (hi <= lo) continue;
        int m = 0;
        for (int i = 0; i < ns; i++) {
            uint64_t v = 0; for (int q = 0; q < 3; q++) { r ^= r << 13; r ^= r >> 7; r ^= r << 17; v += Ws[r % (unsigned)nW]; }
            if ((int64_t)v >= lo && (int64_t)v < hi) smp[m++] = v;
            uint64_t sb = 0; for (int q = 0; q < 2; q++) { r ^= r << 13; r ^= r >> 7; r ^= r << 17; sb += M[r % (unsigned)nM].s; }
            r ^= r << 13; r ^= r >> 7; r ^= r << 17; sb += Zs[r % (unsigned)nZ];
            int64_t x = Sp - jmax - (int64_t)sb;
            if (x >= lo && x < hi) { smp[m++] = (uint64_t)x; smp[m++] = (uint64_t)x; smp[m++] = (uint64_t)x; }
        }
        int nr = (int)((double)nrange_req * m / (4.0 * ns)) + 1; if (nr > nrange_req) nr = nrange_req;
        qsort(smp, m, 8, cmp_u64);
        uint64_t prev = (uint64_t)lo;
        for (int q = 1; q <= nr; q++) {
            uint64_t u = q == nr ? (uint64_t)hi : smp[(size_t)((double)m * q / nr)];
            if (u > prev) { rg[nrg].Sp = (uint64_t)Sp; rg[nrg].U0 = prev; rg[nrg].U1 = u; nrg++; prev = u; }
        }
    }
    printf("NRANGE %d\n", nrg); fflush(stdout);
    fprintf(stderr, "units W %d M %d Z %d, pairs %u %u, index bits %d tag bits %d, band %d entries (max %u), table 2^%d x 2 slots, ranges %d, threads %d, setup %.1fs%s\n", nW, nM, nZ, nPW, nPM, ib, tg, ntarget, cap, dbits, nrg, omp_get_max_threads(), omp_get_wtime() - T0,
#ifdef __AVX512DQ__
            ", avx512"
#else
            ", scalar"
#endif
            );
    char *todo = calloc(nrg + 1, 1); int any = 0;
    for (int q = 0; q < ntodo; q++) for (int x = todo_a[q]; x <= todo_b[q] && x < nrg; x++) if (x >= 0) { todo[x] = 1; any = 1; }
    if (!ntodo) { memset(todo, 1, nrg); any = 1; }
    uint64_t tA = 0, tB = 0, tprobe = 0, thit = 0, tl3 = 0, tc1 = 0, tS = 0, tcand = 0, tband = 0, tsplit = 0, tslow = 0; int ndone = 0; double sbuild = 0, sprobe = 0;
    if (any) {
        #pragma omp parallel reduction(+:tA,tB,tprobe,thit,tl3,tc1,tS,tcand,tband,tsplit,tslow,ndone,sbuild,sprobe)
        {
            thr_t T; memset(&T, 0, sizeof T);
            T.tab = big_alloc(8 * (size_t)nb); T.tid = malloc(8 * ((size_t)cap + 8)); T.ovf = malloc(8 * (2 * (size_t)cap + 8));
            T.cs = malloc(4 * (size_t)nW); T.cn = malloc(4 * (size_t)nW); T.bz = malloc(4 * (size_t)nZ);
            #pragma omp for schedule(dynamic, 1)
            for (int q = 0; q < nrg; q++) {
                if (!todo[q]) continue;
                if (tlimit > 0 && omp_get_wtime() - T0 > tlimit) continue;
                double t0 = omp_get_wtime(); uint64_t a0 = T.nA, b0 = T.nB, p0 = T.nprobe;
                do_range(&T, &rg[q]); ndone++;
                #pragma omp critical
                { printf("DONE %d %llu %llu %llu %.2f\n", q, (unsigned long long)(T.nA - a0), (unsigned long long)(T.nB - b0), (unsigned long long)(T.nprobe - p0), omp_get_wtime() - t0); fflush(stdout); }
            }
            tA += T.nA; tB += T.nB; tprobe += T.nprobe; thit += T.nhit; tl3 += T.nl3; tc1 += T.nc1; tS += T.nS; tcand += T.ncand; tband += T.nband; tsplit += T.nsplit; tslow += T.nslow; sbuild += T.tbuild; sprobe += T.tprobe;
        }
    }
    printf("STAT ranges %d done %d bands %llu splits %llu A %llu B %llu probe %llu slow %llu hit %llu l3 %llu c1 %llu S %llu cand %llu tbuild %.1f tprobe %.1f sec %.1f\n", nrg, ndone, (unsigned long long)tband, (unsigned long long)tsplit,
           (unsigned long long)tA, (unsigned long long)tB, (unsigned long long)tprobe, (unsigned long long)tslow, (unsigned long long)thit, (unsigned long long)tl3, (unsigned long long)tc1, (unsigned long long)tS, (unsigned long long)tcand, sbuild, sprobe, omp_get_wtime() - T0);
    return 0;
}
