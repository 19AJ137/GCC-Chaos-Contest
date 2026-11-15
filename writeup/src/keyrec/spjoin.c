// spjoin: join the solutions of the f3 system (X: multisets of N3) and of the f1 system (Y: multisets of N1) through the
// chunk table of a character set: a pair (N3, N1) is valid when some chunk of the table has those two lattice values.
// build: gcc -O2 -march=native -o spjoin spjoin.c
// usage: spjoin n table.bin X.txt Y.txt   (X lines: tag J3 N.. ; Y lines: tag SJ J1 br N.. ; both sorted by tag)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef struct { uint32_t n3, n1, c; } Ent;
typedef struct { int tag, J; uint32_t N[6]; } XS;
typedef struct { int tag, SJ, J, br; uint32_t N[6]; } YS;
typedef struct { uint32_t y; int id; } YI;
static Ent *tab; static size_t ntab; static int n;
static uint32_t lo3, hi3; static size_t *off3;
static uint64_t npairs = 0, nout = 0;

static int cmp_yi(const void *a, const void *b) { const YI *x = a, *y = b; return x->y < y->y ? -1 : x->y > y->y ? 1 : (x->id - y->id); }

// chunks with (n3, n1): entries of fiber(n3) are sorted by n1
static inline size_t find(uint32_t x, uint32_t y, size_t *end) {
    if (x < lo3 || x >= hi3) { *end = 0; return 0; }
    size_t a = off3[x - lo3], b = off3[x - lo3 + 1], e = b;
    while (a < b) { size_t m = (a + b) >> 1; if (tab[m].n1 < y) a = m + 1; else b = m; }
    size_t s = a; while (a < e && tab[a].n1 == y) a++;
    *end = a; return s;
}

static void rec(const XS *x, const YS *y, int i, unsigned used, uint32_t *acc) {
    if (i == n) {
        printf("%d %d %d %d %d", x->tag, x->J, y->J, y->br, y->SJ);
        for (int k = 0; k < n; k++) printf(" %u", acc[k]);
        printf("\n"); nout++; return;
    }
    for (int j = 0; j < n; j++) {
        if (used >> j & 1) continue;
        if (j > 0 && y->N[j] == y->N[j - 1] && !(used >> (j - 1) & 1)) continue;
        size_t e, s = find(x->N[i], y->N[j], &e);
        for (; s < e; s++) { acc[i] = tab[s].c; rec(x, y, i + 1, used | 1u << j, acc); }
    }
}

int main(int argc, char **argv) {
    if (argc < 5) return 1;
    n = atoi(argv[1]);
    FILE *f = fopen(argv[2], "rb"); fseek(f, 0, SEEK_END); ntab = ftell(f) / sizeof(Ent); fseek(f, 0, SEEK_SET);
    tab = malloc(ntab * sizeof(Ent)); if (fread(tab, sizeof(Ent), ntab, f) != ntab) return 1; fclose(f);
    lo3 = tab[0].n3; hi3 = tab[ntab - 1].n3 + 1;
    off3 = calloc(hi3 - lo3 + 2, sizeof(size_t));
    for (size_t i = 0; i < ntab; i++) off3[tab[i].n3 - lo3 + 1]++;
    for (uint32_t v = 0; v < hi3 - lo3 + 1; v++) off3[v + 1] += off3[v];
    size_t nx = 0, ny = 0, cx = 1 << 20, cy = 1 << 20; XS *X = malloc(cx * sizeof(XS)); YS *Y = malloc(cy * sizeof(YS));
    f = fopen(argv[3], "r");
    for (;;) { XS s; if (fscanf(f, "%d %d", &s.tag, &s.J) != 2) break; for (int k = 0; k < n; k++) if (fscanf(f, "%u", &s.N[k]) != 1) return 1;
        if (nx == cx) { cx *= 2; X = realloc(X, cx * sizeof(XS)); } X[nx++] = s; }
    fclose(f); f = fopen(argv[4], "r");
    for (;;) { YS s; if (fscanf(f, "%d %d %d %d", &s.tag, &s.SJ, &s.J, &s.br) != 4) break; for (int k = 0; k < n; k++) if (fscanf(f, "%u", &s.N[k]) != 1) return 1;
        if (ny == cy) { cy *= 2; Y = realloc(Y, cy * sizeof(YS)); } Y[ny++] = s; }
    fclose(f);
    size_t ix = 0, iy = 0; YI *idx = NULL; size_t cidx = 0;
    while (ix < nx && iy < ny) {
        if (X[ix].tag < Y[iy].tag) { ix++; continue; }
        if (Y[iy].tag < X[ix].tag) { iy++; continue; }
        int tag = X[ix].tag; size_t ex = ix, ey = iy;
        while (ex < nx && X[ex].tag == tag) ex++;
        while (ey < ny && Y[ey].tag == tag) ey++;
        size_t m = (ey - iy) * n; if (m > cidx) { cidx = m; idx = realloc(idx, cidx * sizeof(YI)); }
        size_t k = 0;
        for (size_t j = iy; j < ey; j++) for (int e = 0; e < n; e++) { if (e > 0 && Y[j].N[e] == Y[j].N[e - 1]) continue; idx[k].y = Y[j].N[e]; idx[k].id = (int)j; k++; }
        qsort(idx, k, sizeof(YI), cmp_yi);
        for (size_t i = ix; i < ex; i++) {
            const XS *x = &X[i]; uint32_t x0 = x->N[0];
            if (x0 < lo3 || x0 >= hi3) continue;
            for (size_t e = off3[x0 - lo3]; e < off3[x0 - lo3 + 1]; e++) {
                if (e > off3[x0 - lo3] && tab[e].n1 == tab[e - 1].n1) continue;      // same (n3,n1): handled once, rec() expands chunks
                uint32_t yv = tab[e].n1; size_t a = 0, b = k;
                while (a < b) { size_t mm = (a + b) >> 1; if (idx[mm].y < yv) a = mm + 1; else b = mm; }
                for (; a < k && idx[a].y == yv; a++) {
                    const YS *y = &Y[idx[a].id]; npairs++;
                    if (y->SJ - y->J - x->J < 2) continue;
                    uint32_t acc[6]; rec(x, y, 0, 0, acc);
                }
            }
        }
        ix = ex; iy = ey;
    }
    fprintf(stderr, "X %zu Y %zu pairs %llu out %llu\n", nx, ny, (unsigned long long)npairs, (unsigned long long)nout);
    return 0;
}
