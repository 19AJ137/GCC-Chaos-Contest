# 短い等価鍵の計画: 固定する塊（小数部が 4 つとも大きいもの）を決めたとき、残り 3 塊（乱数）で b0b1 の条件を満たせる割合 τ を測る。
#   python3 short/vshort_plan.py
import sys, os, math, struct, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from gks import K_PI, K_E, C01, chunk_vals, map_of_i2
S, Tlo, Thi = 0x3a9d3481, 0xab0f, 0xabc5
RHO = K_E / K_PI
def clr_np(a): return (a.view(np.uint64) & np.uint64(0xFFFFFFFFFFFFFF00)).view(np.float64)
def vals_np(c):
    c = c.astype(np.float64)
    a1 = clr_np(c * K_PI); i1 = np.floor(a1); a2 = clr_np(i1 * K_PI); i2 = np.floor(a2)
    a3 = clr_np(c * K_E); i3 = np.floor(a3); a4 = clr_np(i3 * K_E); i4 = np.floor(a4)
    return a1 - i1, a2 - i2, a3 - i3, a4 - i4
def tau(fixed, nsamp=400000, a=0x3243F6A8, b=1 << 32, seed=1, w=1):
    """fixed: 固定する塊の値のリスト。未知 3 塊を [a, b) の一様乱数にしたとき、S が合ったとして b0b1 が範囲に入る割合"""
    rng = np.random.default_rng(seed); n = len(fixed) + 3
    PK = [1.0] * 4; FK = [0.0] * 4; R = 0; Q = 0
    for c in fixed:
        f, i, s, t = chunk_vals(c)
        for j in range(4): PK[j] *= 1.0 + f[j]; FK[j] += f[j]
        R += (s - int(i[0])) & 0xFFFFFFFF; Q += (t - int(i[2])) & 0xFFFFFFFF
    R += 9; Q += 6
    P = [np.full(nsamp, PK[j]) for j in range(4)]; F1 = np.full(nsamp, FK[0]); F3 = np.full(nsamp, FK[2])
    for k in range(3):
        f = vals_np(rng.integers(a, b, nsamp, dtype=np.int64))
        for j in range(4): P[j] = P[j] * (1.0 + f[j])
        F1 += f[0]; F3 += f[2]
    SJ = sum(np.floor(p) for p in P)
    Y = R + SJ - F1; Z = Q - F3
    Sw = S + (w << 32)
    v = np.mod(RHO * (Sw - Y) + Z, 65536.0)
    ok = (v >= Tlo - 0.01) & (v <= Thi + 0.01)
    return float(ok.mean()), float(SJ.mean()), PK
def last_table(nbytes):
    out = []
    for c in range(1 << (8 * nbytes)):
        if nbytes > 1 and c < (1 << (8 * (nbytes - 1))): continue      # 最上位バイトが 0 のものは、短い端数と同じ値
        f, i, s, t = chunk_vals(c)
        out.append((c, f, map_of_i2(i[1]), (s - int(i[0])) & 0xFFFFFFFF, (t - int(i[2])) & 0xFFFFFFFF))
    return out

if __name__ == '__main__':
    hi = []
    for line in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hi90.txt')):
        p = line.split(); hi.append((int(p[0], 16), tuple(float(x) for x in p[1:5]), int(p[5])))
    key = lambda x: sum(math.log2(1 + v) for v in x[1])
    hi.sort(key=key, reverse=True)
    print('小数部が 4 つとも 0.90 以上の塊: %d 個。うち 0.95 以上 %d、0.97 以上 %d、0.98 以上 %d、0.99 以上 %d' % (
        len(hi), sum(min(x[1]) >= 0.95 for x in hi), sum(min(x[1]) >= 0.97 for x in hi), sum(min(x[1]) >= 0.98 for x in hi), sum(min(x[1]) >= 0.99 for x in hi)))
    print('上位 5 個:'); 
    for c, f, m in hi[:5]: print('   %08x  f = %.4f %.4f %.4f %.4f  写像 %d' % (c, f[0], f[1], f[2], f[3], m))
    # 端数（最後の塊）1 バイト: 写像 1 のもの
    t1 = last_table(1)
    print('端数 1 バイトで写像 1 になる値: %d 個。Π(1+f) の大きい順に 8 個' % sum(1 for x in t1 if x[2] == 1))
    for c, f, m, r, q in sorted((x for x in t1 if x[2] == 1), key=lambda x: -sum(math.log2(1 + v) for v in x[1]))[:8]:
        print('   %02x  f = %.4f %.4f %.4f %.4f  (1+f) の積 %.3f  r=%d q=%d' % (c, f[0], f[1], f[2], f[3], math.prod(1 + v for v in f), r, q))
    t2 = last_table(2)
    print('端数 2 バイトで写像 1: %d 個。上位 5 個' % sum(1 for x in t2 if x[2] == 1))
    for c, f, m, r, q in sorted((x for x in t2 if x[2] == 1), key=lambda x: -sum(math.log2(1 + v) for v in x[1]))[:5]:
        print('   %04x（鍵のバイト列 %s）  f = %.4f %.4f %.4f %.4f  積 %.3f  r=%d q=%d' % (c, c.to_bytes(2, 'little').hex(), f[0], f[1], f[2], f[3], math.prod(1 + v for v in f), r, q))
    t3 = [x for x in hi if x[0] < (1 << 24) and x[0] >= (1 << 16) and x[2] == 1]
    print('端数 3 バイトで写像 1・小数部 0.90 以上: %d 個。上位 5 個' % len(t3))
    for c, f, m in t3[:5]: print('   %06x（鍵のバイト列 %s）  f = %.4f %.4f %.4f %.4f  積 %.3f' % (c, c.to_bytes(3, 'little').hex(), f[0], f[1], f[2], f[3], math.prod(1 + v for v in f)))
    # τ の見積もり
    big = [x for x in hi if x[0] >= 0x33000000]
    top = big[:3]
    print('24 バイト（固定 3 塊 = 上位 3 個 %s）: τ = %.4f、ΣJ の平均 %.1f' % (' '.join('%08x' % x[0] for x in top), *tau([x[0] for x in top])[:2]))
    for q in (30, 300, 3000, 30000):
        sel = big[q:q + 3]; r = tau([x[0] for x in sel]); print('   上位 %d 番目あたりの 3 個: τ = %.4f、ΣJ の平均 %.1f、最小の小数部 %.3f' % (q, r[0], r[1], min(min(x[1]) for x in sel)))
    for nb, tab in ((3, [(c, f, m, 3, 2) for c, f, m in t3]), (2, [x for x in t2 if x[2] == 1]), (1, [x for x in t1 if x[2] == 1])):
        best = sorted(tab, key=lambda x: -sum(math.log2(1 + v) for v in x[1]))[:3]
        for bl in best:
            r = tau([big[0][0], big[1][0], bl[0]]); print('%d バイト（固定 2 塊＋端数 %x）: τ = %.5f、ΣJ の平均 %.1f' % (20 + nb, bl[0], r[0], r[1]))
    # 塊 7 個（25〜28 バイト）
    for q in (0, 3000, 30000, 300000):
        sel = big[q:q + 4]; r = tau([x[0] for x in sel]); print('28 バイト（固定 4 塊、上位 %d 番目あたり）: τ = %.4f、ΣJ の平均 %.1f' % (q, r[0], r[1]))
