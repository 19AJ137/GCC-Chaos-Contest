# 乱数の鍵（塊の値は一様）で、S が合ったとして b0b1 が鍵 v の範囲に入る割合 τ を塊の数ごとに測る。あわせて 21 バイトの等価鍵の総数を見積もる。
#   python3 short/vshort_count.py
import sys, os, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
from gks import K_PI, K_E, chunk_vals, map_of_i2
from vshort_plan import vals_np, S, Tlo, Thi, RHO
def tau_rand(nfull, last=None, nsamp=4000000, lo=0, hi=1 << 32, seed=3, w=1, batch=1000000):
    rng = np.random.default_rng(seed); n = nfull + (last is not None); ok = 0; tot = 0; sj = 0.0
    PK = [1.0] * 4; FK = [0.0] * 4; R = 3 * nfull; Q = 2 * nfull
    if last is not None:
        f, i, s, t = chunk_vals(last)
        for j in range(4): PK[j] *= 1.0 + f[j]; FK[j] += f[j]
        R += (s - int(i[0])) & 0xFFFFFFFF; Q += (t - int(i[2])) & 0xFFFFFFFF
    Sw = S + (w << 32)
    while tot < nsamp:
        P = [np.full(batch, PK[j]) for j in range(4)]; F1 = np.full(batch, FK[0]); F3 = np.full(batch, FK[2])
        for k in range(nfull):
            f = vals_np(rng.integers(lo, hi, batch, dtype=np.int64))
            for j in range(4): P[j] = P[j] * (1.0 + f[j])
            F1 += f[0]; F3 += f[2]
        SJ = sum(np.floor(p) for p in P); Y = R + SJ - F1; Z = Q - F3
        v = np.mod(RHO * (Sw - Y) + Z, 65536.0)
        ok += int(((v >= Tlo - 0.01) & (v <= Thi + 0.01)).sum()); tot += batch; sj += float(SJ.sum())
    return ok / tot, sj / tot, ok
def irwin_hall(n, x):
    return sum((-1) ** k * math.comb(n, k) * max(x - k, 0.0) ** (n - 1) for k in range(n + 1)) / math.factorial(n - 1)
if __name__ == '__main__':
    print('塊の数ごとの τ（塊の値は 0x20000000〜0xFFFFFFFF の一様乱数。端数なし）')
    for n in range(6, 12):
        t, sj, ok = tau_rand(n, nsamp=10000000 if n <= 7 else 2000000, lo=0x20000000)
        print('  塊 %2d 個（%d〜%d バイト）: τ = %.3g（%d 個）  ΣJ の平均 %.0f' % (n, 4 * n - 3, 4 * n, t, ok, sj))
    # 21 バイト: 塊 5 個＋端数 1 バイト（写像 1 になる 85 通り）
    cells = (2 * 7.276e-12) * 5.82e-11
    xs = (0x3a9d3481 + (1 << 32)) / K_PI / 2.0 ** 32
    tot = 0.0; rows = []
    for b in range(256):
        f, i, s, t = chunk_vals(b)
        if map_of_i2(i[1]) != 1: continue
        tb, sj, ok = tau_rand(5, last=b, nsamp=3000000); tot += tb; rows.append((tb, b, ok))
    rows.sort(reverse=True)
    print('21 バイト: 端数の値ごとの τ（上位 6）: ' + '  '.join('%02x: %.2g（%d 個）' % (b, tb, ok) for tb, b, ok in rows[:6]))
    dens = irwin_hall(5, xs)
    N = (2.0 ** 160) * (math.pi * dens / 2.0 ** 32) * cells * tot
    print('  Σc ÷ 2^32 = %.3f、塊 5 個の和の密度 %.4f、τ の合計 %.3g、セルの割合 %.3g' % (xs, dens, tot, cells))
    print('  21 バイトの等価鍵の総数の見積もり: 2^%.1f 個' % math.log2(N))
    # 22〜24 バイトの数と、塊 6 個の割合（まれな場合に重みをつけた数え方）は vshort_count2.py
