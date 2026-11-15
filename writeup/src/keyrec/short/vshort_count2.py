# 短い等価鍵の総数の見積もり（重みつきの乱数。小数部 f1..f4 は塊ごとに独立な一様乱数とみなす）。
#   python3 short/vshort_count2.py
import sys, os, math
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
from gks import K_PI, K_E, chunk_vals, map_of_i2
from vshort_plan import S, Tlo, Thi, RHO
from vshort_count import irwin_hall
CELLS = (2 * 7.276e-12) * 5.82e-11
XS = (S + (1 << 32)) / K_PI / 2.0 ** 32
def tau_is(nfull, lastf=None, lastrq=(0, 0), theta=3.0, nsamp=2000000, seed=5):
    """nfull 個の塊の f を (1+f)^theta に比例する密度から引き、重みをつけて「b0b1 が範囲に入る割合」を出す。
    lastf: 端数の塊の f（形 (nsamp, 4) の配列か、長さ 4 の組）"""
    rng = np.random.default_rng(seed); a = theta + 1.0; c = 2.0 ** a - 1.0
    u = rng.random((nsamp, nfull, 4)); g = (1.0 + u * c) ** (1.0 / a)          # g = 1 + f
    logw = (np.log(c / a) - theta * np.log(g)).sum(axis=(1, 2))                  # 一様 ÷ 傾けた密度
    P = g.prod(axis=1); F1 = (g[:, :, 0] - 1.0).sum(axis=1); F3 = (g[:, :, 2] - 1.0).sum(axis=1)
    R = 3 * nfull + lastrq[0]; Q = 2 * nfull + lastrq[1]
    if lastf is not None:
        lf = np.broadcast_to(np.asarray(lastf, dtype=np.float64), (nsamp, 4)); P = P * (1.0 + lf); F1 = F1 + lf[:, 0]; F3 = F3 + lf[:, 2]
    SJ = np.floor(P).sum(axis=1); Y = R + SJ - F1; Z = Q - F3
    v = np.mod(RHO * (S + (1 << 32) - Y) + Z, 65536.0); ok = (v >= Tlo - 0.01) & (v <= Thi + 0.01)
    w = np.exp(logw); est = float((w * ok).mean()); se = float((w * ok).std() / math.sqrt(nsamp))
    return est, se, int(ok.sum())
if __name__ == '__main__':
    # 確かめ: 塊 7・8 個は、ふつうの乱数（vshort_count.py）で τ = 0.0024・0.14
    for n, th in ((7, 1.0), (8, 0.3)):
        e, se, k = tau_is(n, theta=th); print('確かめ 塊 %d 個: τ = %.4g ± %.2g' % (n, e, se))
    e6, se6, k = tau_is(6, theta=3.5, nsamp=3000000)
    N24 = 2.0 ** 192 * (math.pi * irwin_hall(6, XS) / 2.0 ** 32) * CELLS * e6 / 3
    print('24 バイト（塊 6 個）: τ = %.3g ± %.1g → 総数 約 2^%.1f' % (e6, se6, math.log2(N24)))
    # 21 バイト: 端数 1 バイト（写像 1 の 85 通り）
    tot = 0.0; var = 0.0; rows = []
    for b in range(256):
        f, i, s, t = chunk_vals(b)
        if map_of_i2(i[1]) != 1: continue
        r = (s - int(i[0])) & 0xFFFFFFFF; q = (t - int(i[2])) & 0xFFFFFFFF
        e, se, k = tau_is(5, lastf=f, lastrq=(r, q), theta=5.0, nsamp=400000, seed=b); tot += e; var += se * se; rows.append((e, b))
    rows.sort(reverse=True)
    N21 = 2.0 ** 160 * (math.pi * irwin_hall(5, XS) / 2.0 ** 32) * CELLS * tot
    print('21 バイト: τ の合計 %.3g ± %.1g（上位: %s）→ 総数 約 2^%.1f' % (tot, math.sqrt(var), '  '.join('%02x %.2g' % (b, e) for e, b in rows[:5]), math.log2(N21)))
    # 22 バイト: 端数 2 バイト（写像 1 のもの全部から乱数で選ぶ）
    rng = np.random.default_rng(11)
    for nb, cnt in ((2, None), (3, None)):
        if nb == 2: vals = [c for c in range(256, 65536) if map_of_i2(chunk_vals(c)[1][1]) == 1]
        else:
            cand = rng.integers(65536, 1 << 24, 300000); vals = [int(c) for c in cand if map_of_i2(chunk_vals(int(c))[1][1]) == 1]
        nmap = len(vals) if nb == 2 else int((1 << 24) * len(vals) / 300000)
        pick = rng.choice(len(vals), 1500000); fv = np.array([chunk_vals(c)[0] for c in vals])[pick]
        e, se, k = tau_is(5, lastf=fv, lastrq=(3, 2), theta=4.5, nsamp=1500000, seed=20 + nb)
        N = 2.0 ** 160 * nmap * (math.pi * irwin_hall(5, XS) / 2.0 ** 32) * CELLS * e
        print('%d バイト: 端数で写像 1 になる値 %d 個、τ の平均 %.3g ± %.1g → 総数 約 2^%.1f' % (20 + nb, nmap, e, se, math.log2(N)))
