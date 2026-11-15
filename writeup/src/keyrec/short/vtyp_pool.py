# 打てる 21 バイトの鍵の探索で使う「固定する 2 塊＋端数」の組を、当たりやすさの見積もりの順に並べる。
#   当たりやすさ ∝ τ（残る 3 塊で b0b1 の条件を満たせる割合。固定した塊の (1+f) の積の和から見積もる）
#                × 解く 3 塊の和になる組の密度 × 解けた鍵が打てる文字になる割合（どちらも和の位置 x の関数）
#   python3 short/vtyp_pool.py        見積もりの確かめ（τ と点数の関係）と、組の数
import sys, os, math, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import numpy as np
from vpart import textlike
from gks import chunk_vals
from vshort import load_hi, SIGMA, A_U, B_U
b4 = lambda c: c.to_bytes(4, 'little'); D = B_U - A_U
LASTS = (0x51, 0x6a)
def tau_est(score): return 0.0100 * np.exp(0.685 * (score - 28.14))      # 実測（下の確かめ）に合わせた式。点数 27.1 で 0.5%、28.4 で 1.2%
def ih3(x): return np.where(x <= 1.0, x * x / 2, np.where(x <= 2.0, (-2 * x * x + 6 * x - 3) / 2, (3 - x) ** 2 / 2))
def typ(x): return np.where(x < 1.95, 0.058, np.where(x < 2.15, 0.052, 0.038))
def chunks_for_pool(hi, fmin=0.93, vmin=0xA0000000):
    return [x for x in hi if min(x[1]) >= fmin and x[0] >= vmin and textlike(b4(x[0])) >= 2]
def ranked_pool(hi, done=(), fmin=0.93, xlo=1.55, xhi=2.25, nkeep=1500000, log=print):
    """戻り値: 既知のバイト列（固定する 2 塊を値の小さい順に並べ、端数を付けたもの）のリスト。当たりやすい順"""
    ch = chunks_for_pool(hi, fmin); N = len(ch)
    C = np.array([x[0] for x in ch], dtype=np.int64); G = 1.0 + np.array([x[1] for x in ch])
    GL = {l: 1.0 + np.array(chunk_vals(l)[0]) for l in LASTS}
    def rows(sel):
        for i in sel:
            b = C[i + 1:]
            if not len(b): continue
            for l in LASTS:
                x = (SIGMA - C[i] - b - l - 3 * A_U) / D; m = (x >= xlo) & (x <= xhi)
                if not m.any(): continue
                sc = ((G[i] * GL[l]) * G[i + 1:][m]).sum(axis=1); y = tau_est(sc) * ih3(x[m]) * typ(x[m])
                yield i, l, np.nonzero(m)[0] + i + 1, y
    rng = random.Random(1); samp = rng.sample(range(N), min(N, 600)); ys = np.concatenate([y for _, _, _, y in rows(samp)] or [np.zeros(1)])
    tot = len(ys) * N / len(samp); thr = 0.0 if tot <= nkeep else float(np.quantile(ys, 1.0 - nkeep / tot))
    out_y = []; out_i = []; out_j = []; out_l = []
    for i, l, j, y in rows(range(N)):
        k = y >= thr
        if k.any(): out_y.append(y[k]); out_j.append(j[k]); out_i.append(np.full(int(k.sum()), i)); out_l.append(np.full(int(k.sum()), l))
    y = np.concatenate(out_y); ii = np.concatenate(out_i); jj = np.concatenate(out_j); ll = np.concatenate(out_l); o = np.argsort(-y)
    log('  固定に使える塊 %d 個（小数部 %.2f 以上・打てる文字）、組 約 %.3g 通りのうち上位 %d 通り（見積もりの比 %.2f〜1）' % (N, fmin, tot, len(o), float(y[o[-1]] / y[o[0]])))
    done = set(done); res = []; ys2 = []
    for k in o:
        a, b = int(C[ii[k]]), int(C[jj[k]]); kb = b''.join(b4(c) for c in sorted((a, b))) + bytes([int(ll[k])])
        if kb.hex() not in done: res.append(kb); ys2.append(float(y[k]))
    return res, ys2
if __name__ == '__main__':
    from vshort_plan import tau
    hi = load_hi(); rng = random.Random(3)
    ch = chunks_for_pool(hi, 0.90); GL = {l: 1.0 + np.array(chunk_vals(l)[0]) for l in LASTS}
    print('点数（(1+f) の積の和）と τ の実測'); pts = []
    for fa, fb in ((0.90, 0.92), (0.92, 0.94), (0.93, 0.95), (0.95, 0.97), (0.97, 1.0)):
        sel = [x for x in ch if fa <= min(x[1]) < fb]
        for _ in range(4):
            a, b = rng.sample(sel, 2); l = rng.choice(LASTS)
            sc = float(((1 + np.array(a[1])) * (1 + np.array(b[1])) * GL[l]).sum()); t = tau([a[0], b[0], l], nsamp=400000)[0]; pts.append((sc, t))
            print('  小数部の下限 %.2f〜%.2f  点数 %.2f  τ 実測 %.4f  見積もり %.4f' % (fa, fb, sc, t, float(tau_est(sc))))
    pool, ys = ranked_pool(hi)
    print('組 %d 通り' % len(pool))
