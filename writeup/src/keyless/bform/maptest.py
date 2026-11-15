"""maptest.py — 解いた鍵ストリームの語が GCC6 の写像の出力として成り立つかを確かめる（2026-11-12）
  鍵ストリームの語 K_n = 出力の倍精度 x_n のバイト (4, 2, 3, 5)。バイト 0・1 は鍵で決まる定数（b0・b1）、バイト 6・7 は未知（連鎖なし ⇔ バイト 6 > 0xE6）。
  写像 0: x' = r·x·(1-x)（r は 3.9〜4.0）。連続する 3 語で、バイト 6 を全通り試し、r = x_{n+1}/(x_n(1-x_n)) = x_{n+2}/(x_{n+1}(1-x_{n+1})) が合うか。
  写像 1: y' = y² - p（p は 1.9〜2.0）、x = (1+y)/2。p = y_n² - y_{n+1} = y_{n+1}² - y_{n+2} が合うか。
  python3 -I maptest.py ks_Ek_Ev.pkl,ks_Dk_Dv.pkl head132_fixed.json"""
import sys, os, pickle, struct
import numpy as np
D = os.path.dirname(os.path.abspath(__file__))
def dbl(kw, b6, b7=0x3F, low=0x8000):
    """kw: 鍵ストリームの語（上位バイトから b4 b2 b3 b5）。b6 は配列"""
    b4 = (kw >> 24) & 0xFF; b2 = (kw >> 16) & 0xFF; b3 = (kw >> 8) & 0xFF; b5 = kw & 0xFF
    base = low | (b2 << 16) | (b3 << 24) | (b4 << 32) | (b5 << 40) | (b7 << 56)
    u = (np.asarray(b6, dtype=np.uint64) << np.uint64(48)) | np.uint64(base)
    return u.view(np.float64)
def cand_b6(flag_chained):
    return np.arange(0x00, 0xE7, dtype=np.uint64) if flag_chained else np.arange(0xE7, 0xF0, dtype=np.uint64)
def test_map0(kws, fls):
    """kws: 連続する 3 語、fls: その連鎖フラグ（1 = 連鎖あり）"""
    ca, cb, cc = (cand_b6(f) for f in fls)
    xa = dbl(kws[0], ca); xb = dbl(kws[1], cb); xc = dbl(kws[2], cc)
    r1 = xb[None, :] / (xa * (1 - xa))[:, None]          # (na, nb)
    r2 = xc[None, :] / (xb * (1 - xb))[:, None]          # (nb, nc)
    hits = []
    ok1 = (r1 > 3.9) & (r1 < 4.0001)
    for ia, ib in zip(*np.nonzero(ok1)):
        d = np.abs(r2[ib] - r1[ia, ib]); ic = int(np.argmin(d))
        if d[ic] < 2e-8: hits.append((float(r1[ia, ib]), float(d[ic]), int(ca[ia]), int(cb[ib]), int(cc[ic])))
    return hits
def test_map1(kws, fls):
    hits = []
    # x は [-0.5, 1.5]。バイト 7 は 0x3F（正）か 0xBF（負）。x ≥ 1 のときバイト 6 は 0xF0〜0xF7（連鎖なし）
    def cands(f):
        out = []
        for b7 in (0x3F, 0xBF):
            rng = (np.arange(0x00, 0xE7, dtype=np.uint64) if f else np.arange(0xE7, 0xF8, dtype=np.uint64))
            out.append((b7, rng))
        return out
    X3 = []
    for kw, f in zip(kws, fls):
        xs = np.concatenate([dbl(kw, rng, b7) for b7, rng in cands(f)]); xs = xs[(xs > -0.5001) & (xs < 1.5001)]; X3.append(2 * xs - 1)
    ya, yb, yc = X3
    p1 = (ya * ya)[:, None] - yb[None, :]; p2 = (yb * yb)[:, None] - yc[None, :]
    ok1 = (p1 > 1.9) & (p1 < 2.0001)
    for ia, ib in zip(*np.nonzero(ok1)):
        d = np.abs(p2[ib] - p1[ia, ib]); ic = int(np.argmin(d))
        if d[ic] < 2e-8: hits.append((float(p1[ia, ib]), float(d[ic]), float(ya[ia]), float(yb[ib]), float(yc[ic])))
    return hits
def main():
    sys.path.insert(0, D); from headdata import load_head
    fl = load_head(sys.argv[2])['flags']
    for fn in sys.argv[1].split(','):
        d = pickle.load(open(fn, 'rb')); det = d['det']
        print('=== %s（B の形 %s）' % (fn, '・'.join(d['form'])))
        for key, which, lab in (('K', 0, '外側の鍵 k'), ('V', 1, '内側の鍵 v')):
            ns = sorted(s[1] for s in det if s[0] == key); trip = [n for n in ns if n + 1 in ns and n + 2 in ns]
            h0 = h1 = 0; rs0 = []; rs1 = []
            for n in trip:
                kws = [det[(key, n + i)] for i in range(3)]; fls = [fl[n + i][which] for i in range(3)]
                a = test_map0(kws, fls); b = test_map1(kws, fls)
                if a: h0 += 1; rs0.append(a[0][0])
                if b: h1 += 1; rs1.append(b[0][0])
            print('  %s: 連続する 3 語の組 %d 個。写像 0 として合う組 %d（r の例 %s）、写像 1 として合う組 %d（p の例 %s）' % (lab, len(trip), h0, ['%.12f' % r for r in rs0[:4]], h1, ['%.12f' % r for r in rs1[:4]]))
if __name__ == '__main__': main()
