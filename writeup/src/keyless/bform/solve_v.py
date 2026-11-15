"""solve_v.py — 内側の鍵 v（写像 1）の状態を、解いた鍵ストリームの語から厳密に復元する（2026-11-12）
  写像 1（gcc.py の step）: y ← set_low16(y*y - p)、出力 x = (1.0 + y)*0.5。鍵ストリームの語は x のバイト (4,2,3,5)、連鎖なし ⇔ x のバイト 6 > 0xE6。
  y・p の下位 16 ビットは鍵で決まる定数 b0・b1。
  手順: 連続する 3 語で b0b1（65,536 通り）と y・p の候補を厳密に当て、先の語で確かめる。見つかったら前後に写像を回して全 256 語を出す。
  python3 -I solve_v.py ks_Ek_Ev.pkl head132_fixed.json [開始語]"""
import sys, os, pickle, struct, math
D = os.path.dirname(os.path.abspath(__file__))
def d2b(x): return struct.pack('<d', x)
def b2d(b): return struct.unpack('<d', bytes(b))[0]
def u2d(u): return struct.unpack('<d', struct.pack('<Q', u))[0]
def d2u(x): return struct.unpack('<Q', struct.pack('<d', x))[0]
def setlow(x, lo): return u2d((d2u(x) & ~0xFFFF) | lo)
def kword(x):
    b = d2b(x); return (b[4] << 24) | (b[2] << 16) | (b[3] << 8) | b[5], b[6]
def step(y, p, lo):
    y = setlow(y * y - p, lo); return y, (1.0 + y) * 0.5
def y_cands(kw, b6s, lo):
    """x のバイト 2〜5 が kw、バイト 6 が b6s のどれか、符号は両方、になる y（下位 16 ビット = lo）を全部"""
    out = []
    b4 = (kw >> 24) & 0xFF; b2 = (kw >> 16) & 0xFF; b3 = (kw >> 8) & 0xFF; b5 = kw & 0xFF
    for b7 in (0x3F, 0xBF):
        for b6 in b6s:
            base = (b2 << 16) | (b3 << 24) | (b4 << 32) | (b5 << 40) | (b6 << 48) | (b7 << 56)
            xlo = u2d(base); xhi = u2d(base | 0xFFFF)
            ya = 2 * xlo - 1; yb = 2 * xhi - 1
            if ya > yb: ya, yb = yb, ya
            if yb < -2.01 or ya > 2.01: continue
            # ya〜yb の間で下位 16 ビットが lo の倍精度を列挙（符号をまたぐ区間は分ける）
            segs = []
            if ya < 0: segs.append((ya, min(yb, -0.0)))
            if yb > 0: segs.append((max(ya, 0.0), yb))
            for a, b in segs:
                ua, ub = sorted((d2u(a), d2u(b)))
                u = (ua & ~0xFFFF) | lo
                if u < ua: u += 0x10000
                cnt = 0
                while u <= ub and cnt < 4096:
                    y = u2d(u); x = (1.0 + y) * 0.5
                    k2, c6 = kword(x)
                    if k2 == kw and c6 == b6 and d2b(x)[7] == b7: out.append(y)
                    u += 0x10000; cnt += 1
    return out
def main():
    sys.path.insert(0, D); from headdata import load_head
    det = pickle.load(open(sys.argv[1], 'rb'))['det']; fl = load_head(sys.argv[2])['flags']
    KV = {s[1]: v for s, v in det.items() if s[0] == 'V'}
    n0 = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    assert all(n0 + i in KV for i in range(8)), '連続する 8 語が要る'
    def b6set(n): return range(0x00, 0xE7) if fl[n][1] else range(0xE7, 0x100)
    # まず近似で b6・符号を絞る（下位を 0x8000 とした近似値で p の一致を見る）
    import itertools
    approx = []
    for n in (n0, n0 + 1, n0 + 2):
        c = []
        for b7 in (0x3F, 0xBF):
            for b6 in b6set(n):
                kw = KV[n]; base = 0x8000 | (((kw >> 16) & 0xFF) << 16) | (((kw >> 8) & 0xFF) << 24) | (((kw >> 24) & 0xFF) << 32) | ((kw & 0xFF) << 40) | (b6 << 48) | (b7 << 56)
                x = u2d(base)
                if -0.5001 < x < 1.5001: c.append((2 * x - 1, b6, b7))
        approx.append(c)
    best = None
    for ya, b6a, b7a in approx[0]:
        for yb, b6b, b7b in approx[1]:
            p = ya * ya - yb
            if not (1.9 < p < 2.0001): continue
            for yc, b6c, b7c in approx[2]:
                dd = abs(yb * yb - yc - p)
                if best is None or dd < best[0]: best = (dd, p, (b6a, b7a), (b6b, b7b), (b6c, b7c), ya, yb, yc)
    print('近似: p ≈ %.12f（ずれ %.2e）、語 %d〜%d の x のバイト 6・7: %s' % (best[1], best[0], n0, n0 + 2, [('%02x' % a, '%02x' % b) for a, b in best[2:5]]))
    (b6a, b7a), (b6b, b7b) = best[2], best[3]; papprox = best[1]
    sols = []
    for lo in range(65536):
        Ya = [y for y in y_cands(KV[n0], [b6a], lo) if (d2b((1.0 + y) * 0.5)[7] == b7a)]
        if not Ya: continue
        Yb = [y for y in y_cands(KV[n0 + 1], [b6b], lo) if (d2b((1.0 + y) * 0.5)[7] == b7b)]
        for ya in Ya:
            q = ya * ya
            for yb in Yb:
                p0 = q - yb; u0 = (d2u(p0) & ~0xFFFF) | lo
                for du in (-0x20000, -0x10000, 0, 0x10000, 0x20000):
                    p = u2d(u0 + du)
                    if abs(p - papprox) > 1e-6: continue
                    if setlow(q - p, lo) != yb: continue
                    # 先の語で確かめる
                    y = yb; ok = 0
                    for i in range(2, 8):
                        y, x = step(y, p, lo)
                        if kword(x)[0] == KV[n0 + i]: ok += 1
                        else: break
                    if ok >= 5:
                        # 分かっている語を最後まで照合する
                        y = yb; ok = 0; bad = 0; n = n0 + 1
                        while n + 1 <= max(KV):
                            y, x = step(y, p, lo); n += 1
                            if n in KV:
                                if kword(x)[0] == KV[n]: ok += 1
                                else: bad += 1
                        sols.append((lo, ya, p, ok, bad))
    print('解: %d 個' % len(sols))
    sols.sort(key=lambda t: -t[3])
    for lo, ya, p, ok, bad in sols[:8]: print('   b0b1 = %04x（b0 = %02x、b1 = %02x）、y_%d = %r（%s）、p = %r（%s）、先の語: 一致 %d、不一致 %d' % (lo, lo & 0xFF, lo >> 8, n0, ya, d2b(ya).hex(), p, d2b(p).hex(), ok, bad))
    if sols:
        lo, ya, p, ok, bad = sols[0]; y = ya; n = n0; mism = []
        while n + 1 <= max(KV):
            y, x = step(y, p, lo); n += 1
            if n in KV and kword(x)[0] != KV[n]: mism.append((n, '%08x' % kword(x)[0], '%08x' % KV[n], '%08x' % (kword(x)[0] ^ KV[n])))
        print('   最良の解で合わない語（語、写像の出力、解いた値、XOR）: %s' % mism)
    pickle.dump(dict(sols=sols, n0=n0), open(os.path.join(D, 'v_state.pkl'), 'wb'))
if __name__ == '__main__': main()
