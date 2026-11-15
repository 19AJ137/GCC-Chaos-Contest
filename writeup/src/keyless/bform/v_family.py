"""v_family.py — 内側の鍵 v（写像 1）: 解いた鍵ストリーム全部に合う (b0b1, y, p) の集まりを調べ、全 256 語の鍵ストリームと連鎖フラグを出す（2026-11-12）
  python3 -I v_family.py"""
import sys, os, pickle, struct, math, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
from solve_v import d2b, u2d, d2u, setlow, kword, step
def back_cands(y, p, lo):
    """y = setlow(yp*yp - p) となる yp（下位 16 ビット = lo）の候補"""
    out = []
    u = d2u(y); ylo = u2d(u & ~0xFFFF); yhi = u2d(u | 0xFFFF)
    a, b = sorted((ylo + p, yhi + p))
    if b < 0: return out
    ra = math.sqrt(max(a, 0.0)); rb = math.sqrt(b)
    for sgn in (1, -1):
        lo_, hi_ = (ra, rb) if sgn > 0 else (-rb, -ra)
        ua, ub = sorted((d2u(lo_), d2u(hi_)))
        uu = ((ua - 0x40000) & ~0xFFFF) | lo
        while uu <= ub + 0x40000:
            yp = u2d(uu)
            if setlow(yp * yp - p, lo) == y: out.append(yp)
            uu += 0x10000
    return out
def main():
    S = pickle.load(open(os.path.join(D, 'v_state.pkl'), 'rb')); n0 = S['n0']; sols = [s for s in S['sols'] if s[4] == 0]
    det = pickle.load(open(os.path.join(D, 'ks_Ek_Ev.pkl'), 'rb'))['det']; KV = {s[1]: v for s, v in det.items() if s[0] == 'V'}
    los = sorted(set(s[0] for s in sols)); print('全部の語に合う解: %d 個。b0b1 の範囲 %04x〜%04x（%d 通り）、p の種類 %d' % (len(sols), los[0], los[-1], len(los), len(set(s[2] for s in sols))))
    fam = collections.defaultdict(list)
    for lo, ya, p, ok, bad in sols:
        # 前向き: 語 n0 から 255 まで
        y = ya; seq = {}; x = (1.0 + y) * 0.5; seq[n0] = (kword(x)[0], 0 if kword(x)[1] > 0xE6 else 1)
        for n in range(n0 + 1, 256):
            y, x = step(y, p, lo); kw, b6 = kword(x); seq[n] = (kw, 0 if b6 > 0xE6 else 1)
        # 後ろ向き: 語 n0-1 から 0 まで。分かっている語に合う候補だけ残す
        front = [(ya, {})]
        for n in range(n0 - 1, -1, -1):
            nf = []
            for y, sq in front:
                for yp in back_cands(y, p, lo):
                    x = (1.0 + yp) * 0.5; kw, b6 = kword(x)
                    if n in KV and kw != KV[n]: continue
                    s2 = dict(sq); s2[n] = (kw, 0 if b6 > 0xE6 else 1); nf.append((yp, s2))
            front = nf[:64]
            if not front: break
        for y0, sq in front:
            full = dict(seq); full.update(sq)
            key = tuple(full[n] for n in range(256))
            fam[key].append((lo, p, y0))
    print('全 256 語の並びの種類: %d' % len(fam))
    out = []
    for key, members in sorted(fam.items(), key=lambda kv: -len(kv[1])):
        los = sorted(set(m[0] for m in members))
        print('  %d 個の (b0b1, y) が同じ並び。b0b1 = %04x〜%04x。連鎖なしの語の数 %d。語 0 の前の y の例 %r、p = %r' % (len(members), los[0], los[-1], sum(1 for kw, f in key if f == 0), members[0][2], members[0][1]))
        out.append(dict(ks=[k for k, f in key], flags=[f for k, f in key], members=members))
    pickle.dump(out, open(os.path.join(D, 'v_family.pkl'), 'wb'))
    # 違いの出る語
    if len(out) > 1:
        for i in range(1, min(len(out), 6)):
            diff = [n for n in range(256) if out[0]['ks'][n] != out[i]['ks'][n]]
            print('  並び 0 と並び %d の違う語: %d 個（最初 %s）' % (i, len(diff), diff[:6]))
if __name__ == '__main__': main()
