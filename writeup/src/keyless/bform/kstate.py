"""kstate.py — 外側の鍵 k（写像 2）と内側の鍵 v（写像 1）の状態から、全 256 語の鍵ストリームと連鎖フラグを作る（2026-11-12）
  k: kexact2 の出力（lo, Q, A, C, u_33, v_33）から前向きに語 255 まで、後ろ向きに語 0 まで（逆写像は 3 次方程式の根のまわりで厳密な候補を探す）。
  v: v_family.pkl。
  python3 -I kstate.py [kexact_out.txt]"""
import sys, os, pickle, struct, math, re, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
from solve_v import d2b, u2d, d2u, setlow, kword
import numpy as np
def step2(Q, A, C, lo, u, v):
    a = Q * u; b = a * u; nu = A * u - a * v; nv = C * v + b
    return setlow(nu, lo), setlow(nv, lo)
def back2(Q, A, C, lo, u1, v1, kw_prev=None, R=48):
    """(u1, v1) の 1 つ前の状態の候補（厳密に一致するもの）"""
    out = []
    roots = np.roots([Q * Q, 0.0, A * C - Q * v1, -C * u1])
    for x in roots:
        if abs(x.imag) > 1e-9 * max(1.0, abs(x.real)): continue
        x = float(x.real); y = (v1 - Q * x * x) / C
        if not (abs(x) < 4 and -1.0 < y < 7.0): continue
        ux = (d2u(x) & ~0xFFFF) | lo; uy = (d2u(y) & ~0xFFFF) | lo
        for du in range(-R, R + 1):
            u0 = u2d(ux + du * 0x10000)
            for dv in range(-R, R + 1):
                v0 = u2d(uy + dv * 0x10000)
                if kw_prev is not None and kword(v0)[0] != kw_prev: continue
                if step2(Q, A, C, lo, u0, v0) == (u1, v1): out.append((u0, v0))
    return out
def parse(fn):
    sols = []
    for ln in open(fn):
        m = re.match(r'lo=([0-9a-f]+) Q=\S+ \(([0-9a-f]+)\) A=\S+ \(([0-9a-f]+)\) C=\S+ \(([0-9a-f]+)\) u_(\d+)=\S+ \(([0-9a-f]+)\) v_\d+=\S+ \(([0-9a-f]+)\)\s+一致 (\d+) 不一致 (\d+)', ln)
        if m: sols.append(dict(lo=int(m.group(1), 16), Q=u2d(int(m.group(2), 16)), A=u2d(int(m.group(3), 16)), C=u2d(int(m.group(4), 16)), n0=int(m.group(5)), u=u2d(int(m.group(6), 16)), v=u2d(int(m.group(7), 16)), ok=int(m.group(8)), bad=int(m.group(9))))
    return sols
def full_k(s, KK, flk=None, verbose=False):
    """解 s から全 256 語。返り値: [(ks, flags, 語 0 の状態)]。後ろ向きは、分かっている語・フラグに合う枝だけ"""
    Q, A, C, lo, n0 = s['Q'], s['A'], s['C'], s['lo'], s['n0']
    u, v = s['u'], s['v']; fw = {n0: (kword(v)[0], 0 if kword(v)[1] > 0xE6 else 1)}
    for n in range(n0 + 1, 256):
        u, v = step2(Q, A, C, lo, u, v); kw, b6 = kword(v); fw[n] = (kw, 0 if b6 > 0xE6 else 1)
    front = [((s['u'], s['v']), {})]
    for n in range(n0 - 1, -1, -1):
        nf = []
        for (u1, v1), sq in front:
            for u0, v0 in back2(Q, A, C, lo, u1, v1, KK.get(n)):
                kw, b6 = kword(v0); f = 0 if b6 > 0xE6 else 1
                s2 = dict(sq); s2[n] = (kw, f); nf.append(((u0, v0), s2))
        if verbose: print('   語 %d: 枝 %d' % (n, len(nf)))
        front = nf[:256]
        if not front: return []
    out = []
    for st0, sq in front:
        full = dict(fw); full.update(sq); out.append(([full[n][0] for n in range(256)], [full[n][1] for n in range(256)], st0))
    return out
def main():
    fn = sys.argv[1] if len(sys.argv) > 1 else os.path.join(D, 'kexact_out.txt')
    sols = parse(fn); print('k の厳密な解: %d 個。不一致の数の内訳 %s、lo の範囲 %04x〜%04x' % (len(sols), dict(collections.Counter(s['bad'] for s in sols)), min(s['lo'] for s in sols), max(s['lo'] for s in sols)))
    det = pickle.load(open(os.path.join(D, 'ks_Ek_Ev.pkl'), 'rb'))['det']; KK = {s[1]: v for s, v in det.items() if s[0] == 'K'}
    from headdata import load_head
    fl = load_head(os.path.join(D, 'head132_fixed.json'))['flags']; flk = {n: fl[n][0] for n in fl if n < 132}
    best = min(s['bad'] for s in sols); fam = collections.defaultdict(list)
    for s in [s for s in sols if s['bad'] == best][::max(1, len(sols) // 40)]:
        # 合わない語を特定
        u, v = s['u'], s['v']; mism = []
        for n in range(s['n0'] + 1, 130):
            u, v = step2(s['Q'], s['A'], s['C'], s['lo'], u, v)
            if n in KK and kword(v)[0] != KK[n]: mism.append((n, '%08x' % (kword(v)[0] ^ KK[n])))
        s['mism'] = mism
        KK2 = dict(KK)
        for n, _ in mism: KK2.pop(n, None)
        for ks, flags, st0 in full_k(s, KK2, flk):
            fam[(tuple(ks), tuple(flags))].append((s['lo'], st0))
    print('不一致の語（語、XOR）: %s' % sols[0].get('mism'))
    print('全 256 語の並びの種類: %d' % len(fam))
    out = []
    for (ks, flags), mem in sorted(fam.items(), key=lambda kv: -len(kv[1])):
        print('  %d 個の解が同じ並び。連鎖なしの語の数 %d。フラグ（語 0〜131）が復号で決めたものと合わない語: %s' % (len(mem), flags.count(0), [n for n in range(132) if flk[n] != flags[n]]))
        out.append(dict(ks=list(ks), flags=list(flags), members=mem))
    pickle.dump(dict(fam=out, sols=sols), open(os.path.join(D, 'k_family.pkl'), 'wb'))
if __name__ == '__main__': main()
