"""kinit.py — 鍵 k（写像 2）の初期値（鍵スケジュールの k780・k788 = 空回しの前の u・v）を、E の 2 ブロック目の最初の状態から逆算し、
  そこから 1 ブロック目の定数で 70 歩（空回し 69 = S_k mod 79 と、語 0 の 1 歩）進めて、語 0 の状態に厳密に着くかを確かめる。2026-11-12
  python3 -I -B kinit.py"""
import sys, os, pickle, struct, math, collections
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
import eblock2 as E2M, eblock2x as X2
V = E2M.V; u2d, d2u, setlow = V.u2d, V.d2u, V.setlow
r = pickle.load(open(os.path.join(D, 'eblock2_exact.pkl'), 'rb')); sols = [s for s in r['sols'] if s['u0'] > 0]
Sk = int.from_bytes(V.K_IV.to_bytes(4, 'big'), 'little'); warm = max(11, Sk % 79); print('S_k = %08x、空回し %d 回' % (Sk, warm))
K, F = V.stream_k()
res = []
for s in sols:
    lo = s['lo']
    if not (0x7ec1 <= lo <= 0x7eef): continue
    Q2, A2, C2, u1, v1 = s['Q'], s['A'], s['C'], s['u0'], s['v0']
    # 2 ブロック目の最初の 1 歩を逆に: Q'^2 u^3 + (A'C' − Q' v1) u − C' u1 = 0、v = (v1 − Q' u^2)/C'
    roots = np.roots([Q2 * Q2, 0.0, A2 * C2 - Q2 * v1, -C2 * u1]); inits = []
    for x in roots:
        if abs(x.imag) > 1e-12: continue
        ui = float(x.real); vi = (v1 - Q2 * ui * ui) / C2
        if 0.0 <= ui < 1.0 and 0.0 <= vi < 1.0: inits.append((ui, vi))
    # 1 ブロック目の定数（lo ごとに下位 16 ビットを置き換える）
    Q = u2d((V.K_Q & ~0xFFFF) | lo); A = u2d((V.K_A & ~0xFFFF) | lo); C = u2d((V.K_C & ~0xFFFF) | lo)
    tu = d2u(abs(u2d(V.K_U0))) >> 16; tv = d2u(u2d(V.K_V0)) >> 16
    for ui, vi in inits:
        a = Q * ui; nu = A * ui - a * vi; nv = C * vi + a * ui          # 1 ブロック目の最初の 1 歩（近似）
        iu = d2u(abs(nu)) >> 16; iv = d2u(nv) >> 16; hit = []
        for du in range(-6, 7):
            for dv in range(-6, 7):
                u = math.copysign(u2d(((iu + du) << 16) | lo), nu); v = u2d(((iv + dv) << 16) | lo)
                for _ in range(warm): u, v = X2.step2(Q, A, C, lo, u, v)       # 最初の 1 歩のあと、さらに 69 歩で語 0 の状態
                if (d2u(abs(u)) >> 16) == tu and (d2u(v) >> 16) == tv: hit.append((du, dv, u, v))
        res.append((lo, ui, vi, hit))
print('lo ごとの結果（1 ブロック目・2 ブロック目の両方に合う lo）: %d 個' % len(res))
ok = [x for x in res if x[3]]
print('  2 ブロック目から逆算した初期値を 70 歩進めて、語 0 の状態（上位 48 ビット）に着いた lo: %d 個（%s）' % (len(ok), ('%04x〜%04x' % (min(x[0] for x in ok), max(x[0] for x in ok))) if ok else '-'))
for lo, ui, vi, hit in res[:3] + res[-3:]:
    print('   lo %04x: 初期値 u = k780 ≈ %.15f、v = k788 ≈ %.15f、着いた候補 %s' % (lo, ui, vi, [(h[0], h[1], '%+.6f' % h[2], '%.6f' % h[3]) for h in hit]))
print('  語 0 の状態（verify.py）: u = %.6f v = %.6f' % (u2d(V.K_U0), u2d(V.K_V0)))
pickle.dump(res, open(os.path.join(D, 'kinit.pkl'), 'wb'))
