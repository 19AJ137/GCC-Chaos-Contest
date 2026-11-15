"""head_enum.py — 先頭の数語のフラグの全通りを試し、MX の候補ごとに 2 本の文章の点を比べる（先頭は最上位ビットの手がかりがないため）"""
import sys, os, itertools, pickle, math
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
import lm2
from lm2 import A7, A7IDX
from fulldec2 import xor4, Z4
A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
Xw = [xor4(A[4 * n:4 * n + 4], E[4 * n:4 * n + 4]) for n in range(256)]
print('X_0〜X_7: ' + ' '.join(x.hex() for x in Xw[:8]))
le = lm2.load('lm2_e.pkl'); la = lm2.load('lm2_a.pkl')
NR = 7
def mxseq(flags):
    jout = -1; rin = Z4; out = []
    for n, (f, g) in enumerate(flags):
        x = Xw[n]; px = xor4(x, Xw[jout] if jout >= 0 else Z4) if f else x; mx = xor4(px, rin) if g else px
        if f: jout = n
        if g: rin = px
        out.append(mx)
    return b''.join(out)
cands = {}
for fl in itertools.product(((0, 0), (0, 1), (1, 0), (1, 1)), repeat=NR):
    if fl[0] != (0, 0) or fl[2] != (0, 0): continue       # 語 0・2 は A = B
    cands.setdefault(mxseq(fl), []).append(fl)
print('MX の候補（語 0〜%d）: %d 通り' % (NR - 1, len(cands)))
cribE = {21: 0x0d, 22: 0x0a, 23: 0x20, 24: 0x20, 25: ord('T'), 26: ord('h'), 27: ord('a')}
def best_split(z, W=300):
    K = 4; beam = {(b'', b''): (0.0, b'', b'')}
    for i, zi in enumerate(z):
        nb = {}
        for (ce, ca), (sc, se, sa) in beam.items():
            ve = le.logp_vec7(ce); va = la.logp_vec7(ca)
            for ei, e in enumerate(A7):
                if i in cribE and e != cribE[i]: continue
                a = e ^ zi; ai = A7IDX.get(a)
                if ai is None: continue
                s = sc + ve[ei] + va[ai]; k = ((ce + bytes([e]))[-K:], (ca + bytes([a]))[-K:])
                if k not in nb or s > nb[k][0]: nb[k] = (s, se + bytes([e]), sa + bytes([a]))
        if not nb: return None
        beam = dict(sorted(nb.items(), key=lambda kv: -kv[1][0])[:W])
    return max(beam.values())
res = []
for z, fls in cands.items():
    r = best_split(z)
    if r: res.append((r[0], z, r[1], r[2], fls[0], len(fls)))
res.sort(reverse=True)
for sc, z, e, a, fl, nf in res[:14]:
    print('%8.1f  E=%r  A=%r  フラグの例 %s（%d 通り）' % (sc, e, a, ''.join('%d%d ' % f for f in fl), nf)); print('          MX ' + ' '.join(z[k:k + 4].hex() for k in range(0, len(z), 4)))
