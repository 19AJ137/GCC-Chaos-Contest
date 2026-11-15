"""edec2.py — 厳密に決めた 2 ブロック目の状態で、E の残り 314 バイトを復号する。"""
import sys, os, pickle, struct, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
import eblock2 as E2M, eblock2x as X2
V = E2M.V; u2d, d2u, setlow = V.u2d, V.d2u, V.setlow
r = pickle.load(open(os.path.join(D, 'eblock2_exact.pkl'), 'rb')); sols, groups = r['sols'], r['groups']
gs = collections.Counter(s['group'] for s in sols); print('解の v の組: %s' % dict(gs))
for gi, (ks, fl, mem) in enumerate(groups):
    print('v の組 %d: 候補 %d 個、lo_v %04x〜%04x、dp %s、idx %s' % (gi, len(mem), min(m[0] for m in mem), max(m[0] for m in mem), sorted(set(m[2] for m in mem)), sorted(set(m[1] for m in mem))))
outs = collections.Counter()
for s in sols:
    ksv, flv, mem = groups[s['group']]
    # k 層の 79 語（78 語と端数用の 1 語）
    u, v = s['u0'], s['v0']; kk = []; fk = []
    for n in range(E2M.NW2 + 1):
        if n: u, v = X2.step2(s['Q'], s['A'], s['C'], s['lo'], u, v)
        w, f = X2.kword(v); kk.append(w); fk.append(f)
    # 復号: 外側 k（レジスタ RK0 から）→ 内側 v（レジスタ RV0 から）
    reg = E2M.RK0; cp = []
    for n in range(E2M.NW2):
        x = E2M.Ew2[n] ^ kk[n]
        if fk[n]: x ^= reg; reg = E2M.Ew2[n]
        cp.append(x)
    reg = E2M.RV0; m = []
    for n in range(E2M.NW2):
        x = cp[n] ^ ksv[n]
        if flv[n]: x ^= reg; reg = cp[n]
        m.append(x)
    body = V.tobytes(m)
    # 端数 2 バイト: c = ks[j] ^ 入力 ^ bl、bl の初期値は IV の 2 バイト目、復号では bl = 入力（暗号文）
    kb = kk[E2M.NW2].to_bytes(4, 'big'); vb = ksv[E2M.NW2].to_bytes(4, 'big'); blk = (V.K_IV >> 16) & 0xFF; blv = 0x34; tail = bytearray()
    for j, e in enumerate(E2M.TAIL):
        c1 = e ^ kb[j] ^ blk; blk = e; mm = c1 ^ vb[j] ^ blv; blv = c1; tail.append(mm)
    outs[body + bytes(tail)] += 1
print('復号結果の種類: %d' % len(outs))
for pt, n in outs.most_common():
    print('--- %d 個の解が同じ結果（%d バイト）---' % (n, len(pt))); print(repr(pt)); print(pt.decode('cp932', 'replace'))
best = outs.most_common(1)[0][0]; open(os.path.join(D, 'plain_E_block2.bin'), 'wb').write(best)
