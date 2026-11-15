"""back1.py — 写像 1（y' = setlow(y*y − p)）を逆向きにたどる（numpy 版）。格子（下位 16 ビット = lo の倍精度の数）の上で、
  厳密に一致する原像を全部、区間の節として持つ。節 = (符号, idx の範囲 [i1, i2])。idx = |y| のビット表現の上位 48 ビット。
  y*y − p は |y| について単調なので、原像は符号ごとに 1 つの連続した区間になる。両端だけ厳密に決める。"""
import struct, math
import numpy as np
def u2d(u): return struct.unpack('<d', struct.pack('<Q', u))[0]
def d2u(x): return struct.unpack('<Q', struct.pack('<d', x))[0]
def setlow(x, lo): return u2d((d2u(x) & ~0xFFFF) | lo)
def tkey(w):
    u = d2u(w); return (-1 if u >> 63 else 1), (u & 0x7FFFFFFFFFFFFFFF) >> 16
def G(y, p, lo): return setlow(y * y - p, lo)
M63 = np.uint64(0x7FFFFFFFFFFFFFFF)
def _f(idx, p, lo):
    """格子の idx（int64 配列）について、w = y*y − p の単調な整数キー（w が大きいほど大きい）"""
    y = ((idx.astype(np.uint64) << np.uint64(16)) | np.uint64(lo)).view(np.float64)
    w = y * y - p; u = w.view(np.uint64); k = ((u & M63) >> np.uint64(16)).astype(np.int64)
    return np.where((u >> np.uint64(63)) != 0, -k - 1, k)
def back_nodes(tsg, t1, t2, p, lo):
    """節の配列（tsg: ±1、t1・t2: int64）の原像。返り値も (sg, i1, i2) の配列。符号 ± の両方。"""
    mlo = (t1.astype(np.uint64) << np.uint64(16)).view(np.float64); mhi = ((t2.astype(np.uint64) << np.uint64(16)) | np.uint64(0xFFFF)).view(np.float64)
    pos = tsg > 0; wlo = np.where(pos, mlo, -mhi); whi = np.where(pos, mhi, -mlo)
    ok = (whi + p) >= 0
    if not ok.any(): z = np.zeros(0, dtype=np.int64); return z, z, z
    wlo = wlo[ok]; whi = whi[ok]; tl = np.where(pos, t1, -t2 - 1)[ok]; th = np.where(pos, t2, -t1 - 1)[ok]
    a = np.sqrt(np.maximum(wlo + p, 0.0)); b = np.sqrt(whi + p)
    ia = (a.view(np.uint64) >> np.uint64(16)).astype(np.int64); ib = (b.view(np.uint64) >> np.uint64(16)).astype(np.int64)
    # 下端: f(i) >= tl となる最小の i（見積もり ia の ±4 の中）
    L = ia + 5
    for d in range(4, -5, -1):
        i = np.maximum(ia + d, 0); L = np.where(_f(i, p, lo) >= tl, i, L)
    U = ib - 5
    for d in range(-4, 5):
        i = np.maximum(ib + d, 0); U = np.where(_f(i, p, lo) <= th, i, U)
    # 見積もりが外れていないかの検算（L−1 は条件を満たさない、U+1 も満たさない）。外れた節は二分法でやり直す
    bad = ((L > 0) & (_f(np.maximum(L - 1, 0), p, lo) >= tl) & (L == np.maximum(ia - 4, 0)) & (ia - 4 > 0)) | ((_f(U + 1, p, lo) <= th) & (U == ib + 4)) | (L == ia + 5) & (_f(np.maximum(ia + 5, 0), p, lo) <= th) & False
    if bad.any():
        for j in np.nonzero(bad)[0]:
            lo_i, hi_i = 0, int(ib[j]) + 64
            Lb, Rb = lo_i, hi_i
            while Lb < Rb:
                m = (Lb + Rb) // 2
                if int(_f(np.array([m]), p, lo)[0]) < int(tl[j]): Lb = m + 1
                else: Rb = m
            L[j] = Lb
            Lb, Rb = lo_i, hi_i
            while Lb < Rb:
                m = (Lb + Rb + 1) // 2
                if int(_f(np.array([m]), p, lo)[0]) > int(th[j]): Rb = m - 1
                else: Lb = m
            U[j] = Lb
    ne = L <= U
    L = L[ne]; U = U[ne]
    sg = np.concatenate([np.ones(len(L), dtype=np.int64), -np.ones(len(L), dtype=np.int64)])
    return sg, np.concatenate([L, L]), np.concatenate([U, U])
def back_n(nodes, p, lo, n, cap=4000000):
    """nodes = [(sg, i1, i2)] または配列 3 つ組。n 歩戻る。節が cap を超えたら None"""
    if isinstance(nodes, list):
        sg = np.array([x[0] for x in nodes], dtype=np.int64); t1 = np.array([x[1] for x in nodes], dtype=np.int64); t2 = np.array([x[2] for x in nodes], dtype=np.int64)
    else: sg, t1, t2 = nodes
    for _ in range(n):
        if len(sg) == 0: break
        sg, t1, t2 = back_nodes(sg, t1, t2, p, lo)
        if len(sg) > cap: return None
    return sg, t1, t2
if __name__ == '__main__':
    import sys, os, random, time
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'gcc6')); import gcc
    random.seed(7); nt = 0; nok = 0; sizes = []; nover = 0; t0 = time.time()
    while nt < 400:
        kb = bytes(random.choice(b'abcdefghijklmnopqrstuvwxyz0123456789 ') for _ in range(random.randint(5, 40)))
        st = gcc.key_schedule(kb)
        if st.mode != 1: continue
        nt += 1; lo = st.b0 | (st.b1 << 8); p = st.r; yinit = gcc.set_low16(st.k780 * 2.0 - 1.0, st.b0, st.b1)
        y = yinit
        for _ in range(st.warm + 1): y = G(y, p, lo)
        s2 = gcc.key_schedule(kb); gcc.step(s2); assert s2.y == y
        sg, k = tkey(y); r = back_n([(sg, k, k)], p, lo, st.warm + 1)
        if r is None: nover += 1; continue
        sgi, ki = tkey(yinit); hit = bool(((r[0] == sgi) & (r[1] <= ki) & (ki <= r[2])).any())
        nok += hit; sizes.append((len(r[0]), int((r[2] - r[1] + 1).sum()), st.warm))
        if not hit: print('  見逃し: 鍵 %r warm %d' % (kb, st.warm))
    print('写像 1 の逆向き: %d 本中 %d 本で真の初期値が候補に入った（節が上限を超えたもの %d 本）、%.0f 秒' % (nt, nok, nover, time.time() - t0))
    sizes.sort(); pts = sorted(s[1] for s in sizes)
    print('候補の節の数（最小・中央・90%%・最大）: %d・%d・%d・%d。点の数（中央・90%%・最大）: %d・%d・%d' % (sizes[0][0], sizes[len(sizes) // 2][0], sizes[int(len(sizes) * 0.9)][0], sizes[-1][0], pts[len(pts) // 2], pts[int(len(pts) * 0.9)], pts[-1]))
