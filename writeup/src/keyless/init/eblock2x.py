"""eblock2x.py — eblock2.py で通った候補から、E の 2 ブロック目の 2 つの層の状態を倍精度で厳密に決め、314 バイトを復号する。2026-11-12
  python3 -I -B eblock2x.py"""
import sys, os, pickle, struct, time, math, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
import eblock2 as E2M
V = E2M.V; u2d, d2u, setlow = V.u2d, V.d2u, V.setlow
Ew2, Mw, RK0, RV0, NW2, TAIL = E2M.Ew2, E2M.Mw, E2M.RK0, E2M.RV0, E2M.NW2, E2M.TAIL
r = pickle.load(open(os.path.join(D, 'eblock2_found.pkl'), 'rb')); found, kp, vp, p48c = r['found'], r['kp'], r['vp'], r['p48c']
print('通った v の候補: %d。dp の内訳 %s、lo の範囲 %04x〜%04x（%d 個）、idx の種類 %d' % (len(found), dict(collections.Counter(f[3] for f in found)), min(f[0] for f in found), max(f[0] for f in found), len(set(f[0] for f in found)), len(set(f[2] for f in found))))
hyp = found[0][4]; assert all(f[4] == hyp for f in found)
NP = len(Mw)                                                     # 平文を予想している語の数
# v の候補を、全 79 語の鍵ストリームが同じものでまとめる
groups = collections.OrderedDict()
for lo, Sv, idx, dp, _ in found:
    yinit = u2d((idx << 16) | lo); p2 = u2d(((p48c + dp) << 16) | lo); ks, fl = E2M.vstream2(yinit, p2, lo, NW2 + 1)
    groups.setdefault((tuple(ks), tuple(fl)), []).append((lo, idx, dp))
print('v の鍵ストリーム（79 語）の種類: %d' % len(groups))
for (ks, fl), mem in groups.items(): print('   %d 個の候補（lo %04x〜%04x、dp %s、idx %s）' % (len(mem), min(m[0] for m in mem), max(m[0] for m in mem), sorted(set(m[2] for m in mem)), sorted(set(m[1] for m in mem))))
def kword(v):
    b = struct.pack('<d', v); return (b[4] << 24) | (b[2] << 16) | (b[3] << 8) | b[5], (0 if b[6] > 0xE6 else 1)
def step2(Q, A, C, lo, u, v):
    a = Q * u; b = a * u; nu = A * u - a * v; nv = C * v + b
    return setlow(nu, lo), setlow(nv, lo)
def vbits(Kw, b7, b6, lo):
    b4 = (Kw >> 24) & 0xFF; b2 = (Kw >> 16) & 0xFF; b3 = (Kw >> 8) & 0xFF; b5 = Kw & 0xFF
    return u2d((b7 << 56) | (b6 << 48) | (b5 << 40) | (b4 << 32) | (b3 << 24) | (b2 << 16) | lo)
Q48c, A48c, C48c = d2u(kp['Q']) >> 16, d2u(kp['A']) >> 16, d2u(kp['C']) >> 16
sols = []; t0 = time.time()
for gi, ((ksv, flv), mem) in enumerate(groups.items()):
    cin, _ = E2M.inner_ct(Mw, list(ksv[:NP]), list(flv[:NP]), RV0); X = [Ew2[n] ^ cin[n] for n in range(NP)]
    # 語 0・1 の出力（上位 48 ビット）は仮説から
    f0 = 0 if hyp[0][1] > 0xE6 else 1; K0 = X[0] ^ (RK0 if f0 else 0); reg1 = Ew2[0] if f0 else RK0
    f1 = 0 if hyp[1][1] > 0xE6 else 1; K1 = X[1] ^ (reg1 if f1 else 0)
    ntry = 0
    for lo in range(0x7e00, 0x7f80):
        v0 = vbits(K0, hyp[0][0], hyp[0][1], lo); v1 = vbits(K1, hyp[1][0], hyp[1][1], lo); v1lo = u2d(d2u(v1) & ~0xFFFF); v1hi = u2d(d2u(v1) | 0xFFFF)
        for dQ in range(-3, 4):
            Q = u2d(((Q48c + dQ) << 16) | lo)
            for dC in range(-3, 4):
                C = u2d(((C48c + dC) << 16) | lo)
                blo = v1lo - C * v0; bhi = v1hi - C * v0
                if bhi <= 0: continue
                ia = (d2u(math.sqrt(max(blo, 0.0) / Q)) >> 16) - 3; ib = (d2u(math.sqrt(bhi / Q)) >> 16) + 3
                for iu in range(ia, ib + 1):
                    for sg in (1.0, -1.0):
                        u0 = sg * u2d((iu << 16) | lo)
                        a = Q * u0; b = a * u0; nv = C * v0 + b
                        if (d2u(nv) >> 16) != (d2u(v1) >> 16): continue
                        for dA in range(-3, 4):
                            A = u2d(((A48c + dA) << 16) | lo); ntry += 1
                            # 厳密に前へ回して、予想している語の全部で鍵ストリームを照合
                            u, v = u0, v0; reg = RK0; ok = True; nmatch = 0
                            for n in range(NP):
                                if n: u, v = step2(Q, A, C, lo, u, v)
                                kw, f = kword(v); want = X[n] ^ (reg if f else 0)
                                if kw != want: ok = False; break
                                nmatch += 1
                                if f: reg = Ew2[n]
                            if ok: sols.append(dict(group=gi, lo=lo, dQ=dQ, dA=dA, dC=dC, Q=Q, A=A, C=C, u0=u0, v0=v0)); print('  ★ 厳密に一致（%d 語全部）: v の組 %d、lo_k = %04x、dQ %+d dA %+d dC %+d、u = %.17g v = %.17g' % (NP, gi, lo, dQ, dA, dC, u0, v0), flush=True)
    print('v の組 %d: 試した (状態, 定数) %d 通り、%.0f 秒' % (gi, ntry, time.time() - t0), flush=True)
print('厳密な解: %d 個。lo_k の範囲 %s' % (len(sols), ('%04x〜%04x' % (min(s['lo'] for s in sols), max(s['lo'] for s in sols))) if sols else '-'))
pickle.dump(dict(sols=sols, groups=[(list(k[0]), list(k[1]), m) for k, m in groups.items()], hyp=hyp), open(os.path.join(D, 'eblock2_exact.pkl'), 'wb'))
