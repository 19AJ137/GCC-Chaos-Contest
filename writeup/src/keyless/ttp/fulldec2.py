"""fulldec2.py — 2 重の連鎖のフラグと 2 本の平文を全 256 語で同時に探す（鍵なし。fulldec.py の改良。2026-11-12）
  モデル: A = E_k(E_v(M_A))、E = E_k(E_v(M_E))（同じ 2 つの鍵）。X = A ⊕ E として、語 n ごとに
      PX_n = X_n ⊕ f_n·Rout、 MX_n = PX_n ⊕ f'_n·Rin、 f_n = 1 なら Rout ← X_n、 f'_n = 1 なら Rin ← PX_n   （最初は Rout = Rin = 0）
      MX_n = M_A ⊕ M_E の語 n
  M_E は英文（7bit）だけ、M_A は英文か Shift_JIS の日本語、として、別々の言語モデルで点を付ける（lm2.py）。
  使う情報は A ⊕ E だけ（B は使わない）。オプションで語のフラグや M_E の文字を固定できる。
  python3 -I fulldec2.py [-w 幅] [-lo 語] [-hi 語] [-o 出力] [-fix FILE] [-crib FILE] [-init FILE]
     -fix FILE:  「語 f f'」の行（フラグを固定）
     -crib FILE: 「バイト位置 文字列（Python のリテラル）」の行（M_E のその位置からの文字を固定）
     -init FILE: 前の実行の出力（pkl）。-lo の語の直前の状態（レジスタと文脈）をそこから取る"""
import sys, os, math, time, pickle, ast
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
import lm2
from lm2 import A7, A7IDX
LF = {(0, 0): math.log(0.07), (1, 0): math.log(0.24), (0, 1): math.log(0.15), (1, 1): math.log(0.54)}
def args():
    a = sys.argv[1:]; o = dict(w=600, lo=0, hi=255, out='fd2_out.pkl', fix=None, crib=None, criba=None, init=None, swap=False)
    i = 0
    while i < len(a):
        k = a[i]
        if k in ('-w', '-lo', '-hi'): o[k[1:]] = int(a[i + 1]); i += 2
        elif k == '-o': o['out'] = a[i + 1]; i += 2
        elif k in ('-fix', '-crib', '-criba', '-init'): o[k[1:]] = a[i + 1]; i += 2
        else: raise SystemExit('?' + k)
    return o
def xor4(a, b): return bytes(x ^ y for x, y in zip(a, b))
Z4 = bytes(4)
def run(Xw, le, la, W, lo, hi, fixed, crib, init=None, verbose=True, criba=None):
    criba = criba or {}
    K = le.order - 1
    st0 = init or (b'', b'', -1, Z4)
    beam = {(st0[0], st0[1], st0[2], st0[3], None, None): (0.0, None)}
    trace = []; t0 = time.time()
    for n in range(lo, hi + 1):
        x = Xw[n]; nb = {}
        for key0, (sc, _) in beam.items():
            ce, ca, jout, rin = key0[:4]
            rout = Xw[jout] if jout >= 0 else Z4
            for f in (0, 1):
                for g in (0, 1):
                    if n in fixed and (f, g) != fixed[n]: continue
                    px = xor4(x, rout) if f else x
                    mx = xor4(px, rin) if g else px
                    key = (ce, ca, n if f else jout, px if g else rin, mx, (f, g))
                    s = sc + LF[(f, g)]; cur = nb.get(key)
                    if cur is None or s > cur[0]: nb[key] = (s, key0)
        trace.append({k: (v[1], None, None) for k, v in nb.items()})
        beam = {k: (v[0], None) for k, v in nb.items()}
        for c in range(4):
            pos = 4 * n + c; nb = {}; cr = crib.get(pos); cra = criba.get(pos)
            for key, (sc, _) in beam.items():
                ce, ca, jout, rin, mx, fl = key; z = mx[c]
                ve = le.logp_vec7(ce)
                if z & 0x80:
                    for ei, e in enumerate(A7):
                        if cr is not None and e != cr: continue
                        a = e ^ z
                        if cra is not None and a != cra: continue
                        s = sc + ve[ei] + la.logp1(ca, a); k2 = ((ce + bytes([e]))[-K:], (ca + bytes([a]))[-K:], jout, rin, mx, fl)
                        cur = nb.get(k2)
                        if cur is None or s > cur[0]: nb[k2] = (s, key, e, a)
                else:
                    va = la.logp_vec7(ca)
                    for ei, e in enumerate(A7):
                        if cr is not None and e != cr: continue
                        a = e ^ z; ai = A7IDX.get(a)
                        if ai is None: continue
                        if cra is not None and a != cra: continue
                        s = sc + ve[ei] + va[ai]; k2 = ((ce + bytes([e]))[-K:], (ca + bytes([a]))[-K:], jout, rin, mx, fl)
                        cur = nb.get(k2)
                        if cur is None or s > cur[0]: nb[k2] = (s, key, e, a)
            if not nb:
                raise SystemExit('バイト %d で候補がなくなった（固定した文字と合わない）' % pos)
            items = sorted(nb.items(), key=lambda kv: -kv[1][0])[:W]
            trace.append({k: (v[1], v[2], v[3]) for k, v in items})
            beam = {k: (v[0], None) for k, v in items}
        if verbose and (n - lo) % 32 == 31: print('   語 %d まで、最良 %.1f、%.0f 秒' % (n, max(v[0] for v in beam.values()), time.time() - t0), flush=True)
    key = max(beam.items(), key=lambda kv: kv[1][0])[0]; best = beam[key][0]; endstate = key[:4]
    e = bytearray(); a = bytearray(); flags = {}; mxs = {}; pos = len(trace) - 1
    for n in range(hi, lo - 1, -1):
        flags[n] = key[5]; mxs[n] = key[4]
        for c in range(3, -1, -1):
            prev, ce_, ca_ = trace[pos][key]; e.append(ce_); a.append(ca_); key = prev; pos -= 1
        prev, _, _ = trace[pos][key]; key = prev; pos -= 1
    return bytes(reversed(e)), bytes(reversed(a)), flags, mxs, best, endstate
def show(t):
    out = []; i = 0
    while i < len(t):
        c = t[i]
        if 0x20 <= c < 0x7f: out.append(chr(c)); i += 1
        elif c == 0x0a: out.append('¶'); i += 1
        elif c == 0x0d: out.append('↵'); i += 1
        elif (0x81 <= c <= 0x9f or 0xe0 <= c <= 0xfc) and i + 1 < len(t):
            try: out.append(t[i:i + 2].decode('cp932')); i += 2
            except Exception: out.append('■'); i += 1
        elif 0xa1 <= c <= 0xdf: out.append(bytes([c]).decode('cp932')); i += 1
        else: out.append('■'); i += 1
    return ''.join(out)
def main():
    o = args()
    A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
    Xw = [xor4(A[4 * n:4 * n + 4], E[4 * n:4 * n + 4]) for n in range(256)]
    same = [n for n in range(256) if A[4 * n:4 * n + 4] == B[4 * n:4 * n + 4]]
    fixed = {}; crib = {}; criba = {}
    def readcrib(fn, d):
        for ln in open(fn, encoding='utf-8'):
            if ln.strip() and not ln.startswith('#'):
                p, t = ln.split(None, 1); t = ast.literal_eval(t.strip())
                if isinstance(t, str): t = t.encode('cp932')
                for i, ch in enumerate(t): d[int(p) + i] = ch
    if o['criba']: readcrib(o['criba'], criba)
    if o['fix']:
        for ln in open(o['fix']):
            p = ln.split()
            if len(p) >= 3 and not ln.startswith('#'): fixed[int(p[0])] = (int(p[1]), int(p[2]))
    if o['crib']:
        for ln in open(o['crib'], encoding='utf-8'):
            if ln.strip() and not ln.startswith('#'):
                p, s = ln.split(None, 1); s = ast.literal_eval(s.strip())
                if isinstance(s, str): s = s.encode('ascii')
                for i, ch in enumerate(s): crib[int(p) + i] = ch
    init = None
    if o['init']:
        d = pickle.load(open(o['init'], 'rb')); init = d['states'][o['lo']]
    le = lm2.load('lm2_e.pkl'); la = lm2.load('lm2_a.pkl')
    e, a, flags, mxs, best, endstate = run(Xw, le, la, o['w'], o['lo'], o['hi'], fixed, crib, init, criba=criba)
    lo = o['lo']; print('幅 %d、語 %d〜%d、点 %.1f' % (o['w'], lo, o['hi'], best))
    for k in range(0, len(a), 64):
        print('  %4d E: %s' % (4 * lo + k, show(e[k:k + 64]))); print('       A: %s' % show(a[k:k + 64]))
    rng = range(lo, o['hi'] + 1)
    print('  外側のフラグ: ' + ''.join(str(flags[n][0]) for n in rng)); print('  内側のフラグ: ' + ''.join(str(flags[n][1]) for n in rng))
    print('  A = B の語のフラグ: ' + ' '.join('%d:%d%d' % (n, flags[n][0], flags[n][1]) for n in same if n in flags))
    print('  外 0・内 0 と決まった語: %s' % [n for n in rng if flags[n] == (0, 0)])
    # 語ごとの直前の状態を作り直す（-init 用）
    states = {}; jout = init[2] if init else -1; rin = init[3] if init else Z4
    for n in rng:
        ce = (init[0] if init else b'') + e[:4 * (n - lo)]; ca = (init[1] if init else b'') + a[:4 * (n - lo)]
        states[n] = (ce[-(le.order - 1):], ca[-(le.order - 1):], jout, rin)
        f, g = flags[n]; x = Xw[n]; px = xor4(x, Xw[jout] if jout >= 0 else Z4) if f else x
        if f: jout = n
        if g: rin = px
    states[o['hi'] + 1] = endstate
    pickle.dump(dict(e=e, a=a, flags=flags, mxs=mxs, lo=lo, hi=o['hi'], best=best, states=states), open(os.path.join(D, o['out']), 'wb'))
if __name__ == '__main__': main()
