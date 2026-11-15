"""fulldec.py — 2 重の連鎖のフラグと 2 本の平文を、全 256 語で同時に探す（鍵なし。2026-11-12）
  モデル: A = E_k(E_v(M_A))、E = E_k(E_v(M_E))（同じ 2 つの鍵）。X = A ⊕ E として、語 n ごとに
      PX_n = X_n ⊕ f_n·Rout、 MX_n = PX_n ⊕ f'_n·Rin、 f_n = 1 なら Rout ← X_n、 f'_n = 1 なら Rin ← PX_n   （最初は Rout = Rin = 0）
      MX_n = M_A ⊕ M_E の語 n
  f・f'（各語 4 通り）と 2 本の文字を、英文の言語モデルの点でビーム探索する。8bit の文字は「高位」という 1 つの記号として罰点つきで許す。
  使う情報は X だけ（B は使わない）。A = B の語を連鎖なし（f = f' = 0）に固定するのは -same を付けたときだけ。
  python3 -I fulldec.py [-w 幅] [-same] [-lo 語] [-hi 語] [-o 出力]"""
import sys, os, math, time, pickle
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
import lm as LMmod
from lm import ALPHA, AIDX
HI = 0x80; LOG_HI = float(os.environ.get("FD_HI", "-10")); LOG_HI_CONT = float(os.environ.get("FD_HIC", "-6"))     # 8bit の文字の罰点（直前も 8bit なら軽い: 2 バイト文字）
LF = {(0, 0): math.log(0.07), (1, 0): math.log(0.24), (0, 1): math.log(0.15), (1, 1): math.log(0.54)}   # フラグの事前確率（窓の中の割合から）
def args():
    a = sys.argv[1:]; o = dict(w=400, same=False, lo=0, hi=255, out='fulldec_out.pkl')
    i = 0
    while i < len(a):
        if a[i] == '-w': o['w'] = int(a[i + 1]); i += 2
        elif a[i] == '-same': o['same'] = True; i += 1
        elif a[i] == '-lo': o['lo'] = int(a[i + 1]); i += 2
        elif a[i] == '-hi': o['hi'] = int(a[i + 1]); i += 2
        elif a[i] == '-o': o['out'] = a[i + 1]; i += 2
        else: raise SystemExit('?' + a[i])
    return o
def xor4(a, b): return bytes(x ^ y for x, y in zip(a, b))
Z4 = bytes(4)
def run(Xw, lm, W, lo, hi, fixed, verbose=True):
    K = lm.order - 1
    # ビームの要素: key = (ctxA, ctxB, jout, Rin, mx, flags) -> (score, back)
    beam = {(b'', b'', -1, Z4, None, None): (0.0, None)}
    trace = []    # 位置ごとの {key: (prevkey, a, b)}
    t0 = time.time()
    for n in range(lo, hi + 1):
        x = Xw[n]
        # 語の頭: フラグを分ける
        nb = {}
        for key0, (sc, _) in beam.items():
            ca, cb, jout, rin = key0[:4]
            rout = Xw[jout] if jout >= 0 else Z4
            for f in (0, 1):
                for g in (0, 1):
                    if n in fixed and (f, g) != fixed[n]: continue
                    px = xor4(x, rout) if f else x
                    mx = xor4(px, rin) if g else px
                    key = (ca, cb, n if f else jout, px if g else rin, mx, (f, g))
                    s = sc + LF[(f, g)]
                    cur = nb.get(key)
                    if cur is None or s > cur[0]: nb[key] = (s, key0)
        trace.append({k: (v[1], None, None) for k, v in nb.items()})
        beam = {k: (v[0], None) for k, v in nb.items()}
        for c in range(4):
            nb = {}
            for key, (sc, _) in beam.items():
                ca, cb, jout, rin, mx, fl = key; z = mx[c]
                va = lm.logp_vec(ca); vb = lm.logp_vec(cb)
                if z & 0x80:
                    # どちらか一方が 8bit
                    pa = LOG_HI_CONT if ca[-1:] == b'\x80' else LOG_HI; pb = LOG_HI_CONT if cb[-1:] == b'\x80' else LOG_HI
                    for ai, a in enumerate(ALPHA):
                        s = sc + va[ai] + pb; k2 = ((ca + bytes([a]))[-K:], (cb + b'\x80')[-K:], jout, rin, mx, fl)
                        cur = nb.get(k2)
                        if cur is None or s > cur[0]: nb[k2] = (s, key, a, a ^ z)
                        s = sc + pa + vb[ai]; k2 = ((ca + b'\x80')[-K:], (cb + bytes([a]))[-K:], jout, rin, mx, fl)
                        cur = nb.get(k2)
                        if cur is None or s > cur[0]: nb[k2] = (s, key, a ^ z, a)
                else:
                    for ai, a in enumerate(ALPHA):
                        b = a ^ z; bi = AIDX.get(b)
                        if bi is None: continue
                        s = sc + va[ai] + vb[bi]; k2 = ((ca + bytes([a]))[-K:], (cb + bytes([b]))[-K:], jout, rin, mx, fl)
                        cur = nb.get(k2)
                        if cur is None or s > cur[0]: nb[k2] = (s, key, a, b)
                    # 両方が 8bit（まれ）
                    if ca[-1:] == b'\x80' and cb[-1:] == b'\x80':
                        s = sc + 2 * LOG_HI_CONT; k2 = ((ca + b'\x80')[-K:], (cb + b'\x80')[-K:], jout, rin, mx, fl)
                        cur = nb.get(k2)
                        if cur is None or s > cur[0]: nb[k2] = (s, key, 0x80, 0x80 ^ z)
            items = sorted(nb.items(), key=lambda kv: -kv[1][0])[:W]
            trace.append({k: (v[1], v[2], v[3]) for k, v in items})
            beam = {k: (v[0], None) for k, v in items}
        if verbose and (n - lo) % 16 == 15: print('   語 %d まで、最良 %.1f、%.0f 秒' % (n, max(v[0] for v in beam.values()), time.time() - t0), flush=True)
    # 逆にたどる
    key = max(beam.items(), key=lambda kv: kv[1][0])[0]; best = beam[key][0]
    a = bytearray(); b = bytearray(); flags = {}; mxs = {}
    pos = len(trace) - 1
    for n in range(hi, lo - 1, -1):
        flags[n] = key[5]; mxs[n] = key[4]
        for c in range(3, -1, -1):
            prev, ca_, cb_ = trace[pos][key]; a.append(ca_); b.append(cb_); key = prev; pos -= 1
        prev, _, _ = trace[pos][key]; key = prev; pos -= 1
    return bytes(reversed(a)), bytes(reversed(b)), flags, mxs, best
def show(t): return ''.join(chr(c) if 0x20 <= c < 0x7f else ('¶' if c == 0x0a else '↵' if c == 0x0d else '■' if c >= 0x80 else '·') for c in t)
def main():
    o = args()
    A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
    Xw = [xor4(A[4 * n:4 * n + 4], E[4 * n:4 * n + 4]) for n in range(256)]
    same = [n for n in range(256) if A[4 * n:4 * n + 4] == B[4 * n:4 * n + 4]]
    fixed = {n: (0, 0) for n in same} if o['same'] else {}
    lm = LMmod.load(os.path.join(D, os.environ.get('TTP_LM', 'lm_dom.pkl')))
    a, b, flags, mxs, best = run(Xw, lm, o['w'], o['lo'], o['hi'], fixed)
    print('幅 %d、語 %d〜%d、点 %.1f' % (o['w'], o['lo'], o['hi'], best))
    lo = o['lo']
    for k in range(0, len(a), 64):
        print('  %4d a: %s' % (4 * lo + k, show(a[k:k + 64]))); print('       b: %s' % show(b[k:k + 64]))
    fo = ''.join(str(flags[n][0]) for n in range(lo, o['hi'] + 1)); fi = ''.join(str(flags[n][1]) for n in range(lo, o['hi'] + 1))
    print('  外側のフラグ: ' + fo); print('  内側のフラグ: ' + fi)
    print('  A = B の語のフラグ: ' + ' '.join('%d:%d%d' % (n, flags[n][0], flags[n][1]) for n in same if lo <= n <= o['hi']))
    print('  外 0・内 0 と決まった語: %s' % [n for n in range(lo, o['hi'] + 1) if flags[n] == (0, 0)])
    pickle.dump(dict(a=a, b=b, flags=flags, mxs=mxs, lo=lo, hi=o['hi'], best=best), open(os.path.join(D, o['out']), 'wb'))
if __name__ == '__main__': main()
