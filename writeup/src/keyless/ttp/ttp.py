"""ttp.py — 2 本の文章の XOR（z = a ⊕ b）から a・b を言語モデルで探す（ビーム探索）。2026-11-12
  python3 -I ttp.py test            合成（自作の英文 2 本）で確かめる
  python3 -I ttp.py FILE [幅]       FILE（バイト列）をほどく"""
import sys, os, math, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lm as LMmod
from lm import ALPHA, AIDX
def decode(z, lm, width=3000, unk=None, verbose=False):
    """z: bytes。unk: z が分からない位置の集合（そこは a・b を独立に選ぶ…重いので、a だけ選んで b は文脈から最良を取る）"""
    K = lm.order - 1
    beam = {(b'', b''): (0.0, None)}   # (ctxA, ctxB) -> (score, backpointer)
    hist = []
    valid = set(ALPHA)
    for i, zi in enumerate(z):
        nb = {}
        for (ca, cb), (sc, _) in beam.items():
            va = lm.logp_vec(ca); vb = lm.logp_vec(cb)
            for ai, a in enumerate(ALPHA):
                b = a ^ zi
                bi = AIDX.get(b)
                if bi is None: continue
                s = sc + va[ai] + vb[bi]
                key = ((ca + bytes([a]))[-K:], (cb + bytes([b]))[-K:])
                cur = nb.get(key)
                if cur is None or s > cur[0]: nb[key] = (s, (ca, cb), a)
        if not nb: raise ValueError('位置 %d で候補がなくなった（z=%02x）' % (i, zi))
        items = sorted(nb.items(), key=lambda kv: -kv[1][0])[:width]
        hist.append({k: (v[1], v[2]) for k, v in items})
        beam = {k: (v[0], None) for k, v in items}
        if verbose and i % 40 == 39: print('   位置 %d、最良 %.1f' % (i + 1, items[0][1][0]), flush=True)
    # 逆にたどる
    best = max(beam.items(), key=lambda kv: kv[1][0]); key = best[0]; a = bytearray()
    for i in range(len(z) - 1, -1, -1):
        prev, ch = hist[i][key]; a.append(ch); key = prev
    a = bytes(reversed(a)); b = bytes(x ^ y for x, y in zip(a, z))
    return a, b, best[1][0]
def show(t): return ''.join(chr(c) if 0x20 <= c < 0x7f else ('¶' if c == 0x0a else '↵' if c == 0x0d else '·') for c in t)
def main():
    lm = LMmod.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), os.environ.get('TTP_LM', 'lm_dom.pkl')))
    if sys.argv[1] == 'test':
        t1 = (b"We are pleased to announce a prize contest for our new cryptosystem.\r\nAnyone who can break the cipher and recover the original message will\r\nreceive one million yen. Please send your answer together with the key\r\nto the address below before the end of December.\r\nGood luck to all of you!")
        t2 = (b"This file was encrypted with the chaos cipher developed by our institute.\r\nThe algorithm uses a nonlinear map and floating point numbers to make\r\na long sequence of random numbers. If you find any weakness in it,\r\nwe would like to hear from you. Our telephone number is 03-1234-5678.")
        n = min(len(t1), len(t2), 248); z = bytes(x ^ y for x, y in zip(t1[:n], t2[:n]))
        for w in (300, 2000):
            t0 = time.time(); a, b, sc = decode(z, lm, w)
            okA = sum(1 for i in range(n) if a[i] in (t1[i], t2[i])); print('幅 %d: %.0f 秒、点 %.1f（正解の点 %.1f）、正解と合う文字 %d / %d' % (w, time.time() - t0, sc, lm.score(t1[:n]) + lm.score(t2[:n]), okA, n))
            print('  a: ' + show(a)); print('  b: ' + show(b))
        return
    z = open(sys.argv[1], 'rb').read(); w = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    t0 = time.time(); a, b, sc = decode(z, lm, w, verbose=True)
    print('幅 %d: %.0f 秒、点 %.1f（1 文字あたり %.2f ビット × 2 本）' % (w, time.time() - t0, sc, -sc / len(z) / 2 / math.log(2)))
    for k in range(0, len(z), 64):
        print('  %4d a: %s' % (k, show(a[k:k + 64]))); print('       b: %s' % show(b[k:k + 64]))
if __name__ == '__main__': main()
