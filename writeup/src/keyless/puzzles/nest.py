"""nest.py — 「8bit の列」は、平文の性質ではなく、内側にもう 1 回かかっている暗号化（別の鍵 v）の連鎖フラグではないか、の検査（鍵なし、2026-11-12）
  仮説: A = E_k(E_v(M_A))、E = E_k(E_v(M_E))（外側の鍵 k、内側の鍵 v は共通。M は本当の平文）。このとき、外側の連鎖を外した平文の XOR は
     PX_n = E_v(M_A)_n ⊕ E_v(M_E)_n = MX_n ⊕ f'_n·PX_{p'(n)}       （MX = M_A⊕M_E、f' は v の連鎖フラグ、p'(n) は n の直前の v の連鎖語）
  で、M_A・M_E がどちらも 7bit なら、最上位ビットが 1 の列は「v の連鎖ありの語で 1、連鎖なしの語で 0」になる（= 8bit の列）。
  そこで、8bit の列の最上位ビットが 1 の行（型 8）を v の連鎖あり、0 の行（型 7）を連鎖なしと読み、内側の連鎖も外す:
     MX_n = PX_n（型 7）、 MX_n = PX_n ⊕ PX_{直前の型 8 の行}（型 8）
  仮説が正しければ、MX は 2 本の文章の XOR らしくなる（小さい値が多い）。8bit の列が平文の性質なら、型 8 の行の PX がすでに平文の XOR で、
  もう 1 回 XOR を取った MX のほうが乱雑になる。
  python3 -I nest.py [フラグ（既定 v2）]"""
import sys, os, collections, math, random
D = os.path.dirname(os.path.abspath(__file__))
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
def W(d, i): return d[4 * i:4 * i + 4]
def X(*a):
    r = bytes(4)
    for x in a: r = bytes(p ^ q for p, q in zip(r, x))
    return r
def read_flags(fn):
    fl = {}
    for line in open(fn):
        line = line.strip()
        if not line: continue
        tag, body = line.split(':', 1)
        if tag == 'z':
            for w in body.split(','):
                if w.strip(): fl[int(w)] = 0
        else:
            b = int(tag)
            for k, c in enumerate(body.strip()): fl[b + k] = 1 if c == '1' else 0
    return fl
def msb(b): return ''.join('1' if c & 0x80 else '0' for c in b)
pr = lambda b: ''.join(chr(c) if 0x20 <= c < 0x7f else '.' for c in b)
def main():
    name = sys.argv[1] if len(sys.argv) > 1 else 'v2'
    A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
    fl = read_flags(os.path.join(D, 'real_flags_%s.txt' % name)); same = set(n for n in range(256) if W(A, n) == W(B, n))
    allst = collections.defaultdict(list)
    for wn, lo, hi, col in (('W2', 64, 119, 2), ('W1', 230, 255, 0)):
        hi = max(n for n in range(lo, hi + 1) if n in fl and all(k in fl for k in range(lo, n + 1)))
        rows = []; prevc = None
        for n in range(lo, hi + 1):
            if fl[n] == 1:
                px = X(W(A, n), W(E, n), W(A, prevc), W(E, prevc)) if prevc is not None else None
                prevc = n
            else: px = X(W(A, n), W(E, n))
            rows.append((n, fl[n], px))
        print('%s（語 %d〜%d、8bit の列は位置 %d）' % (wn, lo, hi, col))
        print('   語  外側  型  PX（A・E の平文の XOR）   内側の連鎖も外した MX      A=B')
        prev8 = None; types = []
        for n, f, px in rows:
            if px is None: print('  %3d  %s  —' % (n, '連' if f else '・')); types.append(None); continue
            t = 8 if px[col] & 0x80 else 7; types.append(t)
            if t == 7: mx = px
            else:
                mx = X(px, prev8) if prev8 is not None else None
                prev8 = px
            print('  %3d  %s   %d  %s %s %-4s   %s' % (n, '連' if f else '・', t, px.hex(), msb(px), pr(px), ('%s %s %-4s' % (mx.hex(), msb(mx), pr(mx))) if mx else '—（窓の最初の型 8）') + ('   一致' if n in same else ''))
            allst[(wn, 'PX', t)].append(px)
            if mx is not None: allst[(wn, 'MX', t)].append(mx)
        tt = [t for t in types if t]; print('   型の並び（7 = v の連鎖なし、8 = 連鎖あり）: ' + ''.join('.' if t == 7 else '1' if t == 8 else '?' for t in types) + '   型 7 の割合 %d / %d' % (tt.count(7), len(tt)))
    def summ(ws, cols=None):
        out = []
        for lab, key in (('型 7 の PX', ('PX', 7)), ('型 8 の PX', ('PX', 8)), ('型 8 の MX（内側の連鎖を外す）', ('MX', 8))):
            by = []
            for wn, col in ws:
                J = [c for c in range(4) if (cols == 'J' and c != col) or (cols == 'col' and c == col) or cols is None]
                for w in allst[(wn,) + key]: by += [w[c] & 0x7f for c in J]
            n = len(by)
            if not n: continue
            out.append('%s: %d バイト、0x20 未満 %d（%.0f%%）、0x10 未満 %d（%.0f%%）、0 のバイト %d' % (lab, n, sum(v < 0x20 for v in by), 100 * sum(v < 0x20 for v in by) / n, sum(v < 0x10 for v in by), 100 * sum(v < 0x10 for v in by) / n, sum(v == 0 for v in by)))
        return out
    print('\n下位 7 ビットの値の小ささ（7bit の一様な乱数なら 0x20 未満は 25%、0x10 未満は 12.5%）')
    for lab, ws in (('W2', [('W2', 2)]), ('W1', [('W1', 0)]), ('両方', [('W2', 2), ('W1', 0)])):
        for cl, cn in ((None, '4 列とも'), ('J', '一致列（3 列）'), ('col', '8bit の列')):
            print('  %s・%s' % (lab, cn))
            for l in summ(ws, cl): print('      ' + l)
if __name__ == '__main__': main()
