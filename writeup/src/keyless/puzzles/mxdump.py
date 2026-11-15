"""mxdump.py — 読み A で内側まで外した MX（= 2 本の平文の XOR と読む値）の並びを見る（鍵なし、2026-11-12）"""
import sys, os
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
from nest4 import decode, XW, same
def stream(lo, hi, oc, ic):
    R = decode(lo, hi, oc, ic); out = {}
    for r in R:
        if r['mx'] is not None: out[r['n']] = r['mx']
    return out
def view(b, k): 
    c = b ^ k
    return chr(c) if 0x20 <= c < 0x7f else '·'
def main():
    for name, lo, hi, oc, ic in (('W2', 55, 118, 3, 2), ('W1', 229, 255, 1, 0)):
        M = stream(lo, hi, oc, ic); ns = sorted(M)
        print('%s: 語 %d〜%d（バイト %d〜%d）' % (name, ns[0], ns[-1], 4 * ns[0], 4 * ns[-1] + 3))
        print('  語（バイト）   MX            相手が空白なら   相手が CR LF なら（位置 2・3）  相手が同じ種類の英字なら大小が違う箇所')
        for n in ns:
            m = M[n]
            print('  %3d（%4d）  %s   %s   %s%s   %s%s' % (n, 4 * n, ' '.join('%02x' % b for b in m), ''.join(view(b, 0x20) if b & 0x40 else ('=' if b == 0 else '.') for b in m),
                  view(m[2], 0x0d), view(m[3], 0x0a), ''.join('^' if (b & 0x60) == 0x20 else ('#' if (b & 0x60) == 0x60 else ' ') for b in m), '   A=B' if n in same else ''))
        by = b''.join(M[n] for n in ns)
        print('  値の区分: 0x00〜1f %d、0x20〜3f %d、0x40〜5f %d、0x60〜7f %d（全 %d バイト）' % (sum(b < 0x20 for b in by), sum(0x20 <= b < 0x40 for b in by), sum(0x40 <= b < 0x60 for b in by), sum(0x60 <= b < 0x80 for b in by), len(by)))
        for c in range(4):
            col = [M[n][c] for n in ns]
            print('   位置 %d: 0x00〜1f %d、0x20〜3f %d、0x40〜5f %d、0x60〜7f %d' % (c, sum(b < 0x20 for b in col), sum(0x20 <= b < 0x40 for b in col), sum(0x40 <= b < 0x60 for b in col), sum(0x60 <= b < 0x80 for b in col)))
        open(os.path.join(D, 'mx_%s.bin' % name), 'wb').write(by); open(os.path.join(D, 'mx_%s.txt' % name), 'w').write('%d %d\n' % (ns[0], ns[-1]))
if __name__ == '__main__': main()
