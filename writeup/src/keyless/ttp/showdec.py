"""showdec.py — fulldec2.py の出力を 16 バイトずつ、MX・E・A を並べて表示する
  python3 -I showdec.py FILE.pkl [バイト lo] [バイト hi]"""
import sys, pickle, os
d = pickle.load(open(sys.argv[1], 'rb')); lo = int(sys.argv[2]) if len(sys.argv) > 2 else 4 * d['lo']; hi = int(sys.argv[3]) if len(sys.argv) > 3 else 4 * d['hi'] + 3
base = 4 * d['lo']; e = d['e']; a = d['a']
mx = b''.join(d['mxs'][n] for n in range(d['lo'], d['hi'] + 1))
def ch(c): return chr(c) if 0x20 <= c < 0x7f else ('¶' if c == 0x0a else '↵' if c == 0x0d else '·')
def sj(t):
    try: return t.decode('cp932')
    except Exception: return t.decode('cp932', 'replace')
for p in range(lo - lo % 16, hi + 1, 16):
    i = p - base; seg_e = e[i:i + 16]; seg_a = a[i:i + 16]; seg_m = mx[i:i + 16]
    fl = ' '.join('%d%d' % d['flags'][n] for n in range(p // 4, min(p // 4 + 4, d['hi'] + 1)))
    print('%4d  %s  [%s]' % (p, ' '.join(seg_m[k:k + 4].hex() for k in range(0, len(seg_m), 4)), fl))
    print('      E: %-18s A: %-18s %s' % (''.join(ch(c) for c in seg_e), ''.join(ch(c) if c < 0x80 else '#' for c in seg_a), seg_a.hex()))
