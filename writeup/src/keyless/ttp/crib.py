"""crib.py — MX（= M_A ⊕ M_E）の並びに、E 側の文字列を当てて A 側を見る（手作業の補助）
  python3 -I crib.py FILE.pkl 位置 '文字列' [位置 '文字列' …]"""
import sys, pickle, ast
d = pickle.load(open(sys.argv[1], 'rb')); base = 4 * d['lo']
mx = b''.join(d['mxs'][n] for n in range(d['lo'], d['hi'] + 1))
def sj(t):
    out = []; i = 0
    while i < len(t):
        c = t[i]
        if c == 0x0d: out.append('↵'); i += 1
        elif c == 0x0a: out.append('¶'); i += 1
        elif c < 0x20 or c == 0x7f: out.append('<%02x>' % c); i += 1
        elif c < 0x80: out.append(chr(c)); i += 1
        elif (0x81 <= c <= 0x9f or 0xe0 <= c <= 0xfc) and i + 1 < len(t):
            try: out.append(t[i:i + 2].decode('cp932')); i += 2
            except Exception: out.append('<%02x>' % c); i += 1
        elif 0xa1 <= c <= 0xdf: out.append(bytes([c]).decode('cp932')); i += 1
        else: out.append('<%02x>' % c); i += 1
    return ''.join(out)
a = sys.argv[2:]
for k in range(0, len(a), 2):
    p = int(a[k]); s = a[k + 1].encode('cp932').decode('unicode_escape').encode('latin1') if '\\' in a[k + 1] else a[k + 1].encode('cp932')
    o = bytes(x ^ y for x, y in zip(s, mx[p - base:p - base + len(s)]))
    print('%4d〜%4d  %r' % (p, p + len(s) - 1, s)); print('           → %s   [%s]' % (sj(o), o.hex()))
