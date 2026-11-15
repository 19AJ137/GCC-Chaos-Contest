# 見つかった鍵 k の中身（かたまりごとの値・格子・整数部・和）を並べる。答え合わせ用。
import math, sys
from gks import key_schedule, chunk_vals, chunks_of, map_of_i2, K_PI, K_E, h2d
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from possum import windows
key = bytes.fromhex(sys.argv[1]) if len(sys.argv) > 1 else '磁石・発光・IISI'.encode('cp932')
r = key_schedule(key)
print('key %s  (%d bytes)  %r' % (key.hex(), len(key), key.decode('cp932')))
print('map %d  S %08x  T16 %04x  warm %d  J %s' % (r['m'], r['S'], r['b0b1'], r['warm'], r['J']))
print('k780 %.16f  k788 %.16f  k7b0 %.16f  k7b8 %.16f' % tuple(r['k']))
cs = chunks_of(key); sig = 0; I1 = I3 = 0; P = [1.0] * 4
def lat(c, K):
    e = math.floor(math.log2(c * K)); l = 44 - e
    return l
for c in cs:
    (f, i, s, t) = chunk_vals(c)
    l1 = lat(c, K_PI); l3 = lat(c, K_E)
    N1 = (1 + f[0]) * 2**l1; N3 = (1 + f[2]) * 2**l3
    b = c.to_bytes(4, 'little')
    print(' c %08x %-6r i1 %10d f1 %.12f l1 %d N1 %-8s | i3 %10d f3 %.12f l3 %d N3 %-8s | i2 %10d i4 %10d f2 %.6f f4 %.6f | s %d t %d map-if-last %d' % (
        c, b.decode('cp932', errors='replace'), i[0], f[0], l1, ('%d' % N1) if N1 == int(N1) else repr(N1), i[2], f[2], l3, ('%d' % N3) if N3 == int(N3) else repr(N3),
        i[1], i[3], f[1], f[3], s, t, map_of_i2(i[1])))
    sig += c; I1 += i[0]; I3 += i[2]
    for j in range(4): P[j] *= 1 + f[j]
n = len(cs)
print('sigma = sum c = %d  (sigma/2^24 = %.4f)   sum i1 = %d  sum i3 = %d   n = %d' % (sig, sig / 2**24, I1, I3, n))
print('products: %s' % ' '.join('%.12f' % p for p in P))
print('S check: sum i1 + sum ftol(i3/i2) + sum J = ?  S - sum i1 - sum J = %d (ratio terms; 3n = %d)' % ((r['S'] - I1 - sum(r['J'])) % 2**32, 3 * n))
bs = [sum((c >> (8 * j)) & 0xff for c in cs) for j in range(4)]
print('byte-position sums B0..B3 = %s' % bs)
W = windows('k', [n])
if n in W:
    sg = W[n][0]; print('sigma window for n=%d: %d..%d (%d values); true sigma is #%d' % (n, sg[0], sg[-1], len(sg), sg.index(sig) if sig in sg else -1))
d = (r['S'] * (math.pi / math.e) - r['b0b1']) % 65536
print('d = S*(pi/e) - T16 mod 65536 = %.1f' % d)
