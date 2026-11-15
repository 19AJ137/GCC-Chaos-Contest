# How precisely does the recovered keystream pin each key-schedule output?  (interval of each constant that reproduces
# block 1 (256 words after the warm-up) and block 2 of E (79 words, k7b0 byte 5 replaced, no warm-up) for some b0b1 in the stated range)
import struct, math
from gks import *
from maps import setlow, ksword

def ks_map2(k780, k788, k7b0, k7b8, b, warm, n):
    u, v = k780, k788
    P = 0.4 * k7b0 + 0.9; Q = (-0.14 * k7b8 - 0.8 * P) + 1.72; A = Q * P + 1.0; C = 1.0 - Q
    Q = setlow(Q, b); A = setlow(A, b); C = setlow(C, b)
    out = []
    for i in range(warm + n):
        a = Q * u; bb = a * u; un = A * u - a * v; vn = C * v + bb
        u = setlow(un, b); v = setlow(vn, b)
        if i >= warm: out.append(ksword(v))
    return out

def ks_map1(k780, k7b0, b, warm, n):
    y = setlow(2.0 * k780 - 1.0, b); p = setlow(k7b0 * C01b + 1.9, b); out = []
    for i in range(warm + n):
        y = setlow(y * y - p, b)
        if i >= warm: out.append(ksword((1.0 + y) / 2.0))
    return out

def rep5(x, byte):
    q = bytearray(struct.pack('<d', x)); q[5] = byte; return struct.unpack('<d', bytes(q))[0]

K = dict(k780=h2d('3fe500287933f358'), k788=h2d('3fe4f462e3890aef'), k7b0=h2d('3fc3ca0d1ca9ebf4'), k7b8=h2d('3fe924aef33bd88d'))
V = dict(k780=h2d('3fe99b19be354000'), k7b0=h2d('3fcf8c8bdb300652'))
BK = [0x7ec1, 0x7ed0, 0x7ed6, 0x7ee0, 0x7eec]; BV = [0xab0f, 0xab40, 0xab6a, 0xaba0, 0xabc5]

def full_k(d, b): return ks_map2(d['k780'], d['k788'], d['k7b0'], d['k7b8'], b, 69, 256) + ks_map2(d['k780'], d['k788'], rep5(d['k7b0'], 0xae), d['k7b8'], b, 0, 79)
def full_v(d, b): return ks_map1(d['k780'], d['k7b0'], b, 23, 256) + ks_map1(d['k780'], rep5(d['k7b0'], 0xed), b, 0, 79)

refk = full_k(K, 0x7ed6); refv = full_v(V, 0xab6a)
print('k: all listed b0b1 give the same 335 words:', all(full_k(K, b) == refk for b in BK), ' block-2 word0 %08x' % refk[256][0])
print('v: all listed b0b1 give the same 335 words:', all(full_v(V, b) == refv for b in BV), ' block-2 word0 %08x' % refv[256][0])
print('k flags block1 (.=no chain):', ''.join('.' if f else '1' for (_, f) in refk[:64]))
print('v flags block1             :', ''.join('.' if f else '1' for (_, f) in refv[:64]))

def ok(full, base, name, val, ref, bs):
    d = dict(base); d[name] = val
    if name == 'k7b0' and struct.pack('<d', val)[5] != struct.pack('<d', base['k7b0'])[5]: return False
    return any(full(d, b) == ref for b in bs)

def edge(full, base, name, ref, bs, sign):
    lo = 0.0; hi = 1e-9                      # offset magnitude: ok at lo, not ok at hi
    if ok(full, base, name, base[name] + sign * hi, ref, bs): return None
    for _ in range(70):
        mid = (lo + hi) / 2
        if ok(full, base, name, base[name] + sign * mid, ref, bs): lo = mid
        else: hi = mid
    return lo

for label, full, base, ref, bs in (('k', full_k, K, refk, BK), ('v', full_v, V, refv, BV)):
    for name in base:
        a = edge(full, base, name, ref, bs, -1); b = edge(full, base, name, ref, bs, +1)
        print('%s %s = %.17g  cell [-%.3e, +%.3e]  width %.3e  = 2^%.1f' % (label, name, base[name], a, b, a + b, math.log2(a + b)))
