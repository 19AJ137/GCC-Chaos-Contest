# Maps and keystream words from a key-schedule output, to cross-check the targets (../TARGETS.json).
import struct, math, sys
from gks import *

def setlow(x, b0b1):
    q = struct.unpack('<Q', struct.pack('<d', x))[0]
    q = (q & ~0xFFFF) | b0b1
    return struct.unpack('<d', struct.pack('<Q', q))[0]

def ksword(x):
    z = struct.pack('<d', x)
    return (z[4] << 24) | (z[2] << 16) | (z[3] << 8) | z[5], (z[6] > 0xE6)

def run_map1(k780, k7b0, b0b1, warm, n, variant=0):
    y = setlow(2.0 * k780 - 1.0, b0b1)
    p = setlow(k7b0 * C01b + 1.9, b0b1)
    out = []
    for i in range(warm + n):
        y = setlow(y * y - p, b0b1)
        if i >= warm:
            out.append(ksword((1.0 + y) / 2.0))
    return out, p

def run_map2(k780, k788, k7b0, k7b8, b0b1, warm, n, variant=0):
    u, v = k780, k788
    P = 0.4 * k7b0 + 0.9
    Q = (-0.14 * k7b8 - 0.8 * P) + 1.72
    if variant == 0:      # A, C from Q before overwrite; then overwrite all
        A = Q * P + 1.0; C = 1.0 - Q
        P = setlow(P, b0b1); Q = setlow(Q, b0b1); A = setlow(A, b0b1); C = setlow(C, b0b1)
    elif variant == 1:    # overwrite P, Q first
        P = setlow(P, b0b1); Q = setlow(Q, b0b1)
        A = setlow(Q * P + 1.0, b0b1); C = setlow(1.0 - Q, b0b1)
    elif variant == 2:    # overwrite P first only
        P = setlow(P, b0b1)
        Q = (-0.14 * k7b8 - 0.8 * P) + 1.72
        A = Q * P + 1.0; C = 1.0 - Q
        Q = setlow(Q, b0b1); A = setlow(A, b0b1); C = setlow(C, b0b1)
    elif variant == 3:
        P = setlow(P, b0b1)
        Q = setlow((-0.14 * k7b8 - 0.8 * P) + 1.72, b0b1)
        A = setlow(Q * P + 1.0, b0b1); C = setlow(1.0 - Q, b0b1)
    out = []
    for i in range(warm + n):
        a = Q * u; b = a * u
        un = A * u - a * v
        vn = C * v + b
        u = setlow(un, b0b1); v = setlow(vn, b0b1)
        if i >= warm:
            out.append(ksword(v))
    return out, (P, Q, A, C)

if __name__ == '__main__':
    tv = '8c8c921c 66070390 96297436 d52c64da 518ae3f8 77476348 e30bd159 275c8044 e877a09c'.split()
    tk = 'cfa907c7 d4425a3b 46c22e5f 3386c848 02cbcf9e d6259668 a41a2c6b 65c96c1c a23788dc'.split()
    k780v = h2d('3fe99b19be354000'); k7b0v = h2d('3fcf8c8bdb300652')
    for kk in (k780v, 1.0 - k780v):
        okb = []
        for b in range(0xab00, 0xac00):
            w, p = run_map1(kk, k7b0v, b, 23, 9)
            if ['%08x' % x[0] for x in w][1:] == tv[1:]: okb.append((b, '%08x' % w[0][0]))
        print('v: k780=%.16g  b0b1 matching words1-8: %d' % (kk, len(okb)), ('%04x..%04x' % (okb[0][0], okb[-1][0]), set(x[1] for x in okb)) if okb else '')
    w, p = run_map1(k780v, k7b0v, 0xab50, 23, 9); print('v p', d2h(p), repr(p), ['%08x%s' % (x[0], '.' if x[1] else '1') for x in w])
    k = [h2d(x) for x in ('3fe500287933f358', '3fe4f462e3890aef', '3fc3ca0d1ca9ebf4', '3fe924aef33bd88d')]
    for var in range(4):
        okb = []
        for b in range(0x7e00, 0x7f00):
            w, par = run_map2(k[0], k[1], k[2], k[3], b, 69, 9, var)
            if ['%08x' % x[0] for x in w][1:] == tk[1:]: okb.append((b, '%08x' % w[0][0]))
        print('k variant', var, 'b0b1 matching words1-8:', len(okb), ('%04x..%04x' % (okb[0][0], okb[-1][0]), set(x[1] for x in okb)) if okb else '')
        w, par = run_map2(k[0], k[1], k[2], k[3], 0x7ed0, 69, 9, var)
        print('   P,Q,A,C', [d2h(x) for x in par], ['%08x%s' % (x[0], '.' if x[1] else '1') for x in w][:4])
