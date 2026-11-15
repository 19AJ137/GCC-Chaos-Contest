# GCC6 key schedule + maps, re-implemented from the spec (docs/algorithm.md) for the key-recovery tools.
# Not the reference ../gcc6/gcc.py; xcheck_gcc.py compares the two (6,000 random keys, no mismatch).
# Also checked against known values (slide key gao1234 and 6 control keys).
import struct, math

def d2h(x): return struct.pack('>d', x).hex()
def h2d(h): return struct.unpack('>d', bytes.fromhex(h))[0]

K_PI = struct.unpack('<d', bytes.fromhex('a535da6d305fd43f'))[0]
K_E  = struct.unpack('<d', bytes.fromhex('81cf1836568bd73f'))[0]
C01  = struct.unpack('<d', bytes.fromhex('9a9999999999b93f'))[0]   # 0.1
C01b = struct.unpack('<d', bytes.fromhex('a09999999999b93f'))[0]   # 0.10000000000000009
M32 = 0xFFFFFFFF

def clr(x):
    b = bytearray(struct.pack('<d', x)); b[0] = 0
    return struct.unpack('<d', bytes(b))[0]

def ftol(x):
    if x != x or x in (math.inf, -math.inf) or abs(x) >= 2.0**63: return 0
    return int(x) & M32

def fdiv(a, b):
    if b == 0:
        if a == 0: return math.nan
        return math.inf if a > 0 else -math.inf
    return a / b

def chunk_vals(c):
    f1, i1 = math.modf(clr(c * K_PI))
    f2, i2 = math.modf(clr(i1 * K_PI))
    f3, i3 = math.modf(clr(c * K_E))
    f4, i4 = math.modf(clr(i3 * K_E))
    s = (ftol(fdiv(i3, i2)) + ftol(i1)) & M32
    t = (ftol(fdiv(i1, i4)) + ftol(i3)) & M32
    return (f1, f2, f3, f4), (i1, i2, i3, i4), s, t

def chunks_of(key):
    return [int.from_bytes(key[i:i+4].ljust(4, b'\0'), 'little') for i in range(0, len(key), 4)]

def map_of_i2(i2):
    x = i2
    while x > 32767.0: x *= C01
    v = ftol(x) & 0xFFFF
    if v >= 0x8000: v -= 0x10000
    return int(math.fmod(v, 3))

def key_schedule(key):
    k = [1.0] * 4; S = 0; T = 0; i2 = 0.0
    for c in chunks_of(key):
        f, i, s, t = chunk_vals(c)
        for j in range(4): k[j] *= (1.0 + f[j])
        S = (S + s) & M32; T = (T + t) & M32
        i2 = i[1]
    J = [0] * 4
    for j in range(4):
        if k[j] > 1.0:
            fr, ip = math.modf(k[j]); k[j] = fr; J[j] = int(ip); S = (S + int(ip)) & M32
    return dict(k=k, S=S, T=T, J=J, warm=max(S % 79, 11), b0b1=T & 0xFFFF, m=map_of_i2(i2))

if __name__ == '__main__':
    print('K_PI', repr(K_PI), 'K_E', repr(K_E))
    tests = [(b'gao1234', 0x0FCCF2DC, 0xdb8f), (bytes.fromhex('37265e45'), 0x16149836, 0xdba9),
             (bytes.fromhex('34685a3146'), 0x0fb5a83d, 0xf065), (bytes.fromhex('484d715d365b31'), 0x1dce1902, 0x4ff7),
             (bytes.fromhex('301b4001d7de'), 0x00662b94, 0x14ad), (bytes.fromhex('dadf560a098ab0'), 0x0382b87f, 0xad05),
             (bytes.fromhex('abbfbb064a495f04'), 0x0388f762, 0xe4e8), (bytes.fromhex('4a495f04abbfbb06'), 0x0388f762, 0xe4e8)]
    for key, S, T16 in tests:
        r = key_schedule(key)
        print(key.hex(), 'S %08x %s' % (r['S'], 'ok' if r['S'] == S else 'NG'), 'T16 %04x %s' % (r['b0b1'], 'ok' if r['b0b1'] == T16 else 'NG'),
              'map', r['m'], 'warm', r['warm'], 'k7b0 %.17g' % r['k'][2])
    r = key_schedule(b'gao1234')
    print('gao1234 k7b0 in slide interval:', 0.49611561009068333 <= r['k'][2] <= 0.49611561023640022)
    a = key_schedule(b'testkey1'); b = key_schedule(b'key1test')
    print('testkey1 vs key1test:', a['k'] == b['k'], a['S'] == b['S'], a['T'] == b['T'], a['m'], b['m'])
