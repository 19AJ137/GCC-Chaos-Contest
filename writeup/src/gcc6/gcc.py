"""GCC (Gao's Chaos Cryptosystem) — Python re-implementation of Gcc6_32.dll (1997-01-23 build).

Reconstructed from the disassembly of Str_cryption (0x10001950), File_cryption (0x10001010),
the key schedule (0x100020e0) and map setup (0x10002600).
The DLL's CRT init sets the x87 precision control to 53 bits (_controlfp(_PC_53, _MCW_PC)),
so every operation below is plain IEEE-754 double arithmetic in the same order as the binary.
"""
import math
import struct

def _d(h):
    return struct.unpack('<d', bytes.fromhex(h))[0]

# Constants exactly as stored in .rdata (note: "1/pi" and the second "0.1" are not the exact values)
C_INV_PI = _d('a535da6d305fd43f')   # 0.3183098862435492   (true 1/pi = 0.3183098861837907)
C_INV_E  = _d('81cf1836568bd73f')   # 0.36787944109823206
C_32767  = 32767.0
C_TENTH  = _d('9a9999999999b93f')   # 0.1
C_X0_DEF = _d('e4cfe0847c8fe93f')   # 0.79876543
C_TENTH2 = _d('a09999999999b93f')   # 0.10000000000000009
C_3_9    = _d('3333333333330f40')
C_1_9    = _d('666666666666fe3f')
C_0_4    = _d('9a9999999999d93f')
C_0_9    = _d('cdccccccccccec3f')
C_M014   = _d('ec51b81e85ebc1bf')
C_0_8    = _d('9a9999999999e93f')
C_1_72   = _d('85eb51b81e85fb3f')


def d2b(x):
    return struct.pack('<d', x)


def b2d(b):
    return struct.unpack('<d', bytes(b))[0]


def clear_low_byte(x):
    b = bytearray(d2b(x)); b[0] = 0
    return b2d(b)


def set_low16(x, b0, b1):
    b = bytearray(d2b(x)); b[0] = b0; b[1] = b1
    return b2d(b)


def ftol(x):
    """MSVC _ftol: fistp qword with truncation, return low dword (eax)."""
    if not math.isfinite(x) or abs(x) >= 2.0**63:
        return 0  # integer indefinite 0x8000000000000000 -> eax = 0
    return int(x) & 0xFFFFFFFF


def fdiv(a, b):
    if b == 0.0:
        if a == 0.0 or math.isnan(a):
            return math.nan
        return math.copysign(math.inf, a) * math.copysign(1.0, b)
    return a / b


class GCCState:
    """Internal state after the key schedule."""
    __slots__ = ('mode', 'x', 'r', 'y', 'u', 'v', 'Q', 'A', 'C', 'b0', 'b1',
                 'S', 'T', 'warm', 'k780', 'k788', 'k7b0', 'k7b8', 'k7b0_saved')

    def copy(self):
        n = GCCState()
        for s in self.__slots__:
            setattr(n, s, getattr(self, s, None))
        return n


def setup_map(st):
    """0x10002600: initialise map variables from k780/k788/k7b0/k7b8."""
    if st.mode == 0:
        x = st.k780
        if x == 0.5 or not (0.0 < x):
            x = C_X0_DEF
            st.k780 = x
        st.x = x
        st.r = st.k7b0 * C_TENTH2 + C_3_9
    elif st.mode == 1:
        st.y = st.k780 * 2.0 - 1.0
        st.r = st.k7b0 * C_TENTH2 + C_1_9          # p
    elif st.mode == 2:
        st.u = st.k780
        st.v = st.k788
        st.r = st.k7b0 * C_0_4 + C_0_9             # P  (c760)
        Q = (st.k7b8 * C_M014 - st.r * C_0_8) + C_1_72
        st.Q = Q                                    # c768
        st.A = Q * st.r + 1.0                       # c770
        st.C = 1.0 - Q                              # c778
    # low 16 bits of c760/c768/c770/c778 are overwritten with b0,b1 (all modes)
    st.r = set_low16(st.r, st.b0, st.b1)
    if st.mode == 1:
        st.y = set_low16(st.y, st.b0, st.b1)
    if st.mode == 2:
        st.Q = set_low16(st.Q, st.b0, st.b1)
        st.A = set_low16(st.A, st.b0, st.b1)
        st.C = set_low16(st.C, st.b0, st.b1)
    # (for modes 0/1 the unused c768/c770/c778 slots are also touched, irrelevant)


def step(st):
    """One map iteration; returns the output double (c750)."""
    b0, b1 = st.b0, st.b1
    if st.mode == 0:
        st.x = set_low16(((1.0 - st.x) * st.x) * st.r, b0, b1)
        return st.x
    if st.mode == 1:
        t = set_low16(st.y * st.y - st.r, b0, b1)
        st.y = t
        st.x = (1.0 + t) * 0.5
        return st.x
    if st.mode == 2:
        a = st.Q * st.u
        b = a * st.u
        nu = st.A * st.u - a * st.v
        nv = st.C * st.v + b
        st.u = set_low16(nu, b0, b1)
        st.v = set_low16(nv, b0, b1)
        st.x = st.v
        return st.x
    raise ValueError(st.mode)


def key_schedule(key: bytes) -> GCCState:
    st = GCCState()
    k780 = k788 = k7b0 = k7b8 = 1.0
    S = 0
    T = 0
    n = len(key) & 0xFFFF
    if n == 0:
        raise ValueError('empty key: DLL returns without initialising the map')
    off = 0
    rem = n
    ip2 = 0.0
    while rem > 0:
        m = min(rem, 4)
        chunk = int.from_bytes(key[off:off + m] + b'\0' * (4 - m), 'little')
        dc = float(chunk)
        f1, ip1 = math.modf(clear_low_byte(dc * C_INV_PI))
        f2, ip2 = math.modf(clear_low_byte(ip1 * C_INV_PI))
        k780 = (1.0 + f1) * k780
        k788 = (f2 + 1.0) * k788
        f3, ip3 = math.modf(clear_low_byte(dc * C_INV_E))
        f4, ip4 = math.modf(clear_low_byte(ip3 * C_INV_E))
        k7b0 = (1.0 + f3) * k7b0
        k7b8 = (f4 + 1.0) * k7b8
        S = (S + ((ftol(fdiv(ip3, ip2)) + ftol(ip1)) & 0xFFFFFFFF)) & 0xFFFFFFFF
        T = (T + ((ftol(fdiv(ip1, ip4)) + ftol(ip3)) & 0xFFFFFFFF)) & 0xFFFFFFFF
        off += m
        rem -= 4
    vals = []
    for k in (k780, k788, k7b0, k7b8):
        if 1.0 < k:
            f, ip = math.modf(k)
            S = (S + ftol(ip)) & 0xFFFFFFFF
            k = f
        vals.append(k)
    k780, k788, k7b0, k7b8 = vals
    warm = S % 79
    if warm < 11:
        warm = 11
    st.b0 = T & 0xFF
    st.b1 = (T >> 8) & 0xFF
    v = ip2
    if not (v <= C_32767):
        while not (v <= C_32767):
            v = C_TENTH * v
    ax = ftol(v) & 0xFFFF
    if ax >= 0x8000:
        ax -= 0x10000
    st.mode = int(math.fmod(ax, 3))  # idiv: remainder has the sign of the dividend
    st.S, st.T, st.warm = S, T, warm
    st.k780, st.k788, st.k7b0, st.k7b8 = k780, k788, k7b0, k7b8
    st.k7b0_saved = k7b0
    if st.mode not in (0, 1, 2):
        raise ValueError('negative mode %d: map stays uninitialised' % st.mode)
    setup_map(st)
    for _ in range(warm):
        step(st)
    return st


def ks_word(x):
    """Keystream dword and 'no-chain' flag from the output double."""
    b = d2b(x)
    return bytes((b[4], b[2], b[3], b[5])), b[6] > 0xE6


def _xor(a, b):
    return bytes(i ^ j for i, j in zip(a, b))


def crypt_block(st, data: bytes, encrypt: bool, chain: bytes, bl: int):
    """Process one buffer exactly like the DLL inner loop. Returns (out, chain, bl)."""
    n = len(data)
    tail = n % 4
    words = (n - tail) // 4
    out = bytearray()
    for i in range(words):
        ks, nochain = ks_word(step(st))
        blk = data[4 * i:4 * i + 4]
        if nochain:
            out += _xor(blk, ks)
        else:
            o = _xor(_xor(blk, ks), chain)
            out += o
            chain = o if encrypt else blk
    if tail:
        ks, _ = ks_word(step(st))
        for j in range(tail):
            c = ks[j] ^ data[n - tail + j] ^ bl
            out.append(c)
            bl = c if encrypt else data[n - tail + j]
    return bytes(out), chain, bl


def str_cryption(key: bytes, data: bytes, encrypt=True):
    """Str_cryption(mode, in, len, out, key). len is a 16-bit count."""
    st = key_schedule(key)
    iv = st.S.to_bytes(4, 'little')
    out, _, _ = crypt_block(st, data, encrypt, iv, iv[1])
    return out


def file_cryption(key: bytes, data: bytes, encrypt=True):
    """File_cryption: 1024-byte blocks; after each block the map is re-seeded (no warm-up) and
    byte 5 of k7b0 is replaced by byte 1 of the current chaining dword."""
    st = key_schedule(key)
    chain = st.S.to_bytes(4, 'little')
    bl = chain[1]
    out = bytearray()
    for off in range(0, len(data), 1024):
        blk = data[off:off + 1024]
        o, chain, bl = crypt_block(st, blk, encrypt, chain, bl)
        out += o
        b = bytearray(d2b(st.k7b0_saved)); b[5] = chain[1]
        st.k7b0 = b2d(b)
        setup_map(st)
    return bytes(out)


if __name__ == '__main__':
    import sys
    k = sys.argv[1].encode() if len(sys.argv) > 1 else b'password'
    st = key_schedule(k)
    print('mode', st.mode, 'warm', st.warm, 'S %08x T %08x' % (st.S, st.T),
          'k780', st.k780, 'k7b0', st.k7b0, 'r', st.r)
    ct = str_cryption(k, b'Hello, chaos world! 0123456789')
    print(ct.hex())
    print(str_cryption(k, ct, encrypt=False))
