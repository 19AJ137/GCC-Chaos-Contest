# What S (=IV, 32 bits) and the b0b1 range say about sigma = sum of chunk values, I1, I3, F1, F3 and the integer parts J.
# Exact integer arithmetic (K_PI, K_E are dyadic rationals).
import math, struct, sys
from gks import K_PI, K_E, h2d

def dyadic(x):
    m, e = math.frexp(x); m = int(m * (1 << 53)); e -= 53
    return m, -e                      # x = m / 2^sh
MP, SHP = dyadic(K_PI); ME, SHE = dyadic(K_E)
def fl_mul(M, sh, v): return (M * v) >> sh              # floor(K*v)
def fr_mul(M, sh, v): return ((M * v) & ((1 << sh) - 1)) / float(1 << sh)   # frac(K*v)

TG = {
 'k': dict(S=0xa678044d, Tlo=0x7ec1, Thi=0x7eec, k780=[h2d('3fe500287933f358')], k788=h2d('3fe4f462e3890aef'),
           k7b0=h2d('3fc3ca0d1ca9ebf4'), k7b8=h2d('3fe924aef33bd88d'), m=2),
 'v': dict(S=0x3a9d3481, Tlo=0xab0f, Thi=0xabc5, k780=[h2d('3fe99b19be354000'), 1.0 - h2d('3fe99b19be354000')], k788=None,
           k7b0=h2d('3fcf8c8bdb300652'), k7b8=None, m=1),
}

def prod_bounds(F, n):
    j = math.floor(F)
    return 2.0**j * (1.0 + (F - j)), (1.0 + F / n)**n

def jrange(F, n, k, slack=2e-3):
    lo, hi = prod_bounds(F, n)
    return [J for J in range(max(1, int(lo - k - 1)), int(hi - k + 2)) if lo - slack <= J + k <= hi + slack]

def sigma_window(S, Tlo, Thi, n, r1=None, r3=None, sjmax=None):
    """all (sigma, I3, T16) consistent with S and the T16 range for n chunks (ratio terms r1=3n, r3=2n by default)"""
    r1 = 3 * n if r1 is None else r1; r3 = 2 * n if r3 is None else r3
    sjmax = 4 * (2**n - 1) if sjmax is None else sjmax
    out = []
    for w in range(0, n + 1):
        Sw = S + (w << 32)
        lo = int((Sw - r1 - sjmax) / K_PI) - 3; hi = int((Sw - r1 - 4 + n) / K_PI) + 3
        if hi < 0 or lo > n * 0xFFFFFFFF: continue
        lo = max(lo, 0)
        # T16 = (I3 + r3) mod 2^16 with I3 = floor(K_E sigma) - j3: step through the window
        for sig in range(lo, hi + 1):
            f3 = fl_mul(ME, SHE, sig)
            for j3 in range(n):
                I3 = f3 - j3; T16 = (I3 + r3) & 0xFFFF
                if Tlo <= T16 <= Thi:
                    out.append((w, sig, I3, T16))
    return out

def analyse(t, n):
    S = t['S']; res = []
    for (w, sig, I3, T16) in sigma_window(S, t['Tlo'], t['Thi'], n):
        Sw = S + (w << 32)
        F3 = (fl_mul(ME, SHE, sig) - I3) + fr_mul(ME, SHE, sig)
        J3s = jrange(F3, n, t['k7b0'])
        if not J3s: continue
        f1 = fl_mul(MP, SHP, sig)
        for j1 in range(n):
            I1 = f1 - j1; F1 = j1 + fr_mul(MP, SHP, sig)
            SJ = Sw - 3 * n - I1
            if SJ < 4: continue
            for bi, k1 in enumerate(t['k780']):
                J1s = jrange(F1, n, k1)
                for a in J1s:
                    for b in J3s:
                        rest = SJ - a - b
                        if rest < 2 or rest > 2 * (2**n - 1): continue
                        res.append(dict(n=n, w=w, sig=sig, T16=T16, I1=I1, I3=I3, F1=F1, F3=F3, J1=a, J3=b, J24=rest, SJ=SJ, br=bi))
    return res

if __name__ == '__main__':
    for name, (a, b) in (('k', (3, 6)),):          # the v key has too many tuples for this listing: use possum.py
        t = TG[name]
        print('==', name, 'S=%d' % t['S'])
        for n in range(a, b + 1):
            r = analyse(t, n)
            if not r: print(' n=%d: none' % n); continue
            for w in sorted(set(x['w'] for x in r)):
                L = [x for x in r if x['w'] == w]
                sigs = sorted(set(x['sig'] for x in L)); t16 = sorted(set(x['T16'] for x in L))
                print(' n=%d wrap=%d: tuples %d, sigma %d values [%d..%d], T16 %d values [%04x..%04x], SJ %d..%d, J1 %s J3 %s' % (
                    n, w, len(L), len(sigs), sigs[0], sigs[-1], len(t16), t16[0], t16[-1], min(x['SJ'] for x in L), max(x['SJ'] for x in L),
                    sorted(set(x['J1'] for x in L)), sorted(set(x['J3'] for x in L))))
