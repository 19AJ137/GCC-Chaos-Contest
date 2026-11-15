# sigma windows of the two contest keys and what they say about the per-position character sums (ASCII hypotheses)
import sys
from sums import *

def windows(name, ns):
    t = TG[name]; out = {}
    for n in ns:
        w = sigma_window(t['S'], t['Tlo'], t['Thi'], n)
        if not w: continue
        sig = sorted(set(x[1] for x in w)); wr = sorted(set(x[0] for x in w))
        sj = []
        for (ww, s, I3, T16) in w:
            for j1 in range(n):
                v = t['S'] + (ww << 32) - 3 * n - (fl_mul(MP, SHP, s) - j1)
                if v >= 4: sj.append(v)
        out[n] = (sig, wr, min(sj), max(sj))
    return out

def decomp(sig, nfull, lo=0x20, hi=0x7E, part=0):
    """(B0,B1,B2,B3) with sigma = sum 256^j B_j; B_j = sum of byte j over nfull full chunks (+ partial chunk bytes for j < part)"""
    res = []
    cnt = [nfull + (1 if j < part else 0) for j in range(4)]
    def rec(j, rem, acc):
        if j == 4:
            if rem == 0: res.append(tuple(acc))
            return
        if j == 3:
            if rem % 1 == 0 and cnt[3] * lo <= rem <= cnt[3] * hi: res.append(tuple(acc + [rem]))
            return
        b = rem & 0xFF
        v = b
        while v <= cnt[j] * hi:
            if v >= cnt[j] * lo: rec(j + 1, (rem - v) >> 8, acc + [v])
            v += 256
    rec(0, sig, [])
    return res

if __name__ == '__main__':
    for name, ns in (('k', range(3, 8)), ('v', range(8, 13))):
        W = windows(name, ns)
        for n, (sig, wr, sjlo, sjhi) in W.items():
            print('%s n=%d: sigma %d..%d (%d values), wrap %s, sum of J %d..%d, sigma/n = 0x%08x' % (name, n, sig[0], sig[-1], len(sig), wr, sjlo, sjhi, sig[len(sig) // 2] // n))
    print()
    for name, n, part in (('k', 5, 0), ('k', 5, 1), ('k', 5, 2), ('k', 5, 3), ('k', 6, 0), ('v', 9, 0), ('v', 10, 0)):
        W = windows(name, [n + (1 if part else 0)])
        if not W: continue
        sig = list(W.values())[0][0]
        b3 = {}; b2 = set(); b1 = set(); tot = 0
        for s in sig:
            for (B0, B1, B2, B3) in decomp(s, n, part=part):
                b3[B3] = b3.get(B3, 0) + 1; b2.add(B2); b1.add(B1); tot += 1
        if not tot: print('%s %d full chunks + %d-byte tail: no printable-ASCII decomposition' % (name, n, part)); continue
        print('%s %d full chunks + %d-byte tail (printable ASCII): B3 (sum of chars 4,8,..) in %s (mean char %.1f); B2 %d values %d..%d; B1 %d values' % (
            name, n, part, sorted(b3), sum(k * v for k, v in b3.items()) / tot / n, len(b2), min(b2), max(b2), len(b1)))
