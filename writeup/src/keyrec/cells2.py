# 2-D feasible regions: (k7b0, k7b8) and (k780, k788) of key k -- bounding boxes (what an attack must allow for each constant)
import sys
from cells import *

def region(n1, n2, r1, r2, N=61):
    pts = []
    for i in range(N):
        d1 = (2 * i / (N - 1) - 1) * r1
        for j in range(N):
            d2 = (2 * j / (N - 1) - 1) * r2
            d = dict(K); d[n1] = K[n1] + d1; d[n2] = K[n2] + d2
            if n1 == 'k7b0' and struct.pack('<d', d[n1])[5] != struct.pack('<d', K[n1])[5]: continue
            if any(full_k(d, b)[:12] + full_k(d, b)[256:268] == refk[:12] + refk[256:268] for b in (0x7ed6,)) and any(full_k(d, b) == refk for b in BK):
                pts.append((d1, d2))
    return pts

for (n1, n2, r1, r2) in (('k7b0', 'k7b8', 4.5e-11, 1.0e-10), ('k780', 'k788', 1.2e-11, 1.6e-11)):
    pts = region(n1, n2, r1, r2)
    a = [p[0] for p in pts]; b = [p[1] for p in pts]
    print('%s x %s: %d feasible grid points of 3721; %s offsets [%.3e, %.3e] (width %.2e), %s offsets [%.3e, %.3e] (width %.2e)' % (
        n1, n2, len(pts), n1, min(a), max(a), max(a) - min(a), n2, min(b), max(b), max(b) - min(b)))
    print('   grid step %.2e / %.2e; touches scan border: %s' % (2 * r1 / 60, 2 * r2 / 60, min(a) <= -r1 * 0.999 or max(a) >= r1 * 0.999 or min(b) <= -r2 * 0.999 or max(b) >= r2 * 0.999))
