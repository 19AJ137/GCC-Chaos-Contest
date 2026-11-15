import sys, os, math, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from vpart import textlike
from vshort import load_hi, SIGMA, A_U, B_U
from vshort_plan import tau
from gks import chunk_vals
b4 = lambda c: c.to_bytes(4, 'little'); D = B_U - A_U
hi = load_hi(); rng = random.Random(9)
typ = [x for x in hi if min(x[1]) >= 0.95 and textlike(b4(x[0])) >= 2 and x[0] >= 0xA8000000]
print('pool', len(typ))
fl = {b: chunk_vals(b)[0] for b in (0x51, 0x6a)}
def score(a, b, last): return sum((1 + a[1][j]) * (1 + b[1][j]) * (1 + fl[last][j]) for j in range(4))
pairs = []
for i in range(len(typ)):
    for j in range(i + 1, len(typ)):
        for last in (0x51, 0x6a):
            su = SIGMA - typ[i][0] - typ[j][0] - last; x = (su - 3 * A_U) / D
            if 1.55 <= x <= 1.95: pairs.append((score(typ[i], typ[j], last), typ[i][0], typ[j][0], last))
pairs.sort(reverse=True); print('valid (pair, last) combos', len(pairs), 'score range', pairs[0][0], pairs[len(pairs) // 2][0], pairs[-1][0])
for name, sel in (('top 2%', pairs[:len(pairs) // 50]), ('top 10%', pairs[:len(pairs) // 10]), ('top 10-30%', pairs[len(pairs) // 10: 3 * len(pairs) // 10]), ('all', pairs)):
    ts = []
    for s, a, b, last in rng.sample(sel, 12): ts.append(tau([a, b, last], nsamp=300000)[0])
    print('%-10s n=%d  mean tau %.4f (min %.4f max %.4f)' % (name, len(sel), sum(ts) / len(ts), min(ts), max(ts)))
