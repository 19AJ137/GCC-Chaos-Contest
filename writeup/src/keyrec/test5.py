# n = 5 (20 characters) scaled-down test: the smallest element of each system is confined to a slice of its range that contains the
# true value (1/SL of the full enumeration); everything else is as in the blind search.
import sys, random, time
import theory as T
from theory import *

SL = int(sys.argv[1]); seed = int(sys.argv[2]); trials = int(sys.argv[3]) if len(sys.argv) > 3 else 1
rng = random.Random(seed)
low = b'abcdefghijklmnopqrstuvwxyz'; last = b'efghijklmnopqrstuvwxyz'
table = ChunkTable(charset_chunks(low, last))
n = 5
for t in range(trials):
    key = b''.join(bytes([rng.choice(low), rng.choice(low), rng.choice(low), rng.choice(last)]) for _ in range(n))
    obs, r = make_obs(key, rng)
    ch = chunks_of(key)
    n3 = sorted(lat(c, K_E)[1] for c in ch); n1 = sorted(lat(c, K_PI)[1] for c in ch)
    Xc, Yc, sigs = candidates(obs['S'], obs['Tlo'], obs['Thi'], n, 15, 15, obs['k7b0'], obs['k780s'])
    def frac_for(Ac, J, nmin):
        amax = max(Ac[J]) // n; lo = 1 << 15
        f = (nmin - lo) / (amax + 1 - lo); a = max(0.0, f - 0.5 / SL); return (a, a + 1.0 / SL)
    J3 = r['J'][2]; J1 = r['J'][0]
    fx = frac_for(Xc, J3, n3[0]); fy = frac_for(Yc, (0, J1), n1[0])
    xb = [(J, obs['k7b0'], obs['w3'], sorted(Xc[J])) for J in sorted(Xc)]
    yb = [(k[1], obs['k780s'][0], obs['w1'], sorted(Yc[k])) for k in sorted(Yc)]
    print('key %s map %d J %s | sigma cand %d | X: %d J, %d A | Y: %d J, %d A' % (key.decode(), r['m'], r['J'], len(sigs), len(xb), sum(len(b[3]) for b in xb), len(yb), sum(len(b[3]) for b in yb)), flush=True)
    xs, xstat, xt = run_spenum(n, 15, xb, 'x5', fx); print('  X slice 1/%d: %d sol, %.0fs | %s' % (SL, len(xs), xt, xstat), flush=True)
    ys, ystat, yt = run_spenum(n, 15, yb, 'y5', fy); print('  Y slice 1/%d: %d sol, %.0fs | %s' % (SL, len(ys), yt, ystat), flush=True)
    truex = tuple(n3) in set(s[2] for s in xs); truey = tuple(n1) in set(s[2] for s in ys)
    ys2 = [((0, J), A, Ns) for (J, A, Ns) in ys]
    t0 = time.time(); cands, npair = join(xs, ys2, Xc, Yc, table, n); jt = time.time() - t0
    keys = set()
    for (chunks, ws, br, j3, j1, SJ) in cands:
        for k in verify(list(chunks), obs): keys.add(k)
    canon = lambda k: tuple(sorted(k[i:i + 4] for i in range(0, len(k) - 4, 4))) + (k[-4:],)
    print('  true N3 set found %s, true N1 set found %s | join: %d pairs (%.0fs) -> %d chunk sets -> %d keys (%d classes) %s' % (
        truex, truey, npair, jt, len(cands), len(keys), len(set(canon(k) for k in keys)), 'HIT' if canon(key) in set(canon(k) for k in keys) else 'MISS'), flush=True)
