# Inversion for keys of n full chunks with arbitrary bytes (every 4th byte >= 0x33): both systems on the 2^-16 lattice,
# join through spjoin2 (fibres + chunk-sum + key-schedule check).  Synthetic controls and the contest key k.
import sys, os, random, time, subprocess, math
from decimal import Decimal
from gks import *
from sums import MP, SHP, ME, SHE, fl_mul, sigma_window, jrange
HERE = os.path.dirname(os.path.abspath(__file__))
L = 16

def cands16(S, Tlo, Thi, n, k7b0, k780s, d3, d1):
    X, Y = {}, {}; sigs = set()
    for (w, sig, I3, T16) in sigma_window(S, Tlo, Thi, n):
        num = ME * sig - (I3 << SHE); F3 = num / 2.0**SHE
        base = (n << L) + ((num << L) >> SHE)
        Js = jrange(F3, n, k7b0)
        if Js: sigs.add((w, sig))
        for J in Js:
            for d in range(d3): X.setdefault(J, {}).setdefault(base - d, set()).add((w, sig))
    for (w, sig) in sigs:
        Sw = S + (w << 32); fl1 = fl_mul(MP, SHP, sig)
        for j1 in range(n):
            I1 = fl1 - j1; SJ = Sw - 3 * n - I1
            if SJ < 4: continue
            num = MP * sig - (I1 << SHP); F1 = num / 2.0**SHP
            base = (n << L) + ((num << L) >> SHP)
            for bi, k1 in enumerate(k780s):
                for J in jrange(F1, n, k1):
                    for d in range(d1): Y.setdefault((bi, J), {}).setdefault(base - d, set()).add((w, sig, SJ))
    return X, Y, sigs

def run_spenum(n, blocks, tag):
    fn = os.path.join(HERE, 'tmp_%s.in' % tag)
    with open(fn, 'w') as f:
        f.write('%d\n' % len(blocks))
        for (J, kobs, half, As) in blocks:
            B = Decimal(J) + Decimal(kobs)
            f.write('%d %s %.6e %d\n' % (J, format(B, '.30f'), half / float(B) + 2.0**-48, len(As))); f.write(' '.join(map(str, As)) + '\n')
    t0 = time.time(); r = subprocess.run([os.path.join(HERE, 'spenum'), str(n), str(L), fn], capture_output=True, text=True)
    sols = []
    for line in r.stdout.split('\n'):
        p = line.split()
        if len(p) == 2 + n: sols.append((int(p[0]), int(p[1]), tuple(int(x) for x in p[2:])))
    return sols, r.stderr.strip().split('\n')[-1], time.time() - t0

def solve(obs, n, d3=13, d1=11, log=print, tagname='ka'):
    """obs: S, Tlo, Thi, m, k7b0, h3 (half-width), k780s, h1, k788, k7b8 (may be None), h2, h4"""
    X, Y, sigs = cands16(obs['S'], obs['Tlo'], obs['Thi'], n, obs['k7b0'], obs['k780s'], d3, d1)
    log('  sigma candidates %d; X: %d J, %d sums; Y: %d J, %d sums' % (len(sigs), len(X), sum(len(v) for v in X.values()), len(Y), sum(len(v) for v in Y.values())))
    xs, xstat, xt = run_spenum(n, [(J, obs['k7b0'], obs['h3'], sorted(X[J])) for J in sorted(X)], tagname + 'x'); log('  X: %d solutions, %.0fs | %s' % (len(xs), xt, xstat))
    ys = []
    for bi in range(len(obs['k780s'])):
        keys = [k for k in sorted(Y) if k[0] == bi]
        s, st, t = run_spenum(n, [(k[1], obs['k780s'][bi], obs['h1'], sorted(Y[k])) for k in keys], tagname + 'y')
        ys += [(bi, J, A, Ns) for (J, A, Ns) in s]; log('  Y (branch %d): %d solutions, %.0fs | %s' % (bi, len(s), t, st))
    tags = {}
    def tid(ws): return tags.setdefault(ws, len(tags))
    xl = sorted((tid(ws), J, Ns) for (J, A, Ns) in xs for ws in X[J][A])
    yl = sorted((tid((w, sig)), SJ, J, bi, Ns) for (bi, J, A, Ns) in ys for (w, sig, SJ) in Y[(bi, J)][A])
    fx = os.path.join(HERE, 'tmp_%s_jx.txt' % tagname); fy = os.path.join(HERE, 'tmp_%s_jy.txt' % tagname); fp = os.path.join(HERE, 'tmp_%s_par.txt' % tagname)
    with open(fx, 'w') as f:
        for (t, J, Ns) in xl: f.write('%d %d %s\n' % (t, J, ' '.join(map(str, Ns))))
    with open(fy, 'w') as f:
        for (t, SJ, J, br, Ns) in yl: f.write('%d %d %d %d %s\n' % (t, SJ, J, br, ' '.join(map(str, Ns))))
    with open(fp, 'w') as f:
        k788 = obs.get('k788'); k7b8 = obs.get('k7b8')
        f.write('%x %x %x %.17g %.17g %.17g %.17g %.3e %d %d\n' % (obs['S'], obs['Tlo'], obs['Thi'], obs['k780s'][0], k788 if k788 is not None else -1.0, obs['k7b0'],
                                                             k7b8 if k7b8 is not None else -1.0, 1e-9, len(obs['k780s']), len(tags)))
        for ws, t in tags.items(): f.write('%d %d\n' % (t, ws[1]))
    t0 = time.time(); r = subprocess.run([os.path.join(HERE, 'spjoin2'), str(n), fp, fx, fy], capture_output=True, text=True)
    log('  join: %.0fs | %s' % (time.time() - t0, r.stderr.strip().split('\n')[-1]))
    keys = set()
    for line in r.stdout.split('\n'):
        p = line.split()
        if len(p) < 1 + n: continue
        chunks = [int(x) for x in p[1:1 + n]]
        for last in sorted(set(chunks)):
            rest = list(chunks); rest.remove(last)
            key = b''.join(c.to_bytes(4, 'little') for c in rest + [last]); rr = key_schedule(key)
            if rr['S'] != obs['S'] or not (obs['Tlo'] <= rr['b0b1'] <= obs['Thi']) or rr['m'] != obs['m']: continue
            if abs(rr['k'][2] - obs['k7b0']) > obs['h3'] or min(abs(rr['k'][0] - k) for k in obs['k780s']) > obs['h1']: continue
            if obs.get('k788') is not None and (abs(rr['k'][1] - obs['k788']) > obs['h2'] or abs(rr['k'][3] - obs['k7b8']) > obs['h4']): continue
            keys.add(key)
    return keys

def synth_key(rng, n, kind):
    out = b''
    for _ in range(n):
        if kind == 'sjis':      # two double-byte characters: lead 0x81-0x9F/0xE0-0xEA, trail 0x40-0x7E/0x80-0xFC
            lead = lambda: rng.choice(list(range(0x81, 0xA0)) + list(range(0xE0, 0xEB))); trail = lambda: rng.choice(list(range(0x40, 0x7F)) + list(range(0x80, 0xFD)))
            out += bytes([lead(), trail(), lead(), trail()])
        elif kind == 'hira':    # hiragana: lead 0x82, trail 0x9F-0xF1
            out += bytes([0x82, rng.randint(0x9F, 0xF1), 0x82, rng.randint(0x9F, 0xF1)])
        elif kind == 'kana':    # half-width kana
            out += bytes(rng.randint(0xA1, 0xDF) for _ in range(4))
        elif kind == 'mix':     # any non-zero bytes, 4th byte >= 0x33
            out += bytes([rng.randint(1, 255), rng.randint(1, 255), rng.randint(1, 255), rng.randint(0x33, 0xFF)])
    return out

def synth_obs(key, rng, tw=44, w1=4.4e-12, w3=4.0e-11):
    r = key_schedule(key); off = rng.randrange(tw)
    return dict(S=r['S'], Tlo=r['b0b1'] - off, Thi=r['b0b1'] - off + tw - 1, m=r['m'], k7b0=r['k'][2] + (rng.random() - 0.5) * w3 * 0.9, h3=w3 / 2,
                k780s=[r['k'][0] + (rng.random() - 0.5) * w1 * 0.9], h1=w1 / 2, k788=r['k'][1], k7b8=r['k'][3], h2=1e-11, h4=5e-11), r

K_OBS = dict(S=0xa678044d, Tlo=0x7ec1, Thi=0x7eec, m=2, k7b0=h2d('3fc3ca0d1ca9ebf4'), h3=2.4e-11, k780s=[h2d('3fe500287933f358')], h1=2.6e-12,
             k788=h2d('3fe4f462e3890aef'), k7b8=h2d('3fe924aef33bd88d'), h2=9e-12, h4=4.8e-11)

if __name__ == '__main__':
    mode = sys.argv[1]; n = int(sys.argv[2])
    canon = lambda k: tuple(sorted(k[i:i + 4] for i in range(0, len(k) - 4, 4))) + (k[-4:],)
    if mode == 'control':
        kind = sys.argv[3]; seed = int(sys.argv[4]); rng = random.Random(seed)
        key = synth_key(rng, n, kind); obs, r = synth_obs(key, rng)
        cl = [(44 - (math.frexp(clr(c * K_E))[1] - 1), 44 - (math.frexp(clr(c * K_PI))[1] - 1)) for c in chunks_of(key)]
        print('control %s key %s map %d J %s classes(f3,f1) %s sigma %d' % (kind, key.hex(), r['m'], r['J'], cl, sum(chunks_of(key))), flush=True)
        keys = solve(obs, n, log=lambda *a: print(*a, flush=True), tagname='kc%d' % seed)
        print('  -> %d keys (%d classes): %s' % (len(keys), len(set(canon(k) for k in keys)), 'HIT' if canon(key) in set(canon(k) for k in keys) else 'MISS'), sorted(k.hex() for k in keys)[:4])
    else:
        print('contest key k, %d full chunks (%d bytes), every 4th byte >= 0x33' % (n, 4 * n), flush=True)
        keys = solve(K_OBS, n, log=lambda *a: print(*a, flush=True), tagname='kreal%d' % n)
        print('  -> keys matching every target: %d' % len(keys))
        for k in sorted(keys): print('     KEY hex %s  cp932 %r' % (k.hex(), k.decode('cp932', 'replace')))
