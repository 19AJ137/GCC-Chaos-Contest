# Decoupled inversion of the GCC6 key schedule: candidates from (S, b0b1 range), per-coordinate sum-product systems, join.
import math, struct, subprocess, itertools, random, sys, time, os, re
from decimal import Decimal, getcontext
import numpy as np
from gks import *
from sums import MP, SHP, ME, SHE, fl_mul, sigma_window, jrange
getcontext().prec = 40
HERE = os.path.dirname(os.path.abspath(__file__))

def lat(c, K):
    """lattice class ell and integer N (1+f = N / 2^ell) of chunk c for multiplier K"""
    a = clr(c * K); e = math.frexp(a)[1] - 1; ell = 44 - e
    f = a - math.floor(a)
    return ell, (1 << ell) + int(round(f * (1 << ell)))

def lat_np(c, K):
    a = (c.astype(np.float64) * K)
    a = (a.view(np.uint64) & np.uint64(0xFFFFFFFFFFFFFF00)).view(np.float64)
    e = np.frexp(a)[1] - 1; ell = 44 - e
    f = a - np.floor(a)
    return ell, (np.int64(1) << ell) + np.rint(np.ldexp(f, ell)).astype(np.int64)

def candidates(S, Tlo, Thi, n, L3, L1, k7b0, k780s):
    """A-candidates per J for the f3 system (X) and the f1 system (Y); each tagged with the (wrap, sigma) it came from"""
    X, Y = {}, {}; sigs = set()
    for (w, sig, I3, T16) in sigma_window(S, Tlo, Thi, n):
        num = ME * sig - (I3 << SHE); F3 = num / 2.0**SHE
        base = (n << L3) + ((num << L3) >> SHE)
        Js = jrange(F3, n, k7b0)
        if Js: sigs.add((w, sig))
        for J in Js:
            for d in range(n): X.setdefault(J, {}).setdefault(base - d, set()).add((w, sig))
    for (w, sig) in sigs:
        Sw = S + (w << 32); fl1 = fl_mul(MP, SHP, sig)
        for j1 in range(n):
            I1 = fl1 - j1; SJ = Sw - 3 * n - I1
            if SJ < 4: continue
            num = MP * sig - (I1 << SHP); F1 = num / 2.0**SHP
            base = (n << L1) + ((num << L1) >> SHP)
            for bi, k1 in enumerate(k780s):
                for J in jrange(F1, n, k1):
                    for d in range(n): Y.setdefault((bi, J), {}).setdefault(base - d, set()).add((w, sig, SJ))
    return X, Y, sigs

def run_spenum(n, L, blocks, tag, frac=None):
    """blocks: list of (J, kobs, width, sorted list of A). returns list of (J, A, [N...]) and stderr stats"""
    fn = os.path.join(HERE, 'tmp_%s.in' % tag)
    with open(fn, 'w') as f:
        f.write('%d\n' % len(blocks))
        for (J, kobs, width, As) in blocks:
            B = Decimal(J) + Decimal(kobs); tol = (width / 2) / float(B) + 2.0**-48
            f.write('%d %s %.6e %d\n' % (J, format(B, '.30f'), tol, len(As)))
            f.write(' '.join(str(a) for a in As) + '\n')
    cmd = [os.path.join(HERE, 'spenum'), str(n), str(L), fn] + ([str(frac[0]), str(frac[1])] if frac else [])
    t0 = time.time(); r = subprocess.run(cmd, capture_output=True, text=True); dt = time.time() - t0
    sols = []
    for line in r.stdout.split('\n'):
        p = line.split()
        if len(p) >= 2 + n: sols.append((int(p[0]), int(p[1]), tuple(int(x) for x in p[2:2 + n])))
    return sols, r.stderr.strip().split('\n')[-1], dt

class ChunkTable:
    """all chunks of a charset in one lattice class pair: lookup (N3, N1) -> chunks"""
    def __init__(self, chunks):
        c = np.asarray(chunks, dtype=np.uint64)
        l3, n3 = lat_np(c, K_E); l1, n1 = lat_np(c, K_PI)
        self.key = (n3.astype(np.int64) << 20) | n1.astype(np.int64)
        o = np.argsort(self.key, kind='stable'); self.key = self.key[o]; self.c = c[o]; self.l3 = l3[o]; self.l1 = l1[o]
        self.binfile = os.path.join(HERE, 'tmp_table_%d.bin' % (len(c)))
        ent = np.empty((len(c), 3), dtype=np.uint32); ent[:, 0] = n3[o]; ent[:, 1] = n1[o]; ent[:, 2] = self.c
        ent.tofile(self.binfile)
    def lookup(self, n3, n1):
        k = (n3 << 20) | n1
        a = np.searchsorted(self.key, k, 'left'); b = np.searchsorted(self.key, k, 'right')
        return [int(x) for x in self.c[a:b]]

def join(xs, ys, Xc, Yc, table, n, prefix='tmp_j', run=True):
    """xs: X solutions (J3, A3, N3 tuple), ys: Y solutions ((br,J1), A1, N1 tuple). returns candidate chunk multisets (through spjoin)"""
    tags = {}
    def tid(ws): return tags.setdefault(ws, len(tags))
    xl = []
    for (J3, A, N3s) in xs:
        for ws in Xc[J3][A]: xl.append((tid(ws), J3, N3s))
    yl = []
    for (key, A, N1s) in ys:
        for (w, sig, SJ) in Yc[key][A]: yl.append((tid((w, sig)), SJ, key[1], key[0], N1s))
    xl.sort(); yl.sort()
    fx = os.path.join(HERE, prefix + 'x.txt'); fy = os.path.join(HERE, prefix + 'y.txt')
    with open(fx, 'w') as f:
        for (t, J, Ns) in xl: f.write('%d %d %s\n' % (t, J, ' '.join(map(str, Ns))))
    with open(fy, 'w') as f:
        for (t, SJ, J, br, Ns) in yl: f.write('%d %d %d %d %s\n' % (t, SJ, J, br, ' '.join(map(str, Ns))))
    if not run: return set(), -1
    r = subprocess.run([os.path.join(HERE, 'spjoin'), str(n), table.binfile, fx, fy], capture_output=True, text=True)
    inv = {v: k for k, v in tags.items()}
    out = set()
    for line in r.stdout.split('\n'):
        p = line.split()
        if len(p) == 5 + n: out.add((tuple(sorted(int(x) for x in p[5:])), inv[int(p[0])], int(p[3]), int(p[1]), int(p[2]), int(p[4])))
    m = re.search(r'pairs (\d+)', r.stderr)
    return out, int(m.group(1)) if m else -1

def verify(chunks, obs):
    """try every chunk as the last one; return keys (bytes) whose schedule matches the observables"""
    good = []
    for last in sorted(set(chunks)):
        rest = list(chunks); rest.remove(last)
        key = b''.join(c.to_bytes(4, 'little') for c in rest + [last])
        r = key_schedule(key)
        if r['S'] != obs['S'] or not (obs['Tlo'] <= r['b0b1'] <= obs['Thi']) or r['m'] != obs['m']: continue
        if abs(r['k'][2] - obs['k7b0']) > obs['w3']: continue
        if min(abs(r['k'][0] - k) for k in obs['k780s']) > obs['w1']: continue
        if obs.get('k788') is not None and (abs(r['k'][1] - obs['k788']) > obs['w1'] or abs(r['k'][3] - obs['k7b8']) > obs['w3']): continue
        good.append(key)
    return good

def make_obs(key, rng, tw=44, w1=1e-11, w3=4e-11, map2=True):
    r = key_schedule(key)
    off = rng.randrange(tw)
    obs = dict(S=r['S'], Tlo=r['b0b1'] - off, Thi=r['b0b1'] - off + tw - 1, m=r['m'],
               k7b0=r['k'][2] + (rng.random() - 0.5) * w3 * 0.9, w3=w3, k780s=[r['k'][0] + (rng.random() - 0.5) * w1 * 0.9], w1=w1)
    if map2: obs['k788'] = r['k'][1]; obs['k7b8'] = r['k'][3]
    return obs, r

def solve_uniform(obs, n, L3, L1, table, frac=None, verbose=True):
    Xc, Yc, sigs = candidates(obs['S'], obs['Tlo'], obs['Thi'], n, L3, L1, obs['k7b0'], obs['k780s'])
    xb = [(J, obs['k7b0'], obs['w3'], sorted(Xc[J])) for J in sorted(Xc)]
    yb = [(key[1], obs['k780s'][key[0]], obs['w1'], sorted(Yc[key])) for key in sorted(Yc)]
    xs, xstat, xt = run_spenum(n, L3, xb, 'x', frac)
    ys = []; ystat = ''; yt = 0
    for bi in range(len(obs['k780s'])):
        keys = [key for key in sorted(Yc) if key[0] == bi]
        s, ystat, t = run_spenum(n, L1, [(key[1], obs['k780s'][bi], obs['w1'], sorted(Yc[key])) for key in keys], 'y', frac)
        ys += [((bi, J), A, Ns) for (J, A, Ns) in s]; yt += t
    cands, npair = join(xs, ys, Xc, Yc, table, n)
    keys = set()
    for (chunks, ws, br, J3, J1, SJ) in cands:
        for k in verify(list(chunks), obs): keys.add(k)
    info = dict(nsig=len(sigs), nXJ=len(xb), nXA=sum(len(b[3]) for b in xb), nYJ=len(yb), nYA=sum(len(b[3]) for b in yb),
                nX=len(xs), nY=len(ys), xt=xt, yt=yt, npair=npair, ncand=len(cands), xstat=xstat, ystat=ystat)
    return keys, info

def charset_chunks(first3, last):
    a = np.array(list(first3), dtype=np.uint64); l = np.array(list(last), dtype=np.uint64)
    g = np.meshgrid(a, a, a, l, indexing='ij')
    return (g[0] | (g[1] << np.uint64(8)) | (g[2] << np.uint64(16)) | (g[3] << np.uint64(24))).ravel()

if __name__ == '__main__':
    n = int(sys.argv[1]); trials = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    frac = (float(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else None
    rng = random.Random(seed)
    low = b'abcdefghijklmnopqrstuvwxyz'; last = b'efghijklmnopqrstuvwxyz'      # 4th char e..z: classes l3 = l1 = 15
    table = ChunkTable(charset_chunks(low, last))
    print('table', len(table.c), 'chunks; classes l3', set(table.l3.tolist()), 'l1', set(table.l1.tolist()))
    ok = 0
    for t in range(trials):
        key = b''.join(bytes([rng.choice(low), rng.choice(low), rng.choice(low), rng.choice(last)]) for _ in range(n))
        obs, r = make_obs(key, rng)
        keys, info = solve_uniform(obs, n, 15, 15, table, frac)
        canon = lambda k: tuple(sorted(k[i:i + 4] for i in range(0, len(k) - 4, 4))) + (k[-4:],)
        hit = canon(key) in set(canon(k) for k in keys)
        ok += hit
        print('%s map %d J %s | sig %d  X: J %d A %d sol %d (%.1fs)  Y: J %d A %d sol %d (%.1fs) | pairs %d cand %d -> keys %d %s' % (
            key.decode(), r['m'], r['J'], info['nsig'], info['nXJ'], info['nXA'], info['nX'], info['xt'], info['nYJ'], info['nYA'], info['nY'], info['yt'],
            info['npair'], info['ncand'], len(keys), 'HIT' if hit else 'MISS'), flush=True)
        if t == 0: print('   ', info['xstat']); print('   ', info['ystat'])
        if len(set(canon(k) for k in keys)) > 1: print('    distinct classes:', [k.decode() for k in keys][:6])
    print('recovered %d / %d' % (ok, trials))
