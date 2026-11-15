# keys of 3 full chunks + a tail of 1..3 bytes (13..15 bytes): tail3 (enumeration) + spjoin2 (chunk reconstruction) + full check
import sys, os, random, subprocess, time, itertools
from gks import *
from kall import K_OBS, synth_key, synth_obs
HERE = os.path.dirname(os.path.abspath(__file__))

def run_tail(obs, taillen, blo=0x20, bhi=0xFC, part=None, log=print):
    cmd = [os.path.join(HERE, 'tail3'), '%x' % obs['S'], '%x' % obs['Tlo'], '%x' % obs['Thi'], str(obs['m']), repr(obs['k7b0']), repr(obs['h3']),
           repr(obs['k780s'][0]), repr(obs['h1']), str(taillen), str(blo), str(bhi)] + ([str(part[0]), str(part[1])] if part else [])
    t0 = time.time(); r = subprocess.run(cmd, capture_output=True, text=True)
    X, Y = [], []
    for line in r.stdout.split('\n'):
        p = line.split()
        if not p: continue
        if p[0] == 'X': X.append((int(p[1]), int(p[2]), int(p[3]), tuple(int(v) for v in p[5:8])))
        elif p[0] == 'Y': Y.append((int(p[1]), int(p[2]), int(p[3]), int(p[4]), tuple(int(v) for v in p[6:9])))
    log('  tail length %d%s: %.0fs | %s | X %d, Y %d solutions' % (taillen, ' part %d/%d' % part if part else '', time.time() - t0, r.stderr.strip().split('\n')[-1], len(X), len(Y)))
    return X, Y

def finish(obs, taillen, X, Y, log=print, tagname='kt'):
    gx = {}; gy = {}
    for (t, sig, J, Ns) in X: gx.setdefault((t, sig), []).append((J, Ns))
    for (t, sig, SJ, J, Ns) in Y: gy.setdefault((t, sig), []).append((SJ, J, Ns))
    common = sorted(set(gx) & set(gy))
    log('  groups (tail, sigma): X %d, Y %d, both %d' % (len(gx), len(gy), len(common)))
    if not common: return set()
    fx = os.path.join(HERE, 'tmp_%s_jx.txt' % tagname); fy = os.path.join(HERE, 'tmp_%s_jy.txt' % tagname); fp = os.path.join(HERE, 'tmp_%s_par.txt' % tagname)
    with open(fx, 'w') as f:
        for i, g in enumerate(common):
            for (J, Ns) in gx[g]: f.write('%d %d %s\n' % (i, J, ' '.join(map(str, Ns))))
    with open(fy, 'w') as f:
        for i, g in enumerate(common):
            for (SJ, J, Ns) in gy[g]: f.write('%d %d %d 0 %s\n' % (i, SJ, J, ' '.join(map(str, Ns))))
    with open(fp, 'w') as f:
        f.write('%x %x %x %.17g %.17g %.17g %.17g %.3e %d %d\n' % (obs['S'], obs['Tlo'], obs['Thi'], obs['k780s'][0], -1.0, obs['k7b0'], -1.0, -1.0, 1, len(common)))
        for i, g in enumerate(common): f.write('%d %d\n' % (i, g[1]))
    r = subprocess.run([os.path.join(HERE, 'spjoin2'), '3', fp, fx, fy], capture_output=True, text=True)
    log('  join: %s' % r.stderr.strip().split('\n')[-1])
    keys = set()
    for line in r.stdout.split('\n'):
        p = line.split()
        if len(p) < 4: continue
        tail = common[int(p[0])][0]; chunks = [int(v) for v in p[1:4]]
        key = b''.join(c.to_bytes(4, 'little') for c in chunks) + tail.to_bytes(taillen, 'little'); rr = key_schedule(key)
        if rr['S'] != obs['S'] or not (obs['Tlo'] <= rr['b0b1'] <= obs['Thi']) or rr['m'] != obs['m']: continue
        if abs(rr['k'][2] - obs['k7b0']) > obs['h3'] or abs(rr['k'][0] - obs['k780s'][0]) > obs['h1']: continue
        if obs.get('k788') is not None and (abs(rr['k'][1] - obs['k788']) > obs['h2'] or abs(rr['k'][3] - obs['k7b8']) > obs['h4']): continue
        keys.add(key)
    return keys

if __name__ == '__main__':
    mode = sys.argv[1]; taillen = int(sys.argv[2])
    log = lambda *a: print(*a, flush=True)
    if mode == 'control':
        seed = int(sys.argv[3]); rng = random.Random(seed)
        while True:
            key = synth_key(rng, 3, sys.argv[4] if len(sys.argv) > 4 else 'sjis') + bytes(rng.randint(0x20, 0xFC) for _ in range(taillen))
            if sum(chunks_of(key)[:3]) > 0: break
        obs, r = synth_obs(key, rng)
        log('control key %s (%d bytes) map %d J %s' % (key.hex(), len(key), r['m'], r['J']))
        if taillen == 3:      # only the slice of tails that contains the true one (1/64 of the tails)
            nb = 0xFC - 0x20 + 1; b = key[12:]; ti = (b[0] - 0x20) + (b[1] - 0x20) * nb + (b[2] - 0x20) * nb * nb
            X, Y = run_tail(obs, taillen, part=(ti % 64, 64), log=log)
        else: X, Y = run_tail(obs, taillen, log=log)
        keys = finish(obs, taillen, X, Y, log=log, tagname='ktc%d' % seed)
        log('  -> %d keys: %s %s' % (len(keys), 'HIT' if key in keys or any(sorted(k[i:i+4] for i in (0, 4, 8)) == sorted(key[i:i+4] for i in (0, 4, 8)) and k[12:] == key[12:] for k in keys) else 'MISS', [k.hex() for k in sorted(keys)][:3]))
    else:
        parts = [(int(sys.argv[3]), int(sys.argv[4]))] if len(sys.argv) > 4 else [None]
        log('contest key k: 3 full chunks + %d-byte tail (%d bytes); tail bytes 0x20..0xFC, every 4th byte of the full chunks >= 0x33' % (taillen, 12 + taillen))
        X, Y = [], []
        for part in parts:
            x, y = run_tail(K_OBS, taillen, part=part, log=log); X += x; Y += y
        keys = finish(K_OBS, taillen, X, Y, log=log, tagname='ktr%d' % taillen)
        log('  -> keys matching every target: %d' % len(keys))
        for k in sorted(keys): log('     KEY hex %s  cp932 %r' % (k.hex(), k.decode('cp932', 'replace')))
