# Equivalent key for a map-1 target, 5 fixed chunks + 5 solved chunks, one product target per system (all its sum candidates),
# optional slices of the enumeration.  Not the original password.
import sys, time, math
from forge import *

def run_blocks_f(nu, L, blocks, tag, frac):
    fn = os.path.join(HERE, 'tmp_%s.in' % tag)
    with open(fn, 'w') as f:
        f.write('%d\n' % len(blocks))
        for (J, B, tol, As) in blocks:
            f.write('%d %s %.6e %d\n' % (J, format(B, '.34f'), tol, len(As))); f.write(' '.join(map(str, As)) + '\n')
    t0 = time.time(); r = subprocess.run([os.path.join(HERE, 'spenum'), str(nu), str(L), fn, str(frac[0]), str(frac[1])], capture_output=True, text=True)
    sols = []
    for line in r.stdout.split('\n'):
        p = line.split()
        if len(p) >= 2 + nu: sols.append((int(p[0]), int(p[1]), tuple(int(x) for x in p[2:2 + nu])))
    return sols, r.stderr.strip().split('\n')[-1], time.time() - t0

def forge_single(tgt, fixed_text, table, n=10, nu=5, fx=(0.0, 0.3), fy=(0.0, 0.3), shrink=0.8, nJ=1, sig_keep=None, log=print, force=None):
    fixed = chunks_of(fixed_text); assert len(fixed) == n - nu
    c3, h3 = tgt['k7b0']; cells1 = tgt['k780']
    X, Y, sigs, (P1f, P3f) = residual_candidates(tgt['S'], tgt['Tlo'], tgt['Thi'], n, fixed, nu, 15, 15, c3, h3, cells1, sig_keep)
    dec = lambda fr: Decimal(fr.numerator) / Decimal(fr.denominator)
    log('  sigma candidates %d; X targets %d, Y targets %d; fixed products P1 %.4f P3 %.4f' % (len(sigs), len(X), len(Y), float(P1f), float(P3f)))
    if not X or not Y: return [], None
    # pick the targets with the most sum candidates
    Jx = sorted(X, key=lambda J: -len(X[J]))[:nJ]
    Jy = sorted(Y, key=lambda k: -len(Y[k]))[:nJ]
    if force: Jx = [force[0]]; Jy = [(0, force[1])]
    sjs = sorted(t[2] for k in Jy for A in Y[k] for t in Y[k][A]); log('  required J2+J4 for these targets: %d..%d (median %d)' % (sjs[0] - Jx[0] - Jy[0][1], sjs[-1] - Jx[0] - Jy[0][1], sjs[len(sjs) // 2] - Jx[0] - Jy[0][1]))
    for J in Jx: log('  X target J3=%d: residual product %.4f, %d sums' % (J, (J + c3) / float(P3f), len(X[J])))
    for k in Jy: log('  Y target branch %d J1=%d: residual product %.4f, %d sums' % (k[0], k[1], (k[1] + cells1[k[0]][0]) / float(P1f), len(Y[k])))
    xs, xstat, xt = run_blocks_f(nu, 15, [(J, (Decimal(J) + Decimal(c3)) / dec(P3f), shrink * h3 / (J + c3), sorted(X[J])) for J in Jx], 'gx', fx)
    log('  X: %d solutions, %.0fs | %s' % (len(xs), xt, xstat))
    ys = []
    for k in Jy:
        kc, kh = cells1[k[0]]
        s, st, t = run_blocks_f(nu, 15, [(k[1], (Decimal(k[1]) + Decimal(kc)) / dec(P1f), shrink * kh / (k[1] + kc), sorted(Y[k]))], 'gy', fy)
        ys += [((k[0], J), A, Ns) for (J, A, Ns) in s]; log('  Y: %d solutions, %.0fs | %s' % (len(s), t, st))
    t0 = time.time(); cands, npair = join(xs, ys, X, Y, table, nu); log('  join: %d index pairs -> %d chunk sets (%.0fs)' % (npair, len(cands), time.time() - t0))
    good = []; stat = dict(sum_ok=0, cells_ok=0, n=len(cands))
    sf = sum(fixed)
    for (chunks, ws, br, J3, J1, SJ) in cands:
        allc = fixed + list(chunks)
        if sum(chunks) + sf == ws[1]: stat['sum_ok'] += 1
        last = [c for c in sorted(set(allc)) if map_of_i2(chunk_vals(c)[1][1]) == tgt['m']]
        if not last: continue
        rest = list(allc); rest.remove(last[0])
        key = b''.join(c.to_bytes(4, 'little') for c in rest + [last[0]])
        r = key_schedule(key)
        inc = abs(r['k'][2] - c3) <= h3 and min(abs(r['k'][0] - kc) - kh for (kc, kh) in cells1) <= 0
        stat['cells_ok'] += inc
        if r['S'] == tgt['S'] and tgt['Tlo'] <= r['b0b1'] <= tgt['Thi'] and r['m'] == tgt['m'] and inc: good.append((key, r))
    return good, stat

if __name__ == '__main__':
    from cells import V, refv, ks_map1, rep5
    rep = V['k780']
    tgt = dict(S=0x3a9d3481, Tlo=0xab0f, Thi=0xabc5, m=1,
               k7b0=(V['k7b0'] + (2.906e-11 - 2.915e-11) / 2, (2.906e-11 + 2.915e-11) / 2),
               k780=[(rep + (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2), (1.0 - rep - (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2)])
    fixed_text = sys.argv[1].encode(); fr = float(sys.argv[2]); sk = (int(sys.argv[3]), int(sys.argv[4]))
    force = (int(sys.argv[5]), int(sys.argv[6])) if len(sys.argv) > 6 else None
    class Tab: binfile = os.path.join(HERE, sys.argv[7] if len(sys.argv) > 7 else 'tmp_table_kana.bin')
    print('fixed part: %r (%d chars), sigma_fix %.4e' % (fixed_text.decode(), len(fixed_text), sum(chunks_of(fixed_text))), flush=True)
    good, stat = forge_single(tgt, fixed_text, Tab, fx=(0.0, fr), fy=(0.0, fr), sig_keep=sk, log=lambda *a: print(*a, flush=True), force=force)
    print('joint candidates:', stat)
    for key, r in good:
        b = r['b0b1']; ks = ks_map1(r['k'][0], r['k'][2], b, r['warm'], 256) + ks_map1(r['k'][0], rep5(r['k'][2], 0xed), b, 0, 79)
        print('KEY hex %s cp932 %r len %d  S %08x b0b1 %04x map %d warm %d k780 %.17g k7b0 %.17g J %s | 335 keystream words identical: %s' % (
            key.hex(), key.decode('cp932', 'replace'), len(key), r['S'], b, r['m'], r['warm'], r['k'][0], r['k'][2], r['J'], ks == refv))
    print('equivalent keys found:', len(good))
