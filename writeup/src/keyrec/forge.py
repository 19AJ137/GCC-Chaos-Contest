# Construct a key that is *equivalent* to a map-1 target (same k780 cell, k7b0 cell, S, b0b1 range, map) -- not the original password.
# n_total chunks: NF chunks are fixed text, the last 4 are solved by the decoupled systems (f3 <-> k7b0, f1 <-> k780) and the join.
import sys, random, time, math
from fractions import Fraction
from decimal import Decimal, getcontext
import numpy as np
import theory as T
from theory import *
from sums import prod_bounds
getcontext().prec = 45

def exact_parts(chunks):
    sig = sum(chunks); I1 = I3 = 0; P1 = Fraction(1); P3 = Fraction(1)
    for c in chunks:
        f, i, s, t = chunk_vals(c)
        assert s == int(i[0]) + 3 and t == int(i[2]) + 2
        I1 += int(i[0]); I3 += int(i[2]); P1 *= Fraction(1.0 + f[0]); P3 *= Fraction(1.0 + f[2])
    return sig, I1, I3, P1, P3

def residual_candidates(S, Tlo, Thi, n, fixed, nu, L3, L1, c3, h3, cells1, sig_keep=None):
    """c3,h3: centre/half-width of the k7b0 cell; cells1: list of (centre, half-width) for k780 branches"""
    sf, I1f, I3f, P1f, P3f = exact_parts(fixed)
    X, Y = {}, {}; sigs = set()
    def targets(F, kc, Pf):
        lo, hi = prod_bounds(F, nu); Pf_ = float(Pf)
        return range(max(1, math.ceil(lo * Pf_ - kc - 1e-6)), math.floor(hi * Pf_ - kc + 1e-6) + 1)
    win = sigma_window(S, Tlo, Thi, n)
    if sig_keep:
        ss = sorted(set(x[1] for x in win))
        if isinstance(sig_keep, tuple): keep = set(ss[sig_keep[0]: sig_keep[0] + sig_keep[1]])
        else: mid = len(ss) // 2; keep = set(ss[max(0, mid - sig_keep // 2): mid + sig_keep // 2])
        win = [x for x in win if x[1] in keep]
    for (w, sig, I3, T16) in win:
        su = sig - sf; I3u = I3 - I3f
        num = ME * su - (I3u << SHE)
        if num < 0 or num >= (nu << SHE): continue
        F = num / 2.0**SHE; base = (nu << L3) + ((num << L3) >> SHE)
        Js = targets(F, c3, P3f)
        if len(Js): sigs.add((w, sig))
        for J in Js:
            for d in range(nu): X.setdefault(J, {}).setdefault(base - d, set()).add((w, sig))
    for (w, sig) in sigs:
        Sw = S + (w << 32); su = sig - sf; fl1 = fl_mul(MP, SHP, sig)
        for j1 in range(n):
            I1 = fl1 - j1; SJ = Sw - 3 * n - I1
            if SJ < 4: continue
            I1u = I1 - I1f; num = MP * su - (I1u << SHP)
            if num < 0 or num >= (nu << SHP): continue
            F = num / 2.0**SHP; base = (nu << L1) + ((num << L1) >> SHP)
            for bi, (kc, kh) in enumerate(cells1):
                for J in targets(F, kc, P1f):
                    for d in range(nu): Y.setdefault((bi, J), {}).setdefault(base - d, set()).add((w, sig, SJ))
    return X, Y, sigs, (P1f, P3f)

def run_blocks(nu, L, blocks, tag):
    """blocks: (label, B as Decimal, relative tol, sorted A list)"""
    fn = os.path.join(HERE, 'tmp_%s.in' % tag)
    with open(fn, 'w') as f:
        f.write('%d\n' % len(blocks))
        for (J, B, tol, As) in blocks:
            f.write('%d %s %.6e %d\n' % (J, format(B, '.34f'), tol, len(As))); f.write(' '.join(map(str, As)) + '\n')
    t0 = time.time(); r = subprocess.run([os.path.join(HERE, 'spenum'), str(nu), str(L), fn], capture_output=True, text=True)
    sols = []
    for line in r.stdout.split('\n'):
        p = line.split()
        if len(p) >= 2 + nu: sols.append((int(p[0]), int(p[1]), tuple(int(x) for x in p[2:2 + nu])))
    return sols, r.stderr.strip().split('\n')[-1], time.time() - t0

def forge(tgt, fixed_text, table, n=10, nu=4, shrink=0.8, log=print):
    fixed = chunks_of(fixed_text); assert len(fixed) == n - nu and len(fixed_text) == 4 * (n - nu)
    c3, h3 = tgt['k7b0']; cells1 = tgt['k780']
    X, Y, sigs, (P1f, P3f) = residual_candidates(tgt['S'], tgt['Tlo'], tgt['Thi'], n, fixed, nu, 15, 15, c3, h3, cells1)
    dec = lambda fr: Decimal(fr.numerator) / Decimal(fr.denominator)
    xb = [(J, (Decimal(J) + Decimal(c3)) / dec(P3f), shrink * h3 / (J + c3), sorted(X[J])) for J in sorted(X)]
    log('  sigma candidates %d; X: %d targets, %d sums; Y: %d targets, %d sums' % (len(sigs), len(xb), sum(len(b[3]) for b in xb), len(Y), sum(len(v) for v in Y.values())))
    if not xb or not Y: return []
    xs, xstat, xt = run_blocks(nu, 15, xb, 'fx'); log('  X: %d solutions, %.0fs' % (len(xs), xt))
    ys = []
    for bi, (kc, kh) in enumerate(cells1):
        keys = [k for k in sorted(Y) if k[0] == bi]
        if not keys: continue
        s, st, t = run_blocks(nu, 15, [(k[1], (Decimal(k[1]) + Decimal(kc)) / dec(P1f), shrink * kh / (k[1] + kc), sorted(Y[k])) for k in keys], 'fy')
        ys += [((bi, J), A, Ns) for (J, A, Ns) in s]; log('  Y branch %d: %d solutions, %.0fs' % (bi, len(s), t))
    cands, npair = join(xs, ys, X, Y, table, nu); log('  join: %d index pairs -> %d chunk sets' % (npair, len(cands)))
    good = []
    for (chunks, ws, br, J3, J1, SJ) in cands:
        allc = fixed + list(chunks)
        for last in sorted(set(allc)):
            if map_of_i2(chunk_vals(last)[1][1]) != tgt['m']: continue
            rest = list(allc); rest.remove(last)
            key = b''.join(c.to_bytes(4, 'little') for c in rest + [last])
            r = key_schedule(key)
            if r['S'] != tgt['S'] or not (tgt['Tlo'] <= r['b0b1'] <= tgt['Thi']) or r['m'] != tgt['m']: continue
            if abs(r['k'][2] - c3) > h3: continue
            if min(abs(r['k'][0] - kc) - kh for (kc, kh) in cells1) > 0: continue
            good.append((key, r)); break
    return good

def text_chunks_table():
    pr = bytes(range(0x21, 0x7F)); last = b'efghijklmnopqrstuvwxyz'       # no spaces inside the solved part
    return ChunkTable(charset_chunks(pr, last))

if __name__ == '__main__':
    from cells import V, refv, full_v, ks_map1, rep5
    rep = V['k780']
    tgt = dict(S=0x3a9d3481, Tlo=0xab0f, Thi=0xabc5, m=1,
               k7b0=(V['k7b0'] + (2.906e-11 - 2.915e-11) / 2, (2.906e-11 + 2.915e-11) / 2),
               k780=[(rep + (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2), (1.0 - rep - (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2)])
    table = text_chunks_table(); print('chunk table', len(table.c))
    fixed_text = sys.argv[1].encode() if len(sys.argv) > 1 else b'equivalent key for GCC-v'
    print('fixed part: %r (%d chars)' % (fixed_text.decode(), len(fixed_text)))
    good = forge(tgt, fixed_text, table)
    for key, r in good:
        b = r['b0b1']; ks = ks_map1(r['k'][0], r['k'][2], b, r['warm'], 256) + ks_map1(r['k'][0], rep5(r['k'][2], 0xed), b, 0, 79)
        print('KEY %r  S %08x b0b1 %04x map %d warm %d k780 %.17g k7b0 %.17g  keystream 335 words identical: %s' % (
            key.decode('latin-1'), r['S'], b, r['m'], r['warm'], r['k'][0], r['k'][2], ks == refv))
    print('equivalent keys found:', len(good))
