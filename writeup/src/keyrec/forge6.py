# Equivalent key for the map-1 target v: 5 fixed + 5 solved chunks.  The pair of product targets (J3, J1) is chosen so that
# many sigma candidates are shared between the two systems and the required J2+J4 is near what forged keys actually have.
import sys, time, math, collections
from forge5 import *

if __name__ == '__main__':
    from cells import V, refv, ks_map1, rep5
    rep = V['k780']
    tgt = dict(S=0x3a9d3481, Tlo=0xab0f, Thi=0xabc5, m=1,
               k7b0=(V['k7b0'] + (2.906e-11 - 2.915e-11) / 2, (2.906e-11 + 2.915e-11) / 2),
               k780=[(rep + (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2), (1.0 - rep - (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2)])
    fixed_text = sys.argv[1].encode(); fr = float(sys.argv[2]); want = int(sys.argv[3]); tolj = int(sys.argv[4]); nsig = int(sys.argv[5])
    mode = sys.argv[6] if len(sys.argv) > 6 else 'run'
    fry = float(sys.argv[7]) if len(sys.argv) > 7 else fr
    n, nu = 10, 5; fixed = chunks_of(fixed_text); c3, h3 = tgt['k7b0']; cells1 = tgt['k780']; shrink = 0.8
    X, Y, sigs, (P1f, P3f) = residual_candidates(tgt['S'], tgt['Tlo'], tgt['Thi'], n, fixed, nu, 15, 15, c3, h3, cells1)
    print('fixed %r; sigma candidates %d; X targets %d; Y targets %d' % (fixed_text.decode(), len(sigs), len(X), len(Y)), flush=True)
    xs_by = {J: collections.Counter(ws for A in X[J] for ws in X[J][A]) for J in X}
    best = None
    for key in Y:
        ycnt = collections.Counter()
        for A in Y[key]:
            for (w, sig, SJ) in Y[key][A]: ycnt[(w, sig, SJ)] += 1
        for J3 in X:
            sc = collections.Counter()
            for (w, sig, SJ), c in ycnt.items():
                if abs(SJ - key[1] - J3 - want) <= tolj and (w, sig) in xs_by[J3]: sc[(w, sig)] += c * xs_by[J3][(w, sig)]
            tot = sum(v for _, v in sc.most_common(nsig))
            if best is None or tot > best[0]: best = (tot, J3, key, set(ws for ws, _ in sc.most_common(nsig)))
    tot, J3, key, keep = best
    print('chosen targets: J3=%d (residual product %.4f), branch %d J1=%d (residual %.4f); %d sigma values, %d (A3,A1) combinations with required J2+J4 within %d of %d' % (
        J3, (J3 + c3) / float(P3f), key[0], key[1], (key[1] + cells1[key[0]][0]) / float(P1f), len(keep), tot, tolj, want), flush=True)
    XA = sorted(A for A in X[J3] if X[J3][A] & keep)
    YA = sorted(A for A in Y[key] if any((w, s) in keep and abs(SJ - key[1] - J3 - want) <= tolj for (w, s, SJ) in Y[key][A]))
    print('sums: X %d, Y %d' % (len(XA), len(YA)), flush=True)
    if mode == 'plan': sys.exit()
    dec = lambda f: Decimal(f.numerator) / Decimal(f.denominator)
    xs, xstat, xt = run_blocks_f(nu, 15, [(J3, (Decimal(J3) + Decimal(c3)) / dec(P3f), shrink * h3 / (J3 + c3), XA)], 'hx', (0.0, fr)); print('X: %d solutions, %.0fs | %s' % (len(xs), xt, xstat), flush=True)
    kc, kh = cells1[key[0]]
    ys, ystat, yt = run_blocks_f(nu, 15, [(key[1], (Decimal(key[1]) + Decimal(kc)) / dec(P1f), shrink * kh / (key[1] + kc), YA)], 'hy', (0.0, fry)); print('Y: %d solutions, %.0fs | %s' % (len(ys), yt, ystat), flush=True)
    ys = [((key[0], J), A, Ns) for (J, A, Ns) in ys]
    class Tab: binfile = os.path.join(HERE, 'tmp_table_kana.bin')
    if mode == 'enum': join(xs, ys, X, Y, Tab, nu, prefix='tmp_h_j', run=False); print('solutions written to tmp_h_jx.txt / tmp_h_jy.txt', flush=True); sys.exit()
    cands, npair = join(xs, ys, X, Y, Tab, nu, prefix='tmp_h_j'); print('join: %d index pairs -> %d chunk sets' % (npair, len(cands)), flush=True)
    good = []; sum_ok = 0; sf = sum(fixed)
    for (chunks, ws, br, j3, j1, SJ) in cands:
        allc = fixed + list(chunks); sum_ok += (sum(chunks) + sf == ws[1])
        last = [c for c in sorted(set(chunks)) if map_of_i2(chunk_vals(c)[1][1]) == 1] or [c for c in sorted(set(fixed)) if map_of_i2(chunk_vals(c)[1][1]) == 1]
        if not last: continue
        rest = list(allc); rest.remove(last[0]); k = b''.join(c.to_bytes(4, 'little') for c in rest + [last[0]]); r = key_schedule(k)
        inc = abs(r['k'][2] - c3) <= h3 and min(abs(r['k'][0] - a) - b for (a, b) in cells1) <= 0
        if r['S'] == tgt['S'] and tgt['Tlo'] <= r['b0b1'] <= tgt['Thi'] and r['m'] == 1 and inc: good.append((k, r))
    print('chunk sets with the right chunk sum: %d; fully matching keys: %d' % (sum_ok, len(good)))
    seenk = set()
    for k, r in good:
        if k in seenk: continue
        seenk.add(k); b = r['b0b1']; ks = ks_map1(r['k'][0], r['k'][2], b, r['warm'], 256) + ks_map1(r['k'][0], rep5(r['k'][2], 0xed), b, 0, 79)
        print('KEY hex %s | cp932 %r | len %d S %08x b0b1 %04x map %d warm %d k780 %.17g k7b0 %.17g J %s | 335 keystream words identical: %s' % (
            k.hex(), k.decode('cp932', 'replace'), len(k), r['S'], b, r['m'], r['warm'], r['k'][0], r['k'][2], r['J'], ks == refv))
