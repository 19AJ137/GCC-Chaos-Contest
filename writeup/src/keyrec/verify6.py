# verify joint candidates (output of spjoin) of the equivalent-key construction against the v target and its keystream
import sys
from gks import *
from cells import V, refv, ks_map1, rep5
fixed = chunks_of(sys.argv[2].encode()); S_t = 0x3a9d3481; rep = V['k780']
c3, h3 = V['k7b0'] + (2.906e-11 - 2.915e-11) / 2, (2.906e-11 + 2.915e-11) / 2
cells1 = [(rep + (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2), (1.0 - rep - (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2)]
seen = set(); n = cells = small = 0; good = []
data = open(sys.argv[1]).read().split('\n')[:-1]
for line in data:
    p = line.split()
    if len(p) != 10: continue
    ch = tuple(sorted(int(x) for x in p[5:]))
    if ch in seen: continue
    seen.add(ch); n += 1
    allc = fixed + list(ch)
    last = [c for c in sorted(set(ch)) if map_of_i2(chunk_vals(c)[1][1]) == 1] or [c for c in sorted(set(fixed)) if map_of_i2(chunk_vals(c)[1][1]) == 1]
    if not last: continue
    rest = list(allc); rest.remove(last[0])
    key = b''.join(c.to_bytes(4, 'little') for c in rest + [last[0]])
    r = key_schedule(key)
    inc = abs(r['k'][2] - c3) <= h3 and min(abs(r['k'][0] - kc) - kh for (kc, kh) in cells1) <= 0
    cells += inc
    d = (r['S'] - S_t) & 0xFFFFFFFF
    if d >= 2**31: d -= 2**32
    small += abs(d) < 400
    if d == 0 and inc and 0xab0f <= r['b0b1'] <= 0xabc5 and r['m'] == 1: good.append((key, r))
print('distinct chunk sets %d, inside both cells %d, right chunk sum (|dS|<400) %d, fully matching keys %d' % (n, cells, small, len(good)))
for key, r in good[:20]:
    b = r['b0b1']; ks = ks_map1(r['k'][0], r['k'][2], b, r['warm'], 256) + ks_map1(r['k'][0], rep5(r['k'][2], 0xed), b, 0, 79)
    print('KEY hex %s | cp932 %r | len %d S %08x b0b1 %04x map %d warm %d k780 %.17g k7b0 %.17g J %s | 335 words identical: %s' % (
        key.hex(), key.decode('cp932', 'replace'), len(key), r['S'], b, r['m'], r['warm'], r['k'][0], r['k'][2], r['J'], ks == refv))
