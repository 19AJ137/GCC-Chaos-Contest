# 対照 1: 実物の鍵 v の目標に対して、等価な鍵（forge で作った 40 バイト）の一部を「分かっている」として残りを解く
import sys, time
import vpart
from vpart import build_hyp, run_vres, TARGET_V, keystream_ok_v, show
FK = bytes.fromhex('6571756976616c656e74204b45593a20474343766e662e694b37586d467a7875df74667668347e69')
L = int(sys.argv[1]) if len(sys.argv) > 1 else 16
cases = [
 ('c0_all',      [FK]),
 ('c1_last4',    [FK[:36], 4]),
 ('c2_last8',    [FK[:32], 8]),
 ('c3_last12',   [FK[:28], 12]),
 ('c4_first12',  [12, FK[12:]]),
 ('c5_mid12',    [FK[:8], 12, FK[20:]]),
 ('c6_mid8',     [FK[:16], 8, FK[24:]]),
 ('c7_8+1',      [FK[:27], 9, FK[36:]]),
 ('c8_12+1',     [FK[:27], 13]),
 ('c9_2byte',    [FK[:38], 2]),
 ('c10_split',   [FK[:8], 4, FK[12:20], 4, FK[24:32], 4, FK[36:]]),
]
hyps = [build_hyp(i, s) for i, s in cases]
for h in hyps: print(h['id'], 'nU', h['nU'], 'E', h['E'], show(h))
keys, stat, dt = run_vres(TARGET_V, hyps, L=L, tag='ctl1', keep=True)
print('time %.1f s' % dt)
ok = True
for h in hyps:
    ks = keys.get(h['id'], set()); good = [k for k in ks if keystream_ok_v(k)]
    hit = FK in ks
    print('%-10s keys %d  keystream-ok %d  original %s  %s' % (h['id'], len(ks), len(good), 'HIT' if hit else 'MISS', stat.get(h['id'])))
    for k in sorted(ks):
        if k != FK: print('     other:', k.hex(), keystream_ok_v(k))
    ok &= hit
print('ALL OK' if ok else 'FAILED')
