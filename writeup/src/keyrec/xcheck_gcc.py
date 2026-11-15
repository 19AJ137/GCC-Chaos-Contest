# cross-check: gks.py (re-implementation from the spec) against the reference ../gcc6/gcc.py
import sys, os, random, struct
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'gcc6')); sys.path.insert(0, HERE)
import gcc, gks
rng = random.Random(20261114); bad = 0; n = 0; modes = [0, 0, 0]; skipped = 0
N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
for t in range(N):
    L = rng.randint(1, 40)
    kind = t % 3
    if kind == 0: key = bytes(rng.randint(0x20, 0x7E) for _ in range(L))
    elif kind == 1: key = bytes(rng.randint(1, 255) for _ in range(L))
    else: key = bytes(rng.choice(b'abcdefghijklmnopqrstuvwxyz0123456789') for _ in range(L))
    try: st = gcc.key_schedule(key)
    except Exception as e: skipped += 1; continue
    r = gks.key_schedule(key); n += 1
    a = (st.mode, st.S, st.T & 0xFFFFFFFF, struct.pack('<d', st.k780), struct.pack('<d', st.k788), struct.pack('<d', st.k7b0_saved), struct.pack('<d', st.k7b8))
    b = (r['m'], r['S'], r['T'], struct.pack('<d', r['k'][0]), struct.pack('<d', r['k'][1]), struct.pack('<d', r['k'][2]), struct.pack('<d', r['k'][3]))
    if 0 <= st.mode <= 2: modes[st.mode] += 1
    if a != b:
        bad += 1
        if bad <= 5: print('MISMATCH', key.hex(), [x == y for x, y in zip(a, b)], a[:3], b[:3])
print('keys compared %d (skipped %d), mismatches %d, modes %s' % (n, skipped, bad, modes))
