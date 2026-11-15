import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import random, collections
from emu import Emu
import gcc
e = Emu()
rnd = random.Random(1997)
cnt = collections.Counter(); bad = 0; N = 600
for t in range(N):
    kl = rnd.choice([1,2,3,4,5,6,7,8,9,12,16,20,31,40])
    if rnd.random() < 0.5:
        key = bytes(rnd.choice(b'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789') for _ in range(kl))
    else:
        key = bytes(rnd.randrange(1,256) for _ in range(kl))
    dl = rnd.choice([1,3,4,7,64,255,1024,1500])
    data = bytes(rnd.randrange(256) for _ in range(dl))
    try:
        st = gcc.key_schedule(key)
    except ValueError as ex:
        cnt['py-err'] += 1; continue
    cnt[st.mode] += 1
    ref = e.str_cryption(key, data)
    mine = gcc.str_cryption(key, data)
    if ref != mine:
        bad += 1
        if bad < 5: print('MISMATCH', key, dl, st.mode, ref[:16].hex(), mine[:16].hex())
    else:
        assert gcc.str_cryption(key, mine, encrypt=False) == data
print('tests', N, 'mismatches', bad, 'modes', dict(cnt))
