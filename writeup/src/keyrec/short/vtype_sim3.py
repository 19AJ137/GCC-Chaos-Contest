# 21 バイト: 固定する 2 塊の選び方（打てる塊ぜんぶ／半角 1 バイト文字 4 つの塊だけ）で、打てる割合がどう変わるか
import sys, os, math, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from vpart import textlike
from vshort import load_hi, SIGMA, A_U, B_U
b4 = lambda c: c.to_bytes(4, 'little'); D = B_U - A_U
hi = load_hi(); rng = random.Random(5)
def single(b): return 0x21 <= b <= 0x7E or 0xA1 <= b <= 0xDF
def trailok(b): return 0x40 <= b <= 0x7E or 0x80 <= b <= 0xFC
typ = [x for x in hi if min(x[1]) >= 0.95 and textlike(b4(x[0])) >= 2]
adp = [x for x in typ if all(single(b) for b in b4(x[0])) and trailok(b4(x[0])[0])]
print('打てる塊 %d 個、うち 1 バイト文字 4 つで先頭が 2 バイト目にもなれるもの %d 個（値が 0xB0000000 以上: %d・%d）' % (len(typ), len(adp), sum(x[0] >= 0xB0000000 for x in typ), sum(x[0] >= 0xB0000000 for x in adp)))
def best_level(chunks, last):
    best = 0
    for p in itertools.permutations(chunks):
        t = textlike(b''.join(b4(c) for c in p) + last)
        if t > best: best = t
    return best
def sim(pool, n=6000, lasts=(b'Q', b'j')):
    ok = 0; m = 0
    while m < n:
        k1, k2 = rng.choice(pool)[0], rng.choice(pool)[0]; last = rng.choice(lasts)
        su = SIGMA - k1 - k2 - last[0]; x = (su - 3 * A_U) / D
        if k1 == k2 or not (1.55 <= x <= 1.95): continue
        c1 = rng.randrange(A_U, B_U); c2 = rng.randrange(A_U, B_U); c3 = su - c1 - c2
        if not (A_U <= c3 <= B_U): continue
        m += 1; ok += best_level([k1, k2, c1, c2, c3], last) >= 2
    return ok / n
print('打てる割合: 打てる塊ぜんぶ %.4f、1 バイト文字 4 つの塊だけ %.4f' % (sim(typ), sim(adp)))
