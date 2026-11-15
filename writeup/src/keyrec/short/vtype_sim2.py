# 21 バイトで「打てる文字だけ」になる割合を、解く 3 塊の和（= 固定する 2 塊の和で決まる）ごとに見る
import sys, os, math, random, itertools, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from vpart import textlike
from vshort import load_hi, SIGMA, A_U, B_U
from vshort_count import irwin_hall
b4 = lambda c: c.to_bytes(4, 'little')
rng = random.Random(11)
def best_level(chunks, last):
    best = 0
    for p in itertools.permutations(chunks):
        t = textlike(b''.join(b4(c) for c in p) + last)
        if t > best: best = t
    return best
D = B_U - A_U
res = collections.defaultdict(lambda: [0, 0])
K = 0xC8C8C8C8                               # 固定する塊の代わり（打てる文字 4 つ）。和だけを変える
for su_x in [1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5]:
    su = int(3 * A_U + su_x * D); n = 0
    while n < 3000:
        c1 = rng.randrange(A_U, B_U); c2 = rng.randrange(A_U, B_U); c3 = su - c1 - c2
        if not (A_U <= c3 <= B_U): continue
        n += 1; res[su_x][0] += best_level([K, K, c1, c2, c3], b'Q') >= 2; res[su_x][1] += 1
    p = res[su_x][0] / res[su_x][1]; dens = irwin_hall(3, su_x)
    print('解く 3 塊の和の位置 x = %.1f（平均の 4 バイト目 0x%02x、固定 2 塊の和 %.2fe9）: 打てる割合 %.4f、密度 %.3f、積 %.4f' % (su_x, (su // 3) >> 24, (SIGMA - su - 0x51) / 1e9, p, dens, p * dens))
