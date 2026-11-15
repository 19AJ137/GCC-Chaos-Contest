# 「打てる文字だけ」（vpart.textlike が 2 以上）になる割合の見積もり。固定する塊は打てる文字のものを選び、解く 3 塊は乱数で代用する。
import sys, os, math, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from vpart import textlike
from vshort import load_hi, SIGMA, A_U, B_U
hi = load_hi()
b4 = lambda c: c.to_bytes(4, 'little')
lv = {}
for c, f, m in hi: lv[c] = textlike(b4(c))
print('小数部が 4 つとも 0.90 以上の塊 %d 個のうち、4 バイトだけで打てる文字（区分 2 以上）: %d 個、区分 3: %d 個' % (len(hi), sum(v >= 2 for v in lv.values()), sum(v >= 3 for v in lv.values())))
for thr in (0.95, 0.97, 0.98):
    s = [c for c, f, m in hi if min(f) >= thr]; print('  %.2f 以上 %d 個のうち 区分 2 以上 %d、区分 3 %d、写像 1 かつ区分 2 以上 %d' % (thr, len(s), sum(lv[c] >= 2 for c in s), sum(lv[c] >= 3 for c in s), sum(1 for c, f, m in hi if min(f) >= thr and m == 1 and lv[c] >= 2)))
def best_level(chunks, last):
    best = 0
    for p in itertools.permutations(chunks):
        t = textlike(b''.join(b4(c) for c in p) + last)
        if t > best: best = t
        if best == 3: break
    return best
rng = random.Random(7)
typ = [x for x in hi[:6000] if lv[x[0]] >= 2]
def sim(L, n=20000):
    nfull = L // 4; nb = L % 4; nk = nfull - 3; cnt = [0, 0, 0, 0]; closed = 0
    last = {0: None, 1: b'Q', 2: bytes.fromhex('4b46'), 3: bytes.fromhex('479b80')}[nb]
    for _ in range(n):
        while True:
            ks = [rng.choice(typ)[0] for _ in range(nk)]
            if nb == 0:
                lastk = rng.choice([x for x in typ[:2000] if x[2] == 1])[0]; lastb = b4(lastk); ks = ks[:-1]
            else: lastb = last; lastk = int.from_bytes(last, 'little')
            su = SIGMA - sum(ks) - lastk
            c1 = rng.randrange(A_U, B_U); c2 = rng.randrange(A_U, B_U); c3 = su - c1 - c2
            if A_U <= c3 <= B_U: break
        us = [c1, c2, c3]
        cnt[best_level(ks + us, lastb)] += 1
        closed += all(textlike(b4(c)) >= 2 for c in us)
    return [x / n for x in cnt], closed / n
for L in (21, 22, 24, 25, 28):
    r, cl = sim(L, 4000 if L < 28 else 1500); print('長さ %d: 並べ替えのどれかで 区分 1 以上 %.3f、2 以上 %.4f、3 %.5f（解く 3 塊がそれぞれ単独で打てる %.4f）' % (L, r[1] + r[2] + r[3], r[2] + r[3], r[3], cl))
