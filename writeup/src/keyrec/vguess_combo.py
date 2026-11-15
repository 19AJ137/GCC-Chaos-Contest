# 鍵 v の候補: 鍵 k の語・社名・製品名・人名などを、GoodLuckTo… の言い回しや「・」でつないだ形（2026-11-14。該当なし）
import itertools
from collections import Counter
from vguess import score, TV
names = ['ChangxingYang','YangChangxing','GaoZhenyu','ZhenyuGao','Gao','Zhenyu','GAO','IISI','GCC','GCC6','GCC7','Chaos','ChaosInforGuard','InforGuard','Japan','China','You','Me','All','Everyone','Furukawa','Lijd','LiJD','Yang','Changxing','Dr.Gao','DrGao','GaoPhD','高振宇','国際情報科学研究所','カオス','磁石','発光','磁石・発光・IISI']
best = []; n = 0
def t(s):
    try: key = s.encode('cp932')
    except Exception: return
    if not (20 <= len(key) <= 48): return
    hit, r, d1, d3 = score(key, TV)
    best.append((hit, s))
    if hit >= 3: print('!!', hit, s)
for a, b in itertools.permutations(names, 2):
    for pat in ('GoodLuckTo%sGoodLuckTo%s', 'GoodLuckTo%sAndGoodLuckTo%s', 'GoodLuckTo%sAnd%s', 'GoodLuckTo%s&%s', 'GoodLuckTo%sGoodLuck%s', 'Good Luck To %s %s', 'GoodLuckTo%s,GoodLuckTo%s', 'GoodLuckTo%s・GoodLuckTo%s', 'GoodLuckTo%s%s'):
        for f in (str, str.lower, str.upper):
            t(f(pat % (a, b))); n += 1
for a, b, c in itertools.permutations(names, 3):
    for pat in ('GoodLuckTo%sGoodLuckTo%sGoodLuckTo%s', 'GoodLuckTo%s%s%s', '%s・%s・%s', '%s%s%s', '%s %s %s', '%s・%s・%s・IISI', 'GoodLuckTo%sGoodLuckTo%s%s'):
        t(pat % (a, b, c)); n += 1
print('made', n, '| tested (20-48 bytes)', len(best), '| matched items (of 5):', sorted(Counter(h for h, _ in best).items()))
