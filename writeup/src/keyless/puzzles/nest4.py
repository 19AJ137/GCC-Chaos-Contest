"""nest4.py — 2 本の連鎖の読み分け（鍵なし、2026-11-12）
  窓の中の A⊕E（= X）の最上位ビットは 2 種類の列にしか立たない（W2: 位置 2 と 3、W1: 位置 0 と 1〜3）。
  2 重の連鎖 (1+S_内)(1+S_外) X = MX（MX は 7bit）なら、
     X_n の最上位ビット = f外_n・ρ_n ⊕ f内_n・σ      （ρ_n = 直前の外側の連鎖語の X の最上位ビット、σ = 内側のレジスタの最上位ビット（一定））
  で、「外側のフラグがそのまま見える列」と「内側のフラグが外側の連鎖語ごとに積み上がる列」が 1 つずつある。どちらの列がどちらかは、最上位ビットだけでは決まらない:
     読み A: 外側 = W2 位置 3・W1 位置 1〜3（これまでの連鎖フラグ）、内側 = W2 位置 2・W1 位置 0（これまでの「8bit の列」）
     読み B: 外側 = W2 位置 2・W1 位置 0、内側 = W2 位置 3・W1 位置 1〜3
  両方で内側まで外し、残り（MX）が文章どうしの XOR らしいかを比べる。フラグのファイルは使わず、X の最上位ビットだけから決める。
  python3 -I nest4.py"""
import sys, os, math, collections, random
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
from nest import W, X, msb, pr
A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
same = set(n for n in range(256) if W(A, n) == W(B, n))
XW = [X(W(A, n), W(E, n)) for n in range(256)]
def decode(lo, hi, oc, ic):
    """語 lo〜hi。外側のフラグ = 列 oc の最上位ビット（X そのまま）。内側のフラグ = 外側を外した PX の列 ic の最上位ビット。
       返り値: 行のリスト dict(n, fo, fi, px, mx)。px・mx は決まらないとき None"""
    rows = []; prevo = None; previ = None
    for n in range(lo, hi + 1):
        x = XW[n]; fo = 1 if x[oc] & 0x80 else 0
        if fo:
            px = X(x, XW[prevo]) if prevo is not None else None; prevo = n
        else: px = x
        fi = None; mx = None
        if px is not None:
            fi = 1 if px[ic] & 0x80 else 0
            if fi:
                mx = X(px, previ) if previ is not None else None; previ = px
            else: mx = px
        rows.append(dict(n=n, fo=fo, fi=fi, px=px, mx=mx, x=x))
    return rows
def stat(ws):
    by = [b for w in ws for b in w]; n = len(by)
    if not n: return None
    return dict(n=n, hi=sum(1 for v in by if v & 0x80), lt20=sum((v & 0x7f) < 0x20 for v in by), b5=sum(not v & 0x20 for v in by), b6=sum(not v & 0x40 for v in by))
def fmt(s):
    if not s: return '—'
    return '%3d バイト: 最上位ビットが 1 %2d、0x20 未満 %3.0f%%、ビット 5 が 0 %3.0f%%、ビット 6 が 0 %3.0f%%' % (s['n'], s['hi'], 100 * s['lt20'] / s['n'], 100 * s['b5'] / s['n'], 100 * s['b6'] / s['n'])
def report(name, lo, hi, oc, ic, show=False):
    R = decode(lo, hi, oc, ic)
    cls = collections.defaultdict(list)
    for r in R:
        if r['mx'] is None: continue
        cls[(r['fo'], r['fi'])].append(r['mx']); cls['all'].append(r['mx'])
        cls['fi%d' % r['fi']].append(r['mx'])
    print('  %s（語 %d〜%d、外側 = 位置 %d、内側 = 位置 %d）: 外側の連鎖なし %d / %d、内側の連鎖なし %d / %d' % (name, lo, hi, oc, ic, sum(1 for r in R if r['fo'] == 0), len(R), sum(1 for r in R if r['fi'] == 0), sum(1 for r in R if r['fi'] is not None)))
    for k, lab in (('all', 'MX 全部'), ('fi0', '内側が連鎖なし（MX = PX）'), ('fi1', '内側が連鎖（MX = PX ⊕ 直前の内側の連鎖語の PX）'), ((0, 0), '  外 0・内 0'), ((1, 0), '  外 1・内 0'), ((0, 1), '  外 0・内 1'), ((1, 1), '  外 1・内 1')):
        print('     %-34s %s' % (lab, fmt(stat(cls[k]))))
    for cols, lab in (((0, 1), '位置 0・1'), ((2, 3), '位置 2・3')):
        for k, kl in (('fi0', '内 0'), ('fi1', '内 1')):
            by = [w[c] for w in cls[k] for c in cols]
            if by: print('     %s・%s: ビット 6 が 1 %2d / %2d、ビット 5 が 1 %2d / %2d' % (lab, kl, sum(1 for v in by if v & 0x40), len(by), sum(1 for v in by if v & 0x20), len(by)))
    nm = [r['n'] for r in R if r['fo'] == 0 and r['fi'] == 0]; print('     外 0・内 0 の語: %s。そのうち A = B: %s。A = B なのに外 0・内 0 でない語: %s' % (nm, [n for n in nm if n in same], [r['n'] for r in R if r['n'] in same and not (r['fo'] == 0 and r['fi'] == 0)]))
    if show:
        print('      語  X        外 内  PX                MX')
        for r in R:
            print('     %3d  %s %s  %d  %s  %-24s %s%s' % (r['n'], r['x'].hex(), msb(r['x']), r['fo'], '-' if r['fi'] is None else r['fi'], (r['px'].hex() + ' ' + msb(r['px'])) if r['px'] else '—', (r['mx'].hex() + ' ' + pr(r['mx'])) if r['mx'] else '—', '   A=B' if r['n'] in same else ''))
    return R
def main():
    show = len(sys.argv) > 1
    print('読み A（外側 = これまでの連鎖フラグの列、内側 = これまでの 8bit の列）')
    report('W2', 55, 118, 3, 2, show); report('W1', 229, 255, 1, 0, show)
    print('\n読み B（外側 = これまでの 8bit の列、内側 = これまでの連鎖フラグの列）')
    report('W2', 55, 118, 2, 3, show); report('W1', 229, 255, 0, 1, show)
if __name__ == '__main__': main()
