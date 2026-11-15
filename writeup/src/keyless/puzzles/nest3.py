"""nest3.py — 内側の連鎖の読み（nest.py）の追加の検査（鍵なし、2026-11-12）
  (1) 全 256 語の A⊕E・A⊕B の最上位ビットの並びと、一致語での A⊕E（= 平文の XOR そのもの）
  (2) 型 8 の行の組を外側のフラグで分けた集計（内側の連鎖が外側のフラグと独立に進むか）
  (3) ビット 5・6 の切り替わりが列どうしで独立か
  (4) フラグの読み（v2・alt・safe・cur）による違い
  python3 -I nest3.py"""
import sys, os, math, collections, itertools
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
from nest import W, X, read_flags, msb, pr
A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
same = [n for n in range(256) if W(A, n) == W(B, n)]
def rows_of(fl, lo, hi, col):
    """窓の行: (語, 外側のフラグ, PX, 型, 直前の外側の連鎖語)"""
    out = []; prevc = None
    for n in range(lo, hi + 1):
        if n not in fl: break
        if fl[n] == 1:
            px = X(W(A, n), W(E, n), W(A, prevc), W(E, prevc)) if prevc is not None else None; pc = prevc; prevc = n
        else: px = X(W(A, n), W(E, n)); pc = None
        if px is not None: out.append(dict(n=n, f=fl[n], px=px, t=8 if px[col] & 0x80 else 7, p=pc))
    return out
def st(ws):
    by = [b & 0x7f for w in ws for b in w]; n = len(by)
    if not n: return '—'
    return '%3d バイト: 0x20 未満 %3.0f%%、ビット 5 が 0 %3.0f%%、ビット 6 が 0 %3.0f%%' % (n, 100 * sum(v < 0x20 for v in by) / n, 100 * sum(not v & 0x20 for v in by) / n, 100 * sum(not v & 0x40 for v in by) / n)
def main():
    fl = read_flags(os.path.join(D, 'real_flags_v2.txt'))
    print('(1) 全 256 語: A⊕E の最上位ビット（X）、A⊕B（Y）、外側のフラグ（v2。1 = 連鎖、. = 連鎖なし、空白 = 不明）、* = A と B が一致')
    for r in range(0, 256, 32):
        print('  語 %3d〜: ' % r + ' '.join(msb(X(W(A, n), W(E, n))) for n in range(r, r + 8)) + ' …')
    print('  16 進で 1 語 1 文字（X の最上位ビット 4 個）:')
    for r in range(0, 256, 64):
        print('   %3d: X %s' % (r, ''.join('%x' % int(msb(X(W(A, n), W(E, n))), 2) for n in range(r, r + 64))))
        print('        Y %s' % ''.join('%x' % int(msb(X(W(A, n), W(B, n))), 2) for n in range(r, r + 64)))
        print('        f %s' % ''.join(('1' if fl[n] else '.') if n in fl else ' ' for n in range(r, r + 64)))
        print('        = %s' % ''.join('*' if n in same else ' ' for n in range(r, r + 64)))
    print('\n  一致語（%d 個）での A⊕E（どちらの読みでも、ここは平文の XOR そのもの）と、その前後の語' % len(same))
    for n in same:
        x = X(W(A, n), W(E, n)); print('   語 %3d（バイト %4d〜）: %s  %s  %-4s   A=%s E=%s' % (n, 4 * n, x.hex(), msb(x), pr(x), W(A, n).hex(), W(E, n).hex()))
    print('\n(2) 型 8 の行の組（その行 ⊕ 直前の型 8 の行）を、外側のフラグで分ける')
    for name in ('v2', 'alt', 'safe', 'cur'):
        fl = read_flags(os.path.join(D, 'real_flags_%s.txt' % name)); cls = collections.defaultdict(list); t7all = []; t8raw = []
        for wn, lo, hi, col in (('W2', 64, 119, 2), ('W1', 230, 255, 0)):
            R = rows_of(fl, lo, hi, col); t8 = [r for r in R if r['t'] == 8]; t7all += [r['px'] for r in R if r['t'] == 7]; t8raw += [r['px'] for r in t8]
            for a, b in zip(t8, t8[1:]):
                mx = X(a['px'], b['px'])
                if b['f'] == 0 and a['f'] == 0: k = 'どちらも外側が連鎖なし'
                elif b['f'] == 0: k = 'その行が外側で連鎖なし（MX = X_n ⊕ 前の行の PX）'
                elif a['f'] == 0: k = '前の行が外側で連鎖なし'
                elif b['p'] == a['n']: k = 'どちらも連鎖で、外側の連鎖でも隣（MX = X_n ⊕ X_{2 つ前の連鎖語}）'
                else: k = 'どちらも連鎖で、間に外側の連鎖語（型 7）をはさむ'
                cls[k].append(mx); cls['全部'].append(mx); cls['全部・' + wn].append(mx)
        print('  フラグ %s: 型 7 の PX %s' % (name, st(t7all))); print('      型 8 の PX（そのまま） %s' % st(t8raw))
        for k in sorted(cls): print('      %-46s %2d 組  %s' % (k, len(cls[k]), st(cls[k])))
    print('\n(3) 列どうしの独立（フラグ v2、型 8 の組ごとに、ビットが切り替わった列の数）')
    fl = read_flags(os.path.join(D, 'real_flags_v2.txt'))
    for bit, nm in ((0x20, 'ビット 5'), (0x40, 'ビット 6')):
        for wn, lo, hi, col in (('W2', 64, 119, 2), ('W1', 230, 255, 0)):
            R = rows_of(fl, lo, hi, col); t8 = [r for r in R if r['t'] == 8]; mxs = [X(a['px'], b['px']) for a, b in zip(t8, t8[1:])]
            cnt = collections.Counter(sum(1 for c in range(4) if w[c] & bit) for w in mxs); pc = [sum(1 for w in mxs if w[c] & bit) / len(mxs) for c in range(4)]
            # 独立なら: ポアソン二項
            exp = [0.0] * 5
            for pat in itertools.product((0, 1), repeat=4):
                p = 1.0
                for c in range(4): p *= pc[c] if pat[c] else 1 - pc[c]
                exp[sum(pat)] += p * len(mxs)
            t7 = [r['px'] for r in R if r['t'] == 7]; c7 = collections.Counter(sum(1 for c in range(4) if w[c] & bit) for w in t7)
            print('   %s・%s: 切り替わった列の数 0〜4 の組の数 %s（列が独立なら %s）。型 7 の PX でそのビットが 1 の列の数 %s' % (wn, nm, [cnt.get(k, 0) for k in range(5)], ['%.1f' % e for e in exp], [c7.get(k, 0) for k in range(5)]))
if __name__ == '__main__': main()
