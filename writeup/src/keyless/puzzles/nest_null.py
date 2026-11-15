"""nest_null.py — nest.py の検査の帰無: 型 8 の行の PX を「直前の型 8 の行」ではなく、でたらめな相手（ほかの型 8 の行）と XOR したら、
  小さい値はどのくらい出るか（値の偏りだけで出る分）。直前の行との XOR だけが特別に小さければ、内側の連鎖の証拠になる。
  python3 -I nest_null.py [フラグ] [回数]"""
import sys, os, random, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
from nest import W, X, read_flags
def rows_of(A, E, fl, lo, hi, col):
    hi = max(n for n in range(lo, hi + 1) if n in fl and all(k in fl for k in range(lo, n + 1)))
    out = []; prevc = None
    for n in range(lo, hi + 1):
        if fl[n] == 1:
            px = X(W(A, n), W(E, n), W(A, prevc), W(E, prevc)) if prevc is not None else None; prevc = n
        else: px = X(W(A, n), W(E, n))
        if px is not None: out.append((n, 8 if px[col] & 0x80 else 7, px))
    return out
def stat(ws): return sum((b & 0x7f) < 0x20 for w in ws for b in w), sum((b & 0x7f) < 0x10 for w in ws for b in w)
def main():
    name = sys.argv[1] if len(sys.argv) > 1 else 'v2'; N = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
    fl = read_flags(os.path.join(D, 'real_flags_%s.txt' % name)); rnd = random.Random(3)
    tot_real = [0, 0]; tot_null = [[0, 0] for _ in range(N)]; nb = 0
    for wn, lo, hi, col in (('W2', 64, 119, 2), ('W1', 230, 255, 0)):
        R = rows_of(A, E, fl, lo, hi, col); t8 = [px for n, t, px in R if t == 8]; t7 = [px for n, t, px in R if t == 7]
        real = [X(t8[i], t8[i - 1]) for i in range(1, len(t8))]
        s20, s10 = stat(real); n = 4 * len(real)
        # 帰無 1: 相手をでたらめに選ぶ（自分と直前の行は除く）
        c20 = c10 = 0; m20 = m10 = 0
        for it in range(N):
            ws = []
            for i in range(1, len(t8)):
                while True:
                    j = rnd.randrange(len(t8))
                    if j != i and j != i - 1: break
                ws.append(X(t8[i], t8[j]))
            a, b = stat(ws); c20 += (a >= s20); c10 += (b >= s10); m20 += a; m10 += b
            tot_null[it][0] += a; tot_null[it][1] += b
        # 離れた相手（2 つ前・3 つ前の型 8 の行）
        lag = {k: stat([X(t8[i], t8[i - k]) for i in range(k, len(t8))]) + (4 * (len(t8) - k),) for k in (2, 3, 4)}
        a7 = stat(t7); a8 = stat(t8)
        print('%s: 型 8 の行 %d、型 7 の行 %d' % (wn, len(t8), len(t7)))
        print('   0x20 未満（0x10 未満）のバイト: 型 7 の PX %d（%d）/ %d = %.0f%%、型 8 の PX %d（%d）/ %d = %.0f%%' % (a7[0], a7[1], 4 * len(t7), 100 * a7[0] / (4 * len(t7)), a8[0], a8[1], 4 * len(t8), 100 * a8[0] / (4 * len(t8))))
        print('   型 8 の行 ⊕ 直前の型 8 の行: %d（%d）/ %d = %.0f%%。でたらめな相手との XOR: 平均 %.1f（%.1f）= %.0f%%、直前の行の値以上になる割合 %.5f（0x10 未満では %.5f）' % (s20, s10, n, 100 * s20 / n, m20 / N, m10 / N, 100 * m20 / N / n, c20 / N, c10 / N))
        print('   2・3・4 つ前の型 8 の行との XOR: ' + '、'.join('%d つ前 %d / %d = %.0f%%' % (k, v[0], v[2], 100 * v[0] / v[2]) for k, v in lag.items()))
        tot_real[0] += s20; tot_real[1] += s10; nb += n
    print('両方の窓: 直前の行との XOR で 0x20 未満 %d / %d = %.0f%%、でたらめな相手でこれ以上になる割合 %.5f（0x10 未満 %d、%.5f）' % (tot_real[0], nb, 100 * tot_real[0] / nb, sum(t[0] >= tot_real[0] for t in tot_null) / N, tot_real[1], sum(t[1] >= tot_real[1] for t in tot_null) / N))
main()
