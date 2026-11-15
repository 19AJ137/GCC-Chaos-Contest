"""nest2.py — nest.py の続き: 内側の連鎖の検査を、ビットごと・相手の選び方ごと・文章らしさの尤度で数える（鍵なし、2026-11-12）
  PX = 外側（鍵 k）の連鎖を外した A・E の平文の XOR。型 8 = 8bit の列の最上位ビットが 1 の行、型 7 = 0 の行。
  相手の選び方: (P1) 直前の型 8 の行（型 7 の行は飛ばす = 選んだ語だけ連鎖する GCC6 の形）、(P2) 直前の行（型によらない）、(P0) でたらめな型 8 の行（帰無）
  python3 -I nest2.py [フラグ（既定 v2）] [帰無の回数]"""
import sys, os, random, math, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
from nest import W, X, read_flags
TEXT = ("the international information science institute announces a prize contest for the gcc chaos cryptosystem . "
        "if you can decrypt the following ciphertext and find the key , please send the key and the plaintext to us by mail or fax . "
        "the prize is ten million yen and the deadline is december thirty first nineteen ninety seven . this software is distributed "
        "on cd rom with several magazines in japan . please read the readme file before you install the program on your computer . "
        "our address is tokyo japan , tel 03 3456 7890 , fax 03 3456 7891 , e-mail info at iisi co jp . thank you for your interest "
        "in our products . the security of information is very important for the network society of the next century . "
        "we believe that chaos theory gives a new method for cryptography which is fast and strong . copyright 1996 all rights reserved . ")
def text_xor_dist(rnd, n=200000):
    """別々の英文 2 本の同じ位置のバイトの XOR の分布（小文字・空白・数字が中心の文）"""
    b = TEXT.encode(); c = collections.Counter()
    for _ in range(n): c[b[rnd.randrange(len(b))] ^ b[rnd.randrange(len(b))]] += 1
    return {v: (c.get(v, 0) + 0.5) / (n + 64) for v in range(128)}
def rows_of(A, E, fl, lo, hi, col):
    hi = max(n for n in range(lo, hi + 1) if n in fl and all(k in fl for k in range(lo, n + 1)))
    out = []; prevc = None
    for n in range(lo, hi + 1):
        if fl[n] == 1:
            px = X(W(A, n), W(E, n), W(A, prevc), W(E, prevc)) if prevc is not None else None; prevc = n
        else: px = X(W(A, n), W(E, n))
        if px is not None: out.append((n, 8 if px[col] & 0x80 else 7, px, fl[n]))
    return out
def stats(ws, PT):
    by = [b & 0x7f for w in ws for b in w]; n = len(by)
    return dict(n=n, lt20=sum(v < 0x20 for v in by), b5=sum((v & 0x20) == 0 for v in by), b6=sum((v & 0x40) == 0 for v in by),
                llr=sum(math.log2(PT[v] * 128) for v in by))
def fmt(s): return '%3d バイト: 0x20 未満 %3.0f%%、ビット 5 が 0 %3.0f%%、ビット 6 が 0 %3.0f%%、文章の XOR らしさ %+6.1f ビット（1 バイトあたり %+.2f）' % (s['n'], 100 * s['lt20'] / s['n'], 100 * s['b5'] / s['n'], 100 * s['b6'] / s['n'], s['llr'], s['llr'] / s['n'])
def main():
    name = sys.argv[1] if len(sys.argv) > 1 else 'v2'; N = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
    A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
    fl = read_flags(os.path.join(D, 'real_flags_%s.txt' % name)); rnd = random.Random(11); PT = text_xor_dist(rnd)
    print('参考: 英文どうしの XOR の分布（見本の文から）: 0x20 未満 %.0f%%、ビット 5 が 0 %.0f%%、ビット 6 が 0 %.0f%%' % (100 * sum(PT[v] for v in range(0x20)), 100 * sum(PT[v] for v in range(128) if not v & 0x20), 100 * sum(PT[v] for v in range(128) if not v & 0x40)))
    agg = collections.defaultdict(list); nullagg = [collections.defaultdict(float) for _ in range(N)]
    for wn, lo, hi, col in (('W2', 64, 119, 2), ('W1', 230, 255, 0)):
        R = rows_of(A, E, fl, lo, hi, col); t8i = [i for i, r in enumerate(R) if r[1] == 8]; t7 = [r[2] for r in R if r[1] == 7]; t8 = [R[i][2] for i in t8i]
        p1 = [X(t8[i], t8[i - 1]) for i in range(1, len(t8))]                       # 直前の型 8 の行
        p2 = [X(R[i][2], R[i - 1][2]) for i in t8i if i >= 1]                      # 直前の行（型によらない）
        p2b = [X(R[i][2], R[i - 1][2]) for i in t8i if i >= 1 and R[i - 1][1] == 7]   # そのうち、直前が型 7 の行だけ
        p1b = [X(t8[k], t8[k - 1]) for k in range(1, len(t8)) if R[t8i[k] - 1][1] == 7]  # 同じ行を、型 7 を飛ばして型 8 の行と組にしたもの
        print('\n%s（型 7 %d 行、型 8 %d 行）' % (wn, len(t7), len(t8)))
        for lab, ws in (('型 7 の PX', t7), ('型 8 の PX（そのまま）', t8), ('型 8 ⊕ 直前の型 8（P1）', p1), ('型 8 ⊕ 直前の行（P2）', p2), ('  うち直前が型 7 の行: P2', p2b), ('  同じ行を P1 で', p1b)):
            if ws: print('   %-28s %s' % (lab, fmt(stats(ws, PT))))
        s1 = stats(p1, PT); cnt = collections.Counter(); mean = collections.Counter()
        for it in range(N):
            ws = []
            for i in range(1, len(t8)):
                while True:
                    j = rnd.randrange(len(t8))
                    if j != i and j != i - 1: break
                ws.append(X(t8[i], t8[j]))
            s0 = stats(ws, PT)
            for k in ('lt20', 'b5', 'llr'):
                cnt[k] += (s0[k] >= s1[k]); mean[k] += s0[k]; nullagg[it][k] += s0[k]
        print('   でたらめな相手（P0、%d 回）: 0x20 未満 平均 %.0f%%、ビット 5 が 0 平均 %.0f%%、文章らしさ 平均 %+.1f。P1 の値以上になる割合: 0x20 未満 %.5f、ビット 5 %.5f、文章らしさ %.5f' % (
            N, 100 * mean['lt20'] / N / s1['n'], 100 * mean['b5'] / N / s1['n'], mean['llr'] / N, cnt['lt20'] / N, cnt['b5'] / N, cnt['llr'] / N))
        for k in ('lt20', 'b5', 'llr'): agg[k].append(s1[k])
        # 列ごと・ビットごと: 型 8 の行を順に並べたとき、PX のビットが切り替わる割合（P1 の MX のビットが 1 の割合）と、型 7 の PX でそのビットが 1 の割合
        print('   ビットごと（上位から 7・6・5・4・3・2・1・0）: P1 の MX で 1 の割合 ／ 型 7 の PX で 1 の割合 ／ 型 8 の PX（そのまま）で 1 の割合')
        for c in range(4):
            a = ['%2.0f' % (100 * sum((w[c] >> b) & 1 for w in p1) / len(p1)) for b in range(7, -1, -1)]
            b7 = ['%2.0f' % (100 * sum((w[c] >> b) & 1 for w in t7) / len(t7)) for b in range(7, -1, -1)]
            r8 = ['%2.0f' % (100 * sum((w[c] >> b) & 1 for w in t8) / len(t8)) for b in range(7, -1, -1)]
            print('      位置 %d: %s ／ %s ／ %s' % (c, ' '.join(a), ' '.join(b7), ' '.join(r8)))
    tot = {k: sum(agg[k]) for k in agg}
    print('\n両方の窓の合計: P1 の文章らしさ %+.1f ビット。でたらめな相手でこれ以上になる割合 %.5f（0x20 未満 %.5f、ビット 5 %.5f）' % (tot['llr'], sum(d['llr'] >= tot['llr'] for d in nullagg) / N, sum(d['lt20'] >= tot['lt20'] for d in nullagg) / N, sum(d['b5'] >= tot['b5'] for d in nullagg) / N))
    import statistics
    v = [d['llr'] for d in nullagg]; print('   帰無の文章らしさ: 平均 %+.1f、標準偏差 %.1f → P1 は %.1f σ' % (statistics.mean(v), statistics.pstdev(v), (tot['llr'] - statistics.mean(v)) / statistics.pstdev(v)))
if __name__ == '__main__': main()
