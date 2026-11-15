"""puzzles.py — 「A は暗号化、B は復号」の推測で残る不可解な点の、鍵なしの数の検算（2026-11-12）
  python3 -I puzzles.py [フラグ: v2（読み α。既定）・alt（読み β）・cur（訂正前）・safe]
  A・E は通常の暗号化（連鎖は暗号文側）、B は復号 1 回の形（連鎖は平文側）として、暗号文だけから出る量を数える。
    X = A⊕E、Y = A⊕B。連鎖なしの語では X = PA⊕PE、Y = PA⊕PB。連鎖ありの語（直前の連鎖語 p）では PX = X_n⊕X_p（A・E の平文の XOR）。
    窓: W2 = 語 64〜119（8bit の列は位置 2）、W1 = 語 224〜255（8bit の列は位置 0）"""
import sys, os, math, itertools, collections
D = os.path.dirname(os.path.abspath(__file__))
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
def W(d, i): return d[4 * i:4 * i + 4]
def X(*a):
    r = bytes(4)
    for x in a: r = bytes(p ^ q for p, q in zip(r, x))
    return r
def msb(b): return ''.join('1' if c & 0x80 else '0' for c in b)
def read_flags(fn):
    fl = {}
    for line in open(fn):
        line = line.strip()
        if not line: continue
        tag, body = line.split(':', 1)
        if tag == 'z':
            for w in body.split(','):
                if w.strip(): fl[int(w)] = 0
        else:
            b = int(tag)
            for k, c in enumerate(body.strip()): fl[b + k] = 1 if c == '1' else 0
    return fl
def prev_chain(fl, i):
    j = i - 1
    while j >= 0:
        if j not in fl: return None
        if fl[j] == 1: return j
        j -= 1
    return None
def binom_le(k, n, p): return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(0, k + 1))
def binom_ge(k, n, p): return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1))
def fisher_ge(a, b, c, d):   # 2×2 表 [[a,b],[c,d]] で a がこれ以上になる確率（片側、超幾何）
    n1, n2, k = a + b, c + d, a + c
    tot = math.comb(n1 + n2, k); return sum(math.comb(n1, x) * math.comb(n2, k - x) for x in range(a, min(n1, k) + 1)) / tot

name = sys.argv[1] if len(sys.argv) > 1 else 'v2'
A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
fl = read_flags(os.path.join(D, 'real_flags_%s.txt' % name))
WIN = {'W2': (64, 119, 2), 'W1': (224, 255, 0)}
def win_of(n):
    for k, (lo, hi, col) in WIN.items():
        if lo <= n <= hi: return k
    return None
print('=== フラグ %s（既知 %d 語: 連鎖あり %d、連鎖なし %d） ===' % (name, len(fl), sum(v == 1 for v in fl.values()), sum(v == 0 for v in fl.values())))

# 1. 一致語と、バイト単位の一致
same = [n for n in range(256) if W(A, n) == W(B, n)]
print('\n[1] A = B の語（一致語） %d 個: %s' % (len(same), same))
eq = collections.Counter(sum(a == b for a, b in zip(W(A, n), W(B, n))) for n in range(256))
print('    語ごとの「A と B で同じバイトの数」: %s （一致語を除く 236 語で同じバイト %d 個、偶然なら %.1f）' % (dict(sorted(eq.items())), sum(k * v for k, v in eq.items() if k < 4), 236 * 4 / 256))
print('    一致語の隣の語（±1）で同じバイト: %d' % sum(sum(a == b for a, b in zip(W(A, m), W(B, m))) for n in same for m in (n - 1, n + 1) if 0 <= m < 256 and m not in same))

# 2. 窓の連鎖なしの語: 共有と行の型
print('\n[2] 窓の連鎖なしの語（X = PA⊕PE、Y = PA⊕PB がそのまま出る）')
rows = []
for n in range(256):
    w = win_of(n)
    if w is None or fl.get(n) != 0: continue
    col = WIN[w][2]; x = X(W(A, n), W(E, n)); y = X(W(A, n), W(B, n)); v = X(W(B, n), W(E, n))
    rows.append((n, w, n in same, (x[col] >> 7) & 1, x, y, v))
    print('    語 %3d %s %s  X %s (%s) 8bit の列 %d  Y %s (%s)  B⊕E %s (%s)' % (n, w, '共有  ' if n in same else '非共有', x.hex(), msb(x), (x[col] >> 7) & 1, y.hex(), msb(y), v.hex(), msb(v)))
a = sum(1 for r in rows if r[2] and r[3] == 0); b = sum(1 for r in rows if not r[2] and r[3] == 0); c = sum(1 for r in rows if r[2] and r[3] == 1); d = sum(1 for r in rows if not r[2] and r[3] == 1)
print('    行の型と共有: 8bit の列が 7bit の行 %d（共有 %d・非共有 %d）、8bit の行 %d（共有 %d・非共有 %d）。規則に合う %d / %d、Fisher 片側 P = %.2g' % (a + b, a, b, c + d, c, d, a + d, len(rows), fisher_ge(a, b, c, d)))
ns = [r for r in rows if not r[2]]
ym = [sum((r[5][k] >> 7) & 1 for r in ns) for k in range(4)]; yb = [bb for r in ns for bb in r[5]]
print('    非共有の %d 語: Y の最上位ビットが 1 のバイト 位置別 %s（計 %d / %d = %.0f%%）、Y が 0 のバイト %d 個（偶然なら %.2f）、Y < 0x20 のバイト %d 個（偶然なら %.1f）' % (
    len(ns), ym, sum(ym), len(yb), 100 * sum(ym) / len(yb), sum(v == 0 for v in yb), len(yb) / 256, sum(v < 0x20 for v in yb), len(yb) / 8))
for k, (lo, hi, col) in WIN.items():
    J = [c for c in range(4) if c != col]; r2 = [r for r in ns if r[1] == k]
    print('      %s: 非共有 %d 語で、一致列（%s）の Y の最上位ビットが 1 のバイト %d / %d、8bit の列（位置 %d）で %d / %d' % (k, len(r2), J, sum((r[5][c] >> 7) & 1 for r in r2 for c in J), 3 * len(r2), col, sum((r[5][col] >> 7) & 1 for r in r2), len(r2)))

# 3. 窓の連鎖ありの語
print('\n[3] 窓の連鎖ありの語（直前の連鎖語 p が窓の中で決まるもの）')
ch = []
for n in range(256):
    w = win_of(n)
    if w is None or fl.get(n) != 1: continue
    p = prev_chain(fl, n)
    if p is None or win_of(p) != w: continue
    col = WIN[w][2]; px = X(W(A, n), W(E, n), W(A, p), W(E, p)); y = X(W(A, n), W(B, n)); q = X(y, W(A, p))
    ch.append((n, w, p, px, y, q, (px[col] >> 7) & 1))
for k, (lo, hi, col) in WIN.items():
    c2 = [r for r in ch if r[1] == k]; J = [c for c in range(4) if c != col]
    pred = [r[0] for r in c2 if r[6] == 0]
    print('  %s: 連鎖語 %d。8bit の列が 7bit の行（規則なら共有） %d 語: %s' % (k, len(c2), len(pred), pred))
    # 一致列の PX（A・E の平文の XOR）の最上位ビット
    bad = [r[0] for r in c2 if any((r[3][c] >> 7) & 1 for c in J)]
    print('      A・E の平文の XOR が一致列で最上位ビット 1 になる語: %s' % bad)
    # (a) 暗号文の最上位ビット: Y = A⊕B
    ymc = [sum((r[4][c] >> 7) & 1 for r in c2) for c in range(4)]
    print('      Y = A⊕B の最上位ビットが 1 の語数（位置 0〜3）: %s / %d  両側 P（二項 1/2）: %s' % (ymc, len(c2), ['%.3f' % min(1, 2 * min(binom_le(v, len(c2), .5), binom_ge(v, len(c2), .5))) for v in ymc]))
    # (b) D = Y⊕A_p（復号 1 回の形なら PA_n⊕PB_n⊕PB_p）: 一致列の型
    pat = collections.Counter(''.join(str((r[5][c] >> 7) & 1) for c in J) for r in c2)
    n000 = pat.get('000', 0)
    print('      D = Y⊕A_p の一致列の最上位ビットの型: %s  000 は %d / %d（偶然なら %.1f、これ以下になる確率 %.4f）' % (dict(sorted(pat.items(), key=lambda kv: -kv[1])), n000, len(c2), len(c2) / 8, binom_le(n000, len(c2), 1 / 8)))
    dm = [sum((r[5][c] >> 7) & 1 for r in c2) for c in range(4)]
    print('      D の最上位ビットが 1 の語数（位置 0〜3）: %s / %d' % (dm, len(c2)))
    # (c) 厳密な痕跡: Y_n = Y_p（B も暗号文側で連鎖なら、共有の語で成り立つ）、Y_n = A_p（B が連鎖なしなら）
    print('      厳密な痕跡: Y_n = Y_p の語 %d、Y_n = A_p の語 %d' % (sum(1 for r in c2 if r[4] == X(W(A, r[2]), W(B, r[2]))), sum(1 for r in c2 if r[4] == W(A, r[2]))))
    # (d) 規則で共有と予想される行で、B に足された値 R = B_n⊕A_n⊕A_p（共有なら B のレジスタそのもの）
    for r in c2:
        if r[6] == 0: print('        語 %3d（p = %3d）: 共有なら B に足された値 R = %s（最上位ビット %s）' % (r[0], r[2], r[5].hex(), msb(r[5])))
allq = sum(1 for r in ch if all(((r[5][c] >> 7) & 1) == 0 for c in range(4) if c != WIN[r[1]][2]))
print('  両方の窓: 連鎖語 %d、D の一致列が 000 の語 %d（偶然なら %.1f）' % (len(ch), allq, len(ch) / 8))
# 全語の走査（フラグによらない）
ys = {n: X(W(A, n), W(B, n)) for n in range(256)}
print('  フラグによらない走査: Y_n が、それより前の A・B・E の語か Y の語（0 を除く）と一致する数: A %d、B %d、E %d、Y %d' % (
    sum(1 for n in range(256) for m in range(n) if ys[n] == W(A, m) and any(ys[n])), sum(1 for n in range(256) for m in range(n) if ys[n] == W(B, m) and any(ys[n])),
    sum(1 for n in range(256) for m in range(n) if ys[n] == W(E, m) and any(ys[n])), sum(1 for n in range(256) for m in range(n) if ys[n] == ys[m] and any(ys[n]))))

# 4. 窓の端（語 116〜121、224〜230）
print('\n[4] 窓の端の語（X の最上位ビット。連鎖ありなら PX = X_n⊕X_p）')
for n in list(range(114, 123)) + list(range(222, 232)):
    f = fl.get(n); p = prev_chain(fl, n) if f == 1 else None
    px = X(W(A, n), W(E, n), W(A, p), W(E, p)) if (f == 1 and p is not None) else (X(W(A, n), W(E, n)) if f == 0 else None)
    print('    語 %3d フラグ %s %s  X %s (%s)  A・E の平文の XOR %s' % (n, {1: '連鎖あり', 0: '連鎖なし', None: '不明    '}[f], 'A=B' if n in same else '   ', X(W(A, n), W(E, n)).hex(), msb(X(W(A, n), W(E, n))), ('%s (%s)' % (px.hex(), msb(px))) if px else '—'))
