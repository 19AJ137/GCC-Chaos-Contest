"""facts.py — 復元した結果の要約（レポート用の数字）: 連鎖フラグ、これまでの読みとの照合、平文の区切り、過去の不可解な点の答え合わせ用の値"""
import sys, os, pickle, hashlib
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D); sys.path.insert(0, os.path.join(D, '..', 'puzzles'))
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
from nest import read_flags
ks = pickle.load(open(os.path.join(D, 'keystreams.pkl'), 'rb')); F = ks['F']; FP = ks['FP']
A = open(os.path.join(CT, 'A.bin'), 'rb').read(); B = open(os.path.join(CT, 'B.bin'), 'rb').read(); E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
ma = open(os.path.join(D, 'plain_A.bin'), 'rb').read(); me = open(os.path.join(D, 'plain_E_1024.bin'), 'rb').read()
W = lambda d, n: d[4 * n:4 * n + 4]
same = [n for n in range(256) if W(A, n) == W(B, n)]
print('外側の鍵 k のフラグ（1 = 連鎖、. = 連鎖なし）:'); 
for r in range(0, 256, 64): print('  %3d: %s' % (r, ''.join('1' if F[n] else '.' for n in range(r, r + 64))))
print('内側の鍵 v のフラグ:')
for r in range(0, 256, 64): print('  %3d: %s' % (r, ''.join('1' if FP[n] else '.' for n in range(r, r + 64))))
print('連鎖なしの語: k %d / 256（%.1f%%）、v %d / 256（%.1f%%）、両方 %d' % (F.count(0), 100 * F.count(0) / 256, FP.count(0), 100 * FP.count(0) / 256, sum(1 for n in range(256) if not F[n] and not FP[n])))
both = [n for n in range(256) if not F[n] and not FP[n]]
print('両方とも連鎖なしの語: %s' % both); print('A = B の語: %s' % same); print('A = B で両方連鎖なしでない語: %s、両方連鎖なしで A ≠ B の語: %s' % ([n for n in same if n not in both], [n for n in both if n not in same]))
print('語 0: k %d、v %d（最初の語は、連鎖があっても IV が同じなので A = B になる）' % (F[0], FP[0]))
for name in ('v2', 'alt', 'safe', 'cur'):
    fl = read_flags(os.path.join(D, '..', 'puzzles', 'real_flags_%s.txt' % name)); bad = [n for n in sorted(fl) if fl[n] != F[n]]
    print('これまでのフラグの読み %s: %d 語のうち外側の鍵のフラグと違う語 %s' % (name, len(fl), bad))
for n in (117, 118, 119, 120, 224, 225, 226, 227, 228, 229, 230):
    print('  語 %d: k %d、v %d、A の平文 %s、E の平文 %s%s' % (n, F[n], FP[n], W(ma, n).hex(), W(me, n).hex(), '  A=B' if n in same else ''))
# 8bit のバイトの位置
def runs(pred):
    out = []; s = None
    for i in range(1025):
        ok = i < 1024 and pred(i)
        if ok and s is None: s = i
        if not ok and s is not None: out.append((s, i - 1)); s = None
    return out
print('A の平文で 8bit のバイトがある範囲（8 バイト以内の切れ目はつなぐ）:')
def merged(rs, gap=8):
    out = []
    for a, b in rs:
        if out and a - out[-1][1] <= gap: out[-1] = (out[-1][0], b)
        else: out.append((a, b))
    return out
print('  A: %s' % merged(runs(lambda i: ma[i] >= 0x80))); print('  E: %s' % merged(runs(lambda i: me[i] >= 0x80)))
print('両方とも 7bit のバイトが続く範囲（16 バイト以上）: %s' % [(a, b) for a, b in runs(lambda i: ma[i] < 0x80 and me[i] < 0x80) if b - a >= 15])
print('A の平文 SHA-256 %s' % hashlib.sha256(ma).hexdigest()); print('E の平文（先頭 1,024 B）SHA-256 %s' % hashlib.sha256(me).hexdigest())
# 文字数など
ta = ma.decode('cp932'); te = me.decode('cp932', 'replace')
print('A の平文: %d バイト、%d 文字、行数 %d' % (len(ma), len(ta), ta.count('\n'))); print('E の平文（先頭）: %d バイト' % len(me))
# 2 回出てくる文
s1 = b"Is true? Do you think you can decrypte GCC? No, you cann't! Because this is GCC6 only!\r\nYou can not decrypte GCC7!\r\n"
print('A の冒頭の英文（%d バイト）の出現位置: %s' % (len(s1), [i for i in range(1024) if ma.startswith(s1, i)]))
s2 = 'GCCカオス暗号は、共通キー型のデジタル型カオス暗号であり、'.encode('cp932')
print('「GCCカオス暗号は、共通キー型の…」の位置: A %s、E %s' % ([i for i in range(1024) if ma.startswith(s2, i)], [i for i in range(1024) if me.startswith(s2, i)]))
s3 = b'GCC is a digital chaos cryptosystem; it is conventional cryptosystem '
print('「GCC is a digital chaos cryptosystem; …」の位置: E %s' % [i for i in range(1024) if me.startswith(s3, i)])
