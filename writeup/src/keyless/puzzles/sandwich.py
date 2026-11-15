"""sandwich.py — 「懸賞の鍵 k の処理を、別の鍵 v の処理 2 回で挟んだ形」（B = V(K(V(P)))）の検査（鍵なし。v は候補の鍵、2026-11-12）
  この形では、v の鍵ストリームは v の連鎖なしの語でだけ打ち消し合う。したがって A = B になれるのは「k でも v でも連鎖なし」の語だけで、
  v が分かっていれば、一致語（語 0 = ファイルの最初の連鎖語を除く）がすべて v の連鎖なしの語でなければならない。
  python3 -I sandwich.py                      既定の金庫の鍵「123」などの候補で、一致語と v の連鎖フラグを突き合わせる
  python3 -I sandwich.py null [本数]          乱数の鍵で、一致語 19 語のうち v の連鎖なしに入る数の分布（帰無）"""
import sys, os, random, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D, '..', '..', 'gcc6'))
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
import gcc
def W(d, i): return d[4 * i:4 * i + 4]
A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]
same = [n for n in range(256) if W(A, n) == W(B, n)]
NONSH = [67, 77, 85, 97, 103, 109, 117, 237, 251]          # 窓の連鎖なし（k）で一致しない語（読み α）
def flags(key):
    st = gcc.key_schedule(key); f = []
    for n in range(256):
        x = gcc.step(st); k, nc = gcc.ks_word(x); f.append(0 if nc else 1)
    return st.mode, f
def report(name, key):
    try: mode, f = flags(key)
    except ValueError as e: print('%-28s 鍵として使えない（%s）' % (name, e)); return
    nc = [n for n in range(256) if f[n] == 0]; first = f.index(1) if 1 in f else -1
    m = [n for n in same if n != first]; ok = [n for n in m if f[n] == 0]
    ns = [n for n in NONSH if f[n] == 1]
    print('%-28s 写像 %d、連鎖なし %3d 語（%.0f%%）、最初の連鎖語 %d。一致語 %d 語のうち v の連鎖なし %d（すべてなら矛盾なし）。窓の非共有 9 語のうち v の連鎖あり %d' % (name, mode, len(nc), 100 * len(nc) / 256, first, len(m), len(ok), len(ns)))
    return len(ok), len(m)
if len(sys.argv) > 1 and sys.argv[1] == 'null':
    rnd = random.Random(5); AL = b'abcdefghijklmnopqrstuvwxyz0123456789'; N = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    cnt = collections.Counter(); frac = []; bymode = collections.defaultdict(list)
    for _ in range(N):
        key = bytes(rnd.choice(AL) for _ in range(rnd.randrange(3, 13)))
        try: mode, f = flags(key)
        except ValueError: continue
        first = f.index(1) if 1 in f else -1; m = [n for n in same if n != first]; k = sum(1 for n in m if f[n] == 0)
        cnt[k] += 1; frac.append(f.count(0) / 256); bymode[mode].append(k)
    tot = sum(cnt.values())
    print('乱数の鍵 %d 本: 一致語（最初の連鎖語を除く 19〜20 語）のうち v の連鎖なしに入る数の分布: %s' % (tot, dict(sorted(cnt.items()))))
    print('   最大 %d、19 以上の鍵 %d 本。連鎖なしの割合 平均 %.2f（%.2f〜%.2f）。写像ごとの平均: %s' % (max(cnt), sum(v for k, v in cnt.items() if k >= 19), sum(frac) / len(frac), min(frac), max(frac), {m: round(sum(v) / len(v), 1) for m, v in sorted(bymode.items())}))
else:
    print('一致語 %d 語: %s' % (len(same), same))
    cands = [('123（既定の金庫の鍵）', b'123'), ('gao1234（1998 年スライド）', b'gao1234'), ('gao12345（同 DES）', b'gao12345'),
             ('読み込み中・・・（2.2 のキーFD）', '読み込み中・・・'.encode('cp932')), ('GoodLuckToIISI（4.02a）', b'GoodLuckToIISI'),
             ('GoodLuckToIISIGoodLuckToChangxingYang', b'GoodLuckToIISIGoodLuckToChangxingYang'), ('testkey1（対照。自分の試験用）', b'testkey1')]
    for name, key in cands: report(name, key)
