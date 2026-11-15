"""en_test.py — B が同じ鍵での多重暗号化（E^n、n = 1・3・5・7・9）の場合の、鍵なしの式と検査（2026-11-12）。dm_test.py の対になるもの
  連鎖語の順に j を振る。 (1+S)^n B = PB ⊕ G_n K、G_n = ((1+S)^n + 1)/S、K_j = s_j ⊕ PA_j（s_j = A_j ⊕ A_(j−1)）から
    Ξ_j := [(1+S)^n B]_j ⊕ [G_n s]_j = PB_j ⊕ ⊕{PA_(j−i) : i ∈ G_n}          （Ξ は暗号文だけで出る）
  共有の語 j（PB_j = PA_j）では Ξ_j = ⊕{PA_(j−i) : i ∈ G_n、i ≠ 0}。n = 1 なら厳密に 0（差分の繰り返し）、n ≥ 3 なら一致列が 7bit（A の平文が一致列で 7bit のとき）。
  python3 -I en_test.py real [フラグ] ／ python3 -I en_test.py syn [本数]"""
import sys, os, random
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D); sys.path.insert(0, os.path.join(D, '..', '..', 'gcc6'))
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
from dm_test import W, X, Z4, read_flags, onepS, msbJ
def gset(n):
    num = onepS(n) ^ {0}; assert 0 not in num; return sorted(i - 1 for i in num)
def xi(Aw, Bw, n):
    M = sorted(onepS(n)); G = gset(n); N = len(Aw); s = [None] + [X(Aw[j], Aw[j - 1]) for j in range(1, N)]
    out = [None] * N
    for j in range(max(max(M), max(G) + 1), N):
        v = Z4
        for i in M: v = X(v, Bw[j - i])
        for i in G: v = X(v, s[j - i])
        out[j] = v
    return out
def real(name):
    A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
    fl = read_flags(os.path.join(D, 'real_flags_%s.txt' % name)); tot = {n: [0, 0, 0, 0] for n in (1, 3, 5, 7, 9)}
    for wn, lo, hi, col in (('W2', 64, 119, 2), ('W1', 230, 255, 0)):
        hi = max(n for n in range(lo, hi + 1) if n in fl and all(k in fl for k in range(lo, n + 1)))
        ch = [n for n in range(lo, hi + 1) if fl[n] == 1]; J = [c for c in range(4) if c != col]
        pred = [j for j in range(1, len(ch)) if not (X(W(A, ch[j]), W(E, ch[j]), W(A, ch[j - 1]), W(E, ch[j - 1]))[col] & 0x80)]
        for n in (1, 3, 5, 7, 9):
            x = xi([W(A, m) for m in ch], [W(B, m) for m in ch], n)
            ok = [j for j in pred if x[j] is not None and msbJ(x[j], J) == '000']; val = [j for j in pred if x[j] is not None]
            allj = [j for j in range(len(ch)) if x[j] is not None]; ok2 = [j for j in allj if msbJ(x[j], J) == '000']
            z = sum(1 for j in allj if x[j] == Z4)
            print('  %s E^%d: 規則で共有と予想される語 %d 個のうち、Ξ の一致列が 7bit のもの %d（偶然なら %.1f）。全連鎖語では %d / %d（偶然なら %.1f）、Ξ = 0 の語 %d' % (wn, n, len(val), len(ok), len(val) / 8, len(ok2), len(allj), len(allj) / 8, z))
            t = tot[n]; t[0] += len(ok); t[1] += len(val); t[2] += len(ok2); t[3] += len(allj)
    print('両方の窓（予想される語）: ' + '、'.join('E^%d %d / %d' % (n, tot[n][0], tot[n][1]) for n in tot))
def syn(nkeys):
    import gcc
    rnd = random.Random(41); AL = b'abcdefghijklmnopqrstuvwxyz0123456789'; agg = {nt: {n: [0, 0] for n in (1, 3, 5, 7, 9)} for nt in (1, 3, 5, 7, 9)}; ident = {n: [0, 0] for n in (1, 3, 5, 7, 9)}
    for trial in range(nkeys):
        while True:
            key = bytes(rnd.choice(AL) for _ in range(rnd.randrange(6, 13)))
            try: st = gcc.key_schedule(key)
            except ValueError: continue
            if st.mode == 2: break
        st = gcc.key_schedule(key); fl = []
        for n in range(256):
            x = gcc.step(st); k, nc = gcc.ks_word(x); fl.append(0 if nc else 1)
        ch = [n for n in range(256) if fl[n] == 1]
        PA = bytearray()
        for w in range(256):
            b2 = rnd.randrange(128, 256) if rnd.random() < 0.7 else rnd.randrange(32, 127)
            PA += bytes([rnd.randrange(32, 127), rnd.randrange(32, 127), b2, rnd.randrange(32, 127)])
        PA = bytes(PA); PB = bytearray(PA)
        for w in range(256):
            if PA[4 * w + 2] & 0x80: PB[4 * w:4 * w + 4] = bytes(rnd.randrange(256) for _ in range(4))
        PB = bytes(PB); A = gcc.file_cryption(key, PA); Aw = [W(A, n) for n in ch]
        pred = [j for j in range(40, len(ch)) if not (PA[4 * ch[j] + 2] & 0x80)]
        for nt in (1, 3, 5, 7, 9):
            Bf = PB
            for _ in range(nt): Bf = gcc.file_cryption(key, Bf, encrypt=True)
            Bw = [W(Bf, n) for n in ch]
            x = xi(Aw, Bw, nt); G = gset(nt)
            for j in range(40, len(ch)):
                v = W(PB, ch[j])
                for i in G: v = X(v, W(PA, ch[j - i]))
                ident[nt][1] += 1; ident[nt][0] += (v == x[j])
            for n in (1, 3, 5, 7, 9):
                xx = xi(Aw, Bw, n)
                agg[nt][n][0] += sum(1 for j in pred if msbJ(xx[j], [0, 1, 3]) == '000'); agg[nt][n][1] += len(pred)
    print('合成 %d 本。式 Ξ_j = PB_j⊕⊕PA の確かめ（真の n）: ' % nkeys + '、'.join('E^%d %d / %d' % (n, ident[n][0], ident[n][1]) for n in ident))
    print('規則どおりに共有する合成で、予想される語の Ξ が一致列で 7bit になる割合。行 = 作り方、列 = 検査の n')
    for nt in agg: print('   B=E^%d | ' % nt + ' | '.join('%4d/%4d=%.2f' % (agg[nt][n][0], agg[nt][n][1], agg[nt][n][0] / max(1, agg[nt][n][1])) for n in (1, 3, 5, 7, 9)))
if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'syn': syn(int(sys.argv[2]) if len(sys.argv) > 2 else 20)
    else: real(sys.argv[2] if len(sys.argv) > 2 else 'v2')
