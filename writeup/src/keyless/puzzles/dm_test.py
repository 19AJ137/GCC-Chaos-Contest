"""dm_test.py — B が同じ鍵での多重復号（D^m、m = 1・3・5・7・9）の場合の、鍵なしの式と検査（2026-11-12）
  連鎖ありの語を連鎖語の順に j = 0, 1, 2, … と番号を振る（連鎖なしの語は連鎖に関わらない）。S は「1 つ前の連鎖語へ」の移動。
    A_j = PA_j ⊕ K_j ⊕ A_(j−1)                      （通常の暗号化）       → K_j = s_j ⊕ PA_j、s_j = A_j ⊕ A_(j−1)（鍵なしで出る）
    B   = (1+S)^m PB ⊕ G_m K、 G_m = ((1+S)^m + 1)/S   （復号を m 回）
  (1+S)^(L−m) を掛けると (1+S)^L = 1 + S^L（L = 1・4・8・8・16）になり、
    Θ_j := [(1+S)^(L−m) B]_j ⊕ [T_m s]_j = PB_j ⊕ PB_(j−L) ⊕ ⊕{PA_(j−i) : i ∈ T_m}        （Θ は暗号文だけで出る）
    T_m = ((1+S)^L + (1+S)^(L−m))/S の項: T_1 = {0}、T_3 = {0, 3}、T_5 = {0, 1, 2, 7}、T_7 = {0, 7}、T_9 = {0, …, 6, 15}
  帰結:
    (a) A の平文が一致列 J で 7bit なら、MSB_J(Θ_j) = b_j ⊕ b_(j−L)（b = B の平文の一致列の最上位ビット）。L 語おきの組ごとに積み上げると、
        共有の連鎖語（b = 0）は、同じ組の中で同じしるしになる。
    (b) 同じ組（j の差が L の倍数）の共有の語 j1 < j2 では、⊕{Θ_j : j1 < j ≤ j2、同じ組} の一致列の最上位ビットが 0。
  python3 -I dm_test.py real [フラグ]       実データ（規則で予想される連鎖語を共有として (b) を数える。しるしの最頻値も出す）
  python3 -I dm_test.py syn [本数]          合成（写像 2 の鍵、A = Enc(PA)、B = D^m(PB)）で式と検査を確かめる"""
import sys, os, random, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D); sys.path.insert(0, os.path.join(D, '..', '..', 'gcc6'))
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
def W(d, i): return d[4 * i:4 * i + 4]
def X(*a):
    r = bytes(4)
    for x in a: r = bytes(p ^ q for p, q in zip(r, x))
    return r
Z4 = bytes(4)
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
def polymul(a, b):   # GF(2) の多項式（指数の集合）
    r = set()
    for i in a:
        for j in b: r ^= {i + j}
    return r
def onepS(n):
    r = {0}
    for _ in range(n): r = polymul(r, {0, 1})
    return r
LM = {1: 1, 3: 4, 5: 8, 7: 8, 9: 16}
def coeffs(m):
    L = LM[m]; M = onepS(L - m)                      # B に掛ける多項式
    num = onepS(L) ^ onepS(L - m)                    # (1+S)^L + (1+S)^(L−m)
    assert 0 not in num; T = set(i - 1 for i in num)  # /S
    return L, sorted(M), sorted(T)
def theta(Aw, Bw, m):
    """連鎖語の順に並べた A・B の語の列から Θ_j（出せない j は None）"""
    L, M, T = coeffs(m); n = len(Aw); s = [None] + [X(Aw[j], Aw[j - 1]) for j in range(1, n)]
    out = [None] * n; j0 = max(L, max(T) + 1, max(M))
    for j in range(j0, n):
        v = Z4
        for i in M: v = X(v, Bw[j - i])
        for i in T: v = X(v, s[j - i])
        out[j] = v
    return out
def msbJ(w, J): return ''.join('1' if w[c] & 0x80 else '0' for c in J)

def analyze(Aw, Bw, J, pred, label, quiet=False):
    """pred: 共有と決めた語の、連鎖語の順の番号。戻り値 {m: (成り立った関係, 関係の数)}"""
    res = {}
    for m in (1, 3, 5, 7, 9):
        L, M, T = coeffs(m); th = theta(Aw, Bw, m)
        ok = tot = 0; det = []
        for r in range(L):
            cls = [j for j in pred if j % L == r]
            for j1, j2 in zip(cls, cls[1:]):
                js = list(range(j1 + L, j2 + 1, L))
                if any(th[j] is None for j in js): continue
                v = Z4
                for j in js: v = X(v, th[j])
                tot += 1; z = msbJ(v, J) == '0' * len(J); ok += z; det.append('%d→%d:%s' % (j1, j2, msbJ(v, J)))
        # しるし（組ごとの積み上げ）の最頻値: 組ごとに、Θ が出る語をたどる
        modes = []
        for r in range(L):
            acc = Z4; sig = collections.Counter(); first = True
            for j in range(r, len(Aw), L):
                if th[j] is None:
                    if j + L < len(Aw) and th[j + L] is not None: sig[msbJ(Z4, J)] += 1   # 組の基準の語（しるし 0 とする）
                    continue
                acc = X(acc, th[j]); sig[msbJ(acc, J)] += 1
            if sig: modes.append((sum(sig.values()), max(sig.values())))
        res[m] = (ok, tot, det, modes)
        if not quiet:
            print('  %s D^%d（L = %d、B の係数 %s、T = %s）: 規則で共有と予想される語どうしの関係 %d 個のうち、一致列が 7bit のもの %d（偶然なら %.1f）  %s' % (label, m, L, M, T, tot, ok, tot / 8, ' '.join(det)))
            print('        組ごとの（語数, しるしの最頻値の語数）: %s' % modes)
    return res

def real(name):
    A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(CT, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
    fl = read_flags(os.path.join(D, 'real_flags_%s.txt' % name))
    tot = {m: [0, 0] for m in (1, 3, 5, 7, 9)}
    for wn, lo, hi, col in (('W2', 64, 119, 2), ('W1', 230, 255, 0)):
        hi = max(n for n in range(lo, hi + 1) if n in fl and all(k in fl for k in range(lo, n + 1)))
        ch = [n for n in range(lo, hi + 1) if fl[n] == 1]; J = [c for c in range(4) if c != col]
        pred = [j for j in range(1, len(ch)) if not (X(W(A, ch[j]), W(E, ch[j]), W(A, ch[j - 1]), W(E, ch[j - 1]))[col] & 0x80)]
        print('%s 語 %d〜%d: 連鎖語 %d、規則で共有と予想される語 %s（連鎖語の順の番号 %s）' % (wn, lo, hi, len(ch), [ch[j] for j in pred], pred))
        res = analyze([W(A, n) for n in ch], [W(B, n) for n in ch], J, pred, wn)
        for m in res: tot[m][0] += res[m][0]; tot[m][1] += res[m][1]
    print('両方の窓: ' + '、'.join('D^%d %d / %d' % (m, tot[m][0], tot[m][1]) for m in tot))

def syn(nkeys):
    import gcc
    rnd = random.Random(31); AL = b'abcdefghijklmnopqrstuvwxyz0123456789'
    agg = {mt: {m: [0, 0] for m in (1, 3, 5, 7, 9)} for mt in (1, 3, 5, 7, 9)}; ident = {m: [0, 0] for m in (1, 3, 5, 7, 9)}; rnd_rel = [0, 0]
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
        # A の平文: 一致列（0・1・3）は 7bit、位置 2 は 7 割 8bit。B の平文: 位置 2 が 7bit の行だけ A と同じ（規則）、ほかは 8bit の乱数
        PA = bytearray(); 
        for w in range(256):
            b2 = rnd.randrange(128, 256) if rnd.random() < 0.7 else rnd.randrange(32, 127)
            PA += bytes([rnd.randrange(32, 127), rnd.randrange(32, 127), b2, rnd.randrange(32, 127)])
        PA = bytes(PA); PB = bytearray(PA)
        for w in range(256):
            if PA[4 * w + 2] & 0x80: PB[4 * w:4 * w + 4] = bytes(rnd.randrange(256) for _ in range(4))
        PB = bytes(PB)
        A = gcc.file_cryption(key, PA)
        Aw = [W(A, n) for n in ch]; pred = [j for j in range(1, len(ch)) if not (PA[4 * ch[j] + 2] & 0x80)]
        for mt in (1, 3, 5, 7, 9):
            Bf = PB
            for _ in range(mt): Bf = gcc.file_cryption(key, Bf, encrypt=False)
            Bw = [W(Bf, n) for n in ch]
            # 式の確かめ（真の m で）
            L, M, T = coeffs(mt); th = theta(Aw, Bw, mt)
            for j in range(len(ch)):
                if th[j] is None or j < 2 * L + 20: continue      # 先頭の IV の影響が残る所は除く
                v = X(W(PB, ch[j]), W(PB, ch[j - L]))
                for i in T: v = X(v, W(PA, ch[j - i]))
                ident[mt][1] += 1; ident[mt][0] += (v == th[j])
            p2 = [j for j in pred if j >= 2 * LM[9] + 20]
            res = analyze(Aw, Bw, [0, 1, 3], p2, '', quiet=True)
            for m in res: agg[mt][m][0] += res[m][0]; agg[mt][m][1] += res[m][1]
    print('合成 %d 本（写像 2）。式 Θ_j = PB_j⊕PB_(j−L)⊕⊕PA の確かめ（真の m）: ' % nkeys + '、'.join('D^%d %d / %d' % (m, ident[m][0], ident[m][1]) for m in ident))
    print('規則どおりに共有する合成で、関係（一致列が 7bit）が成り立つ割合。行 = 作り方、列 = 検査の m')
    print('   作り方 | ' + ' | '.join('検査 D^%d' % m for m in (1, 3, 5, 7, 9)))
    for mt in agg:
        print('   B=D^%d  | ' % mt + ' | '.join('%4d/%4d=%.2f' % (agg[mt][m][0], agg[mt][m][1], agg[mt][m][0] / max(1, agg[mt][m][1])) for m in (1, 3, 5, 7, 9)))

if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'syn': syn(int(sys.argv[2]) if len(sys.argv) > 2 else 20)
    else: real(sys.argv[2] if len(sys.argv) > 2 else 'v2')
