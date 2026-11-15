"""layers2.py — 2 本の鍵（懸賞の鍵 k と別の鍵 v）の暗号化・復号を 1〜4 回重ねた形を全部作り、A = E_k(PA) との関係を合成で数える（2026-11-12）
  B = L_n(…L_2(L_1(PB))…)、L_i ∈ {E_k, D_k, E_v, D_v}。PB は (a) PA と全部同じ、(b) 語の 4 割だけ同じ（残りは乱数）の 2 通り。
  数えるもの（実データの値）: 一致語（20）、そのうち k の連鎖ありの語（先頭の連鎖語を除いて 0〜2）、
     連鎖ありの語で Y_n = Y_p（差分の繰り返し。0）、Y_n = A_p（0）。 Y = A⊕B、p は k の直前の連鎖語
  python3 -I layers2.py [鍵の組の数（既定 12）] [最大の層の数（既定 4）]"""
import sys, os, random, itertools, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(D, '..', '..', 'gcc6'))
import gcc
def W(d, i): return d[4 * i:4 * i + 4]
def X(*a):
    r = bytes(4)
    for x in a: r = bytes(p ^ q for p, q in zip(r, x))
    return r
def flags(key):
    st = gcc.key_schedule(key); f = []
    for n in range(256):
        x = gcc.step(st); k, nc = gcc.ks_word(x); f.append(0 if nc else 1)
    return st.mode, f
def reduce_seq(seq):
    """隣り合う E_x・D_x（同じ鍵の逆の処理）を消す。残りが空なら恒等"""
    out = []
    for op in seq:
        if out and out[-1][1] == op[1] and out[-1][0] != op[0]: out.pop()
        else: out.append(op)
    return tuple(out)
def main():
    nk = int(sys.argv[1]) if len(sys.argv) > 1 else 12; maxl = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    rnd = random.Random(77); AL = b'abcdefghijklmnopqrstuvwxyz0123456789'
    def newkey(mode=None):
        while True:
            key = bytes(rnd.choice(AL) for _ in range(rnd.randrange(6, 13)))
            try: m, f = flags(key)
            except ValueError: continue
            if (mode is None or m == mode) and 0.12 < f.count(0) / 256 < 0.6: return key, f
    ops = [('E', 'k'), ('D', 'k'), ('E', 'v'), ('D', 'v')]
    seqs = set()
    for L in range(1, maxl + 1):
        for s in itertools.product(ops, repeat=L):
            r = reduce_seq(s)
            if r and r == s: seqs.add(s)          # 消せる並びは、短い並びとして数える
    seqs = sorted(seqs, key=lambda s: (len(s), s))
    stat = {(s, pb): collections.Counter() for s in seqs for pb in 'ab'}
    for t in range(nk):
        k, fk = newkey(2); v, fv = newkey()
        PA = bytes(rnd.randrange(256) for _ in range(1024))
        PBb = bytearray(PA)
        for w in range(256):
            if rnd.random() >= 0.4: PBb[4 * w:4 * w + 4] = bytes(rnd.randrange(256) for _ in range(4))
        A = gcc.file_cryption(k, PA); ch = [n for n in range(256) if fk[n] == 1]; prev = {ch[i]: ch[i - 1] for i in range(1, len(ch))}
        for pb, P in (('a', PA), ('b', bytes(PBb))):
            for s in seqs:
                Bf = P
                for op, key in s: Bf = gcc.file_cryption(k if key == 'k' else v, Bf, encrypt=(op == 'E'))
                same = [n for n in range(1, 256) if W(A, n) == W(Bf, n)]
                Y = [X(W(A, n), W(Bf, n)) for n in range(256)]
                c = stat[(s, pb)]
                c['n'] += 1; c['same'] += len(same); c['same_kch'] += sum(1 for n in same if fk[n] == 1 and n != ch[0]); c['same_vch'] += sum(1 for n in same if fv[n] == 1)
                c['rep'] += sum(1 for n in prev if Y[n] == Y[prev[n]] and any(Y[n])); c['ap'] += sum(1 for n in prev if Y[n] == W(A, prev[n]))
                c['knc'] += fk.count(0); c['both'] += sum(1 for n in range(256) if fk[n] == 0 and fv[n] == 0)
    def name(s): return ' '.join('%s_%s' % (op, key) for op, key in reversed(s))      # 外側（最後にかけた処理）を左に書く
    print('鍵の組 %d、並び %d 通り（同じ鍵の E・D が隣り合って消えるものは除く）。値は 1 本あたりの平均。実データ: 一致語 20（語 0 を除いて 19）、k の連鎖ありの一致語 0〜2、繰り返し 0、Y=A_p 0' % (nk, len(seqs)))
    print('   連鎖なしの語（k）の平均 %.1f、k・v とも連鎖なしの語の平均 %.1f' % (stat[(seqs[0], 'a')]['knc'] / nk, stat[(seqs[0], 'a')]['both'] / nk))
    rows = []
    for s in seqs:
        a = stat[(s, 'a')]; b = stat[(s, 'b')]; n = a['n']
        rows.append((name(s), len(s), sum(1 for o in s if o[1] == 'v'), a['same'] / n, a['same_kch'] / n, a['rep'] / n, a['ap'] / n, b['same'] / n, b['same_kch'] / n, b['rep'] / n, b['ap'] / n))
    def verdict(r):
        nm, L, nv, sa, ska, repa, apa, sb, skb, repb, apb = r
        out = []
        for lab, s_, sk, rep, ap in (('平文が同じ', sa, ska, repa, apa), ('4 割共有', sb, skb, repb, apb)):
            if s_ < 3: out.append('%s: 一致語が出ない' % lab)
            elif sk > 3: out.append('%s: k の連鎖ありの語でも一致' % lab)
            elif rep > 3: out.append('%s: 繰り返しが出る' % lab)
            elif ap > 3: out.append('%s: Y=A_p が出る' % lab)
            elif s_ > 45: out.append('%s: 一致語が多すぎる（%.0f）' % (lab, s_))
            else: out.append('%s: 合う（一致語 %.1f）' % (lab, s_))
        return '、'.join(out)
    groups = collections.defaultdict(list)
    for r in rows: groups[verdict(r)].append(r)
    for vd, rs in sorted(groups.items(), key=lambda kv: -len(kv[1])):
        print('\n[%d 通り] %s' % (len(rs), vd))
        ex = sorted(rs, key=lambda r: (r[1], r[0]))
        for r in ex[:6]:
            print('     %-34s 平文が同じ: 一致 %5.1f（k 連鎖あり %4.1f）繰り返し %4.1f Y=A_p %4.1f ｜ 4 割共有: 一致 %5.1f（%4.1f）繰り返し %4.1f Y=A_p %4.1f' % (r[0], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10]))
        if len(ex) > 6: print('     … ほか %d 通り' % (len(ex) - 6))
main()
