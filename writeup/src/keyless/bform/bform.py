"""bform.py — B が 2 つの鍵（外側 k、内側 v）のどの順序・向きの処理かを、鍵なしの恒等式で判定する（2026-11-12）
  前提: A = E_k(E_v(M))（M は A の平文）。語ごとの連鎖フラグ f（k）・f'（v）と M が、語 0〜N-1 で分かっている。
  B の形の候補（E_k・D_k・E_v・D_v を 1〜3 回重ねたもの）ごとに、A・B の各語を
      鍵ストリーム K_n・K'_n、初期値 IV_k・IV_v、平文 M_n の XOR の式
  として記号で作り、未知数（K・K'・IV）が消える A・B の語の組み合わせ（= 鍵なしで成り立つ恒等式）を GF(2) の掃き出しで全部求める。
  それぞれの恒等式に実データ（32 ビットの語）を入れ、成り立つ数を数える。形が違えば 1 本ごとに確率 2^-32 でしか成り立たない。
  python3 -I bform.py head132_fixed.json [N] [-syn 形]   （-syn: 合成データ（gcc.py）でその形の B を作って同じ判定をする）"""
import sys, os, pickle, itertools, random
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D); sys.path.insert(0, os.path.join(D, '..', '..', 'gcc6'))
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
def W(d, i): return d[4 * i:4 * i + 4]
def X(*a):
    r = 0
    for x in a: r ^= x
    return r
def enc(I, fl, key):
    """暗号化の向き: O_n = I_n ⊕ K_n ⊕ f_n·(直前の連鎖語の出力、なければ IV)"""
    O = []; prev = None
    for n, s in enumerate(I):
        o = set(s) ^ {(key, n)}
        if fl[n]:
            o ^= ({('IV' + key,)} if prev is None else O[prev]); prev = n
        O.append(frozenset(o))
    return O
def dec(I, fl, key):
    """復号の向き: O_n = I_n ⊕ K_n ⊕ f_n·(直前の連鎖語の入力、なければ IV)"""
    O = []; prev = None
    for n, s in enumerate(I):
        o = set(s) ^ {(key, n)}
        if fl[n]:
            o ^= ({('IV' + key,)} if prev is None else I[prev]); prev = n
        O.append(frozenset(o))
    return O
OPS = {'Ek': (enc, 'K', 0), 'Dk': (dec, 'K', 0), 'Ev': (enc, 'V', 1), 'Dv': (dec, 'V', 1)}
def apply(form, N, flags):
    """form: 先に掛ける処理から順に並べたタプル（例 ('Ev','Ek') = E_k(E_v(M)))"""
    cur = [frozenset({('M', n)}) for n in range(N)]
    for op in form:
        fn, key, which = OPS[op]; cur = fn(cur, [flags[n][which] for n in range(N)], key)
    return cur
def identities(Asym, Bsym):
    """未知数（K・V・IV）が消える組み合わせ。返り値: [(A の語の集合, B の語の集合, M の語の集合)]"""
    N = len(Asym); rows = []
    for tag, S in (('A', Asym), ('B', Bsym)):
        for n in range(N):
            unk = frozenset(x for x in S[n] if x[0] != 'M'); ms = frozenset(x for x in S[n] if x[0] == 'M')
            rows.append([set(unk), {(tag, n)}, set(ms)])
    piv = {}; out = []
    for unk, obs, ms in rows:
        unk = set(unk); obs = set(obs); ms = set(ms)
        while unk:
            p = min(unk)
            if p in piv:
                u2, o2, m2 = piv[p]; unk ^= u2; obs ^= o2; ms ^= m2
            else:
                piv[p] = (set(unk), set(obs), set(ms)); break
        else:
            out.append((obs, ms))
    return out
def test(form, N, flags, Aw, Bw, Mw, known):
    """known: M が分かっている語の集合。恒等式のうち、分かっている M だけで書けるものを実データで確かめる"""
    Asym = apply(('Ev', 'Ek'), N, flags); Bsym = apply(form, N, flags)
    ids = identities(Asym, Bsym); ok = 0; tot = 0; triv = 0
    for obs, ms in ids:
        if any(m[1] not in known for m in ms): continue
        v = 0
        for tag, n in obs: v ^= (Aw if tag == 'A' else Bw)[n]
        for _, n in ms: v ^= Mw[n]
        tot += 1; ok += (v == 0)
        if len(obs) == 2 and not ms and {t for t, _ in obs} == {'A', 'B'} and len({n for _, n in obs}) == 1: triv += 1
    return len(ids), tot, ok, triv
def name(form): return ''.join({'Ek': 'E_k', 'Dk': 'D_k', 'Ev': 'E_v', 'Dv': 'D_v'}[o] + '∘' for o in reversed(form)).rstrip('∘')
def main():
    from headdata import load_head
    d = load_head(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 and not sys.argv[2].startswith('-') else 132
    flags = d['flags']; A = open(os.path.join(CT, 'A.bin'), 'rb').read(); B = open(os.path.join(CT, 'B.bin'), 'rb').read()
    Aw = [int.from_bytes(W(A, n), 'big') for n in range(N)]; Bw = [int.from_bytes(W(B, n), 'big') for n in range(N)]
    Mb = d['MA']; Mw = [int.from_bytes(W(Mb, n), 'big') for n in range(N)]; known = set(d.get('known', range(N)))
    forms = [f for L in (1, 2, 3) for f in itertools.product(('Ek', 'Dk', 'Ev', 'Dv'), repeat=L)]
    if '-syn' in sys.argv:
        import gcc
        fs = tuple(sys.argv[sys.argv.index('-syn') + 1].split(',')); rnd = random.Random(5)
        while True:
            kk = bytes(rnd.choice(b'abcdefghijklmnopqrstuvwxyz0123456789') for _ in range(8)); kv = bytes(rnd.choice(b'abcdefghijklmnopqrstuvwxyz0123456789') for _ in range(8))
            try: gcc.key_schedule(kk); gcc.key_schedule(kv); break
            except Exception: continue
        def layer(data, op):
            key = kk if op[1] == 'k' else kv; return gcc.file_cryption(key, data, encrypt=(op[0] == 'E'))
        M = Mb[:4 * N] + bytes(1024 - 4 * N)
        a = layer(layer(M, 'Ev'), 'Ek'); b = M
        for op in fs: b = layer(b, op)
        # 合成の本当のフラグ
        def fl(key):
            st = gcc.key_schedule(key); out = []
            for _ in range(256):
                ks, nc = gcc.ks_word(gcc.step(st)); out.append(0 if nc else 1)
            return out
        fk = fl(kk); fv = fl(kv); flags = {n: (fk[n], fv[n]) for n in range(256)}
        Aw = [int.from_bytes(W(a, n), 'big') for n in range(N)]; Bw = [int.from_bytes(W(b, n), 'big') for n in range(N)]; known = set(range(N))
        print('合成: 鍵 k=%s v=%s、B = %s、A = B の語 %d 個' % (kk.decode(), kv.decode(), name(fs), sum(1 for n in range(N) if Aw[n] == Bw[n])))
    print('語 0〜%d。M が分かっている語 %d 個。A = B の語: %s' % (N - 1, len(known), [n for n in range(N) if Aw[n] == Bw[n]]))
    print('  形（右から順に掛ける）            恒等式の数  うち確かめられる数  成り立つ数  （うち「A_n = B_n」だけの式）')
    res = []
    for f in forms:
        nid, tot, ok, triv = test(f, N, flags, Aw, Bw, Mw, known); res.append((f, nid, tot, ok, triv))
    for f, nid, tot, ok, triv in sorted(res, key=lambda r: (-(r[3] == r[2] and r[2] > 0), -r[3], r[1])):
        mark = '  ← 全部成り立つ' if tot > 0 and ok == tot else ''
        print('  %-34s %5d %10d %10d   （%d）%s' % (name(f), nid, tot, ok, triv, mark))
if __name__ == '__main__': main()
