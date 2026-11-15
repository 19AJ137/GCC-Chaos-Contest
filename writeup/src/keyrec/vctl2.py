# 対照 2: 乱数で作った「鍵 v に似た鍵」（29〜40 バイト、写像 1）で、一部を未知にして解けるかを見る。
#   観測値は鍵 v と同じ条件にする: b0b1 は幅 183 の範囲、k780 は 2 通り（鏡像）、セルの幅は鍵 v と同じ。
#   python3 vctl2.py L(16|17) 本数 種 [形 …]     形 = nU,E の組（例 3,0 3,1 2,2）
import sys, random, time
from gks import key_schedule
from vpart import build_hyp, run_vres, synth_target, show
SJ2 = [x for x in range(0x40, 0xFD) if x != 0x7F]
def gen(kind, L, rng):
    if kind == 'ascii': return bytes(rng.randrange(0x20, 0x7F) for _ in range(L))
    if kind == 'alnum': return bytes(rng.choice(b'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789') for _ in range(L))
    if kind == 'sjis':
        out = b''
        while len(out) + 2 <= L: out += bytes([rng.choice(list(range(0x81, 0xA0)) + list(range(0xE0, 0xEB))), rng.choice(SJ2)])
        return out + bytes(rng.randrange(0x30, 0x7B) for _ in range(L - len(out)))
    if kind == 'klike':
        out = b''
        while len(out) + 6 + 4 <= L: out += bytes([rng.randrange(0x88, 0x98), rng.choice(SJ2), rng.randrange(0x88, 0x98), rng.choice(SJ2)]) + b'\x81\x45'
        return out + bytes(rng.randrange(0x30, 0x7B) for _ in range(L - len(out) - 4)) + b'IISI'
    if kind == 'mix':
        out = b''
        while len(out) < L:
            r = rng.random()
            if r < 0.5 and len(out) + 2 <= L: out += bytes([rng.choice(list(range(0x81, 0xA0))), rng.choice(SJ2)])
            elif r < 0.6: out += bytes([rng.randrange(0xA1, 0xE0)])
            else: out += bytes([rng.randrange(0x20, 0x7F)])
        return out
    if kind == 'any': return bytes(rng.randrange(0x20, 0xFD) for _ in range(L))
def pattern(key, nU, E, rng):
    """4 バイトそろった塊から nU 個を未知にし、ほかの塊の E バイトを未知にする（E バイトは 1 つの塊か 2 つの塊に置く）"""
    L = len(key); nfull = L // 4; nch = (L + 3) // 4
    us = set(rng.sample(range(nfull), nU)); unk = [False] * L
    for u in us:
        for b in range(4): unk[4 * u + b] = True
    rest = [c for c in range(nch) if c not in us]; e = E
    while e > 0:
        c = rng.choice(rest); pos = [4 * c + b for b in range(4) if 4 * c + b < L and not unk[4 * c + b]]
        if len(pos) <= 1: continue                       # 塊のうち少なくとも 1 バイトは既知に残す（全部未知にすると U か別の形になる）
        t = rng.randrange(1, min(e, len(pos) - 1) + 1)
        for p in rng.sample(pos, t): unk[p] = True
        e -= t
    segs = []; i = 0
    while i < L:
        j = i
        while j < L and unk[j] == unk[i]: j += 1
        segs.append((j - i) if unk[i] else key[i:j]); i = j
    return segs, us
if __name__ == '__main__':
    LAT = int(sys.argv[1]); NT = int(sys.argv[2]); seed = int(sys.argv[3]); shapes = [tuple(int(x) for x in a.split(',')) for a in sys.argv[4:]] or [(3, 0)]
    rng = random.Random(seed); kinds = ['ascii', 'alnum', 'sjis', 'klike', 'mix', 'any']
    cmin4 = 0x33 if LAT == 16 else 0x20
    tot = hit = cov = covhit = 0; T0 = time.time()
    for (nU, E) in shapes:
        for t in range(NT):
            kind = kinds[t % len(kinds)]
            while True:
                L = rng.randrange(29, 41); key = gen(kind, L, rng)
                if len(key) == L and key_schedule(key)['m'] == 1 and 0 not in key: break
            segs, us = pattern(key, nU, E, rng)
            tgt, r = synth_target(key, rng); h = build_hyp('s', segs, lat=LAT)
            ulist = [i for i, ln in enumerate(h['lines']) if ln == 'U']          # 入口が格子で解くことにした塊（未知のバイトを含む塊を U に上げることがある）
            incov = (len(ulist) <= 1 or all(key[4 * u + 3] >= cmin4 for u in ulist)) and all(0x20 <= b <= 0xFC for b in key)
            keys, stat, dt = run_vres(tgt, [h], L=LAT, tag='c2_%d' % seed)
            ks = keys.get('s', set()); ok = key in ks
            tot += 1; hit += ok; cov += incov; covhit += (ok and incov)
            print('未知: 塊 %d＋%d バイト → 解き方 U%d E%d  %-6s L %2d  %-58s %s%s  鍵 %d 本  %.1f s' % (nU, E, h['nU'], h['E'], kind, L, show(h)[:58], 'HIT ' if ok else 'MISS', '' if incov else '（範囲外）', len(ks), dt)); sys.stdout.flush()
    print('合計 %d 本: 当たり %d 本。範囲内 %d 本のうち当たり %d 本。%.0f 秒' % (tot, hit, cov, covhit, time.time() - T0))
