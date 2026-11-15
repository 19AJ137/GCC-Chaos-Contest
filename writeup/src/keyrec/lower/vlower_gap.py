#!/usr/bin/env python3
# 英小文字だけ・36〜38 文字の鍵は S が鍵 v の値に届かないことの、鍵スケジュールでの確かめ（docs/password.md §8.9）
#   python3 lower/vlower_gap.py [乱数の鍵の本数]
# (1) 乱数の英小文字の鍵で、S と目標の差（2^32 を法とした近いほう）の最小を見る。S が一様なら、N 本での最小は 2^31/N 前後になる。
# (2) σ が範囲にいちばん近くなる鍵（範囲のすぐ上・すぐ下）を作り、S の差を見る。差は K_PI × 隔たり の前後になるはず。
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gks import key_schedule, K_PI
from vlower import T, colcounts, wide, sigma_of

def dist(S): d = (S - T['S']) & 0xFFFFFFFF; return min(d, (1 << 32) - d)

def build(L, cols):
    """cols[j] = 列 j の文字の並び。鍵に組み立てる"""
    return bytes(cols[i & 3][i >> 2] for i in range(L))

def split_sum(total, m):
    """m 文字（a〜z）で和が total になる並びを 1 つ"""
    out = [97] * m; rem = total - 97 * m; i = 0
    assert 0 <= rem <= 25 * m
    while rem > 0: d = min(25, rem); out[i] += d; rem -= d; i += 1
    return out

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 1000000
    rng = random.Random(20261115)
    for L in (36, 37, 38, 39, 40):
        best = 1 << 32; hit = 0
        for _ in range(N):
            k = bytes(rng.randrange(97, 123) for _ in range(L))
            d = dist(key_schedule(k)['S'])
            if d < best: best = d
            if d == 0: hit += 1
        print('%d 文字・乱数の英小文字の鍵 %d 本: S と目標の差の最小 %d（S が一様なら %d 前後）、一致 %d 本' % (L, N, best, (1 << 31) // N, hit))
    for L in (36, 37, 38):
        cnt, n = colcounts(L); lo, hi = wide(n)
        res = []
        for side in ('上', '下'):
            cand = []
            for B3 in range(97 * cnt[3], 122 * cnt[3] + 1):
                for B2 in range(97 * cnt[2], 122 * cnt[2] + 1):
                    if side == '上':
                        s = (B3 << 24) + (B2 << 16) + (97 * cnt[1] << 8) + 97 * cnt[0]
                        if s > hi: cand.append((s, B3, B2))
                    else:
                        s = (B3 << 24) + (B2 << 16) + (122 * cnt[1] << 8) + 122 * cnt[0]
                        if s < lo: cand.append((s, B3, B2))
            s, B3, B2 = min(cand) if side == '上' else max(cand)
            c = 97 if side == '上' else 122
            key = build(L, [[c] * cnt[0], [c] * cnt[1], split_sum(B2, cnt[2]), split_sum(B3, cnt[3])])
            assert sigma_of(key) == s
            r = key_schedule(key)
            res.append('範囲のすぐ%s: %s（σ の隔たり %d、S の差 %d。K_PI × 隔たり = %d）' % (side, key.decode(), s - hi if side == '上' else lo - s, dist(r['S']), round(K_PI * (s - hi if side == '上' else lo - s))))
        print('%d 文字: ' % L + ' / '.join(res))
