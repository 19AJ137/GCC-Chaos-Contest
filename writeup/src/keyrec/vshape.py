# 鍵 v の「位置ごとのバイトの和」と、鍵の形（末尾・文字の種類）の制約。python3 vshape.py
import io, contextlib, math
with contextlib.redirect_stdout(io.StringIO()):
    from possum import windows
from gks import chunk_vals, map_of_i2

W = windows('v', range(8, 13))
print('sigma の候補（鍵 v）')
for n, (sig, wr, sjlo, sjhi) in W.items():
    print('  n=%2d: %d..%d (%d 通り)  sigma/2^24 = %.3f..%.3f' % (n, sig[0], sig[-1], len(sig), sig[0] / 2**24, sig[-1] / 2**24))

def b3_set(sig_list, nfull, part, rng):
    """4・8・…バイト目の和 B3 の候補。rng[j] = (lo, hi) は j 番目のバイトの値の範囲（塊ごと）。端数の塊は part バイト"""
    out = set()
    cnt = [nfull + (1 if j < part else 0) for j in range(4)]
    for s in sig_list:
        for B0 in range(cnt[0] * rng[0][0], cnt[0] * rng[0][1] + 1):
            if (s - B0) % 256: continue
            r1 = (s - B0) // 256
            for B1 in range(cnt[1] * rng[1][0], cnt[1] * rng[1][1] + 1):
                if (r1 - B1) % 256: continue
                r2 = (r1 - B1) // 256
                for B2 in range(cnt[2] * rng[2][0], cnt[2] * rng[2][1] + 1):
                    if (r2 - B2) % 256: continue
                    B3 = (r2 - B2) // 256
                    if cnt[3] * rng[3][0] <= B3 <= cnt[3] * rng[3][1]: out.add(B3)
    return sorted(out)

kinds = {
 '印字可能な ASCII（0x20〜0x7E）': [(0x20, 0x7E)] * 4,
 '全角だけ（第 1 バイト 0x81〜0x9F、第 2 バイト 0x40〜0xFC）': [(0x81, 0x9F), (0x40, 0xFC), (0x81, 0x9F), (0x40, 0xFC)],
 '任意のバイト（0x20〜0xFC）': [(0x20, 0xFC)] * 4,
}
print('\n4・8・12…バイト目の和 B3 の候補（4 バイトそろった塊だけを数える）')
for nfull, part in ((7, 1), (7, 2), (7, 3), (8, 0), (8, 1), (8, 2), (8, 3), (9, 0), (9, 1), (9, 2), (9, 3), (10, 0)):
    n = nfull + (1 if part else 0)
    if n not in W: continue
    sig = W[n][0]
    L = 4 * nfull + part
    for name, rng in kinds.items():
        if '全角' in name and L % 2: continue
        b3 = b3_set(sig, nfull, part, rng)
        if b3: print('  %2d バイト（塊 %d 個＋端数 %d）%s: B3 = %s、平均 %.1f' % (L, nfull, part, name, ('%d〜%d' % (b3[0], b3[-1])) if len(b3) > 1 else str(b3[0]), sum(b3) / len(b3) / nfull))

# 漢字の第 2 バイトが一様（0x40〜0x7E・0x80〜0xFC の 188 通り）としたときの、和が 987 以下になる確率
vals = list(range(0x40, 0x7F)) + list(range(0x80, 0xFD))
mu = sum(vals) / len(vals); var = sum(v * v for v in vals) / len(vals) - mu * mu
print('\n漢字の第 2 バイト: 平均 %.1f、標準偏差 %.1f' % (mu, math.sqrt(var)))
for nfull in (7, 8, 9, 10):
    z = (987.5 - nfull * mu) / math.sqrt(nfull * var)
    p = 0.5 * math.erfc(-z / math.sqrt(2))
    print('  塊 %d 個の 4 バイト目が全部漢字の第 2 バイト: 和の平均 %.0f、987 以下になる確率 %.2g（z = %.2f）' % (nfull, nfull * mu, p, z))
print('  参考: カタカナの第 2 バイト 0x40〜0x96（平均 %.1f）、ひらがな 0x9F〜0xF1（平均 %.1f）、「・」0x45 = 69、英小文字の平均は 100〜110' % ((0x40 + 0x96) / 2, (0x9F + 0xF1) / 2))

print('\n鍵が IISI で終わるときの、最後の塊と写像の番号')
def m_of(b):
    c = int.from_bytes(b.ljust(4, b'\0'), 'little'); (f, i, s, t) = chunk_vals(c); return map_of_i2(i[1])
for r, last in ((0, b'IISI'), (1, b'I'), (2, b'SI'), (3, b'ISI')):
    print('  長さ ≡ %d (mod 4): 最後の塊 %-6s → 写像 %d' % (r, last.decode(), m_of(last)))
