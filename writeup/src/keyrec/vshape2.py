# 形ごとの「和の条件に合う割合」の見積もり（11-14、3 回目）。
#   鍵の形（語 n 個・IISI・区切り）を決め、語を辞書の上位から乱数で入れたとき、塊ごとの s の和が鍵 v の目標（S − SJ、SJ は 90〜330）の
#   近くに来る割合を、一様な場合との比で出す。比が 0 に近い形は、その辞書の語では S に合わせられない（位置ごとのバイトの和が合わない）。
#   中間一致（kmitm）の「A が合う組の数」と同じ量を、乱数で粗く見積もるもの。近くの幅は ±half（既定 50,000）。
#   python3 vshape2.py n N2 N3 [区切り…]
import sys, os, random
import numpy as np
from gks import chunk_vals
from vpart import TARGET_V
from vdict import word_bytes, struct_list, last_slot_filter, chunk_slots, show_struct

def slot_arrays(tpl, slots, wbs):
    L = len(tpl); cs = chunk_slots(L, slots); afix = 0; arrs = []
    for i, s in enumerate(cs):
        if not s: afix += chunk_vals(int.from_bytes(tpl[4 * i:4 * i + 4], 'little'))[2]
        if len(s) > 1: return None, None
    for k, (o, n) in enumerate(slots):
        own = [i for i, s in enumerate(cs) if s == [k]]; a = np.zeros(len(wbs[k]), dtype=np.uint64)
        for x, w in enumerate(wbs[k]):
            key = bytearray(tpl); key[o:o + n] = w; t = 0
            for i in own: t += chunk_vals(int.from_bytes(bytes(key[4 * i:4 * i + 4]), 'little'))[2]
            a[x] = t
        arrs.append(a)
    return afix, arrs

def ratio(tpl, slots, wbs, M=20_000_000, half=50_000, rng=None):
    afix, arrs = slot_arrays(tpl, slots, wbs)
    if afix is None: return None
    rng = rng or np.random.default_rng(1); tot = np.full(M, afix, dtype=np.uint64)
    for a in arrs: tot += a[rng.integers(0, len(a), M)]
    d = (np.uint64(TARGET_V['S']) + np.uint64(1 << 34) - tot) & np.uint64(0xFFFFFFFF)          # 要る SJ（mod 2^32）
    dd = (d.astype(np.int64) - 210 + (1 << 31)) % (1 << 32) - (1 << 31)                      # 210 からのずれ（符号つき）
    out = []
    for h in (half, half // 10):
        cnt = int(np.count_nonzero(np.abs(dd) <= h)); exp = M * (2 * h + 1) / 2.0**32; out.append((cnt, exp, cnt / exp))
    return out

if __name__ == '__main__':
    n = int(sys.argv[1]); N2 = int(sys.argv[2]); N3 = int(sys.argv[3]); seps = sys.argv[4:] or ['dot']
    W2 = word_bytes(2, N2); W3 = word_bytes(3, N3); rows = []
    for sizes, pos, sn, tpl, slots in struct_list(n, seps):
        wbs = [W2 if s == 4 else W3 for s in sizes]
        if any(not w for w in wbs): continue
        wbs, mok = last_slot_filter(tpl, slots, wbs)
        if not mok: continue
        r = ratio(tpl, slots, wbs, M=int(float(os.environ.get('VSHAPE_M', '2e7'))))
        if r is None: continue
        rows.append((len(tpl), show_struct(sizes, pos, sn), r))
        print('%2d  %-40s  比 %.3f（±5 万: %d / %.0f）  比 %.3f（±5 千: %d / %.0f）' % (len(tpl), show_struct(sizes, pos, sn), r[0][2], r[0][0], r[0][1], r[1][2], r[1][0], r[1][1])); sys.stdout.flush()
