#!/usr/bin/env python3
# 英小文字だけの鍵（a〜z）が、鍵 v の S・b0b1 の条件に合うかを長さごとに調べる（docs/password.md §8.9）
#   python3 lower/vlower.py            … 塊の値の和 σ の範囲、既知の等価な鍵での確かめ、文字の種類 × 長さの表、英小文字 39・40 文字の形
# 鍵の 4 バイトの塊（little-endian）を c、その和を σ と書く。S は K_PI·σ でほぼ決まる（K_PI·σ = S + w·2^32 − Y）。
# σ = B0 + 256·B1 + 65536·B2 + 2^24·B3（B_j = 各塊の j 番目のバイトの和。端数の塊は j < 端数の長さ の列にだけ入る）。
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from fractions import Fraction as Fr
from sums import sigma_window, TG
from gks import K_PI, chunks_of

T = TG['v']
KF = Fr(K_PI)

def refined(n):
    """S と b0b1 の両方に合う σ の範囲（塊 n 個。sums.sigma_window）"""
    w = sigma_window(T['S'], T['Tlo'], T['Thi'], n)
    sig = sorted(set(x[1] for x in w))
    return (sig[0], sig[-1]) if sig else None

def wide(n, w=1):
    """S だけから言える σ の範囲（必要条件）。Y = R + ΣJ − Σf1 −（切り捨て）は 0 より大きく、3n + 4(2^n − 1) + 8 以下"""
    ymax = 3 * n + 4 * (2 ** n - 1) + 8
    lo = math.ceil(Fr(T['S'] + (w << 32) - ymax) / KF); hi = math.floor(Fr(T['S'] + (w << 32)) / KF)
    return lo, hi

def colcounts(L):
    nf, r = divmod(L, 4)
    return [nf + (1 if j < r else 0) for j in range(4)], nf + (1 if r else 0)

def ways(pm, m):
    """pm: バイト値ごとの重み（長さ 256）。m 文字の和の分布（重みの畳み込み）"""
    out = np.array([1.0])
    for _ in range(m): out = np.convolve(out, pm)
    return out

def count_window(pm, L, win):
    """長さ L の鍵で σ が win=(lo,hi) に入るものの重みの合計と、(B3,B2,B1) ごとの内訳"""
    cnt, n = colcounts(L)
    W = [ways(pm, m) for m in cnt]
    lo, hi = win
    P0 = np.concatenate([[0.0], np.cumsum(W[0])])        # P0[k] = B0 < k の重み
    def r0(a, b):                                          # B0 が a..b の重み
        a = max(a, 0); b = min(b, len(W[0]) - 1)
        return P0[b + 1] - P0[a] if b >= a else 0.0
    tot = 0.0; parts = {}
    nz = [np.nonzero(w)[0] for w in W]
    for B3 in nz[3]:
        base3 = int(B3) << 24
        if base3 > hi: break
        b2lo = max(0, (lo - base3 - 256 * (len(W[1]) - 1) - (len(W[0]) - 1)) // 65536); b2hi = (hi - base3) // 65536
        for B2 in range(b2lo, min(b2hi, len(W[2]) - 1) + 1):
            if W[2][B2] == 0: continue
            base2 = base3 + (B2 << 16)
            for B1 in nz[1]:
                base1 = base2 + (int(B1) << 8)
                if base1 > hi: break
                v = r0(lo - base1, hi - base1)
                if v > 0:
                    p = W[3][B3] * W[2][B2] * W[1][B1] * v
                    tot += p; parts[(int(B3), B2, int(B1))] = p
    return tot, parts, cnt, n

def pm_of(chars):
    pm = np.zeros(256)
    for c in chars: pm[c] = 1.0
    return pm

LOW = list(range(97, 123)); UP = list(range(65, 91)); DIG = list(range(48, 58))
ALPH = [('a〜z', LOW), ('a〜z と空白', LOW + [32]), ('a〜z と数字', LOW + DIG), ('A〜Z と a〜z', UP + LOW),
        ('英数字', UP + LOW + DIG), ('英数字と空白', UP + LOW + DIG + [32]), ('印字可能な ASCII', list(range(32, 127)))]

def sigma_of(key): return sum(chunks_of(key))

if __name__ == '__main__':
    print('== 塊の値の和 σ の範囲（鍵 v）')
    for n in range(8, 12):
        r = refined(n); w = wide(n)
        print('  塊 %2d 個: S と b0b1 から %d〜%d（%d 通り）、S だけから %d〜%d（%d 通り）' % (n, r[0], r[1], r[1] - r[0] + 1, w[0], w[1], w[1] - w[0] + 1))
    print('  折り返し 0 回・2 回なら σ はおよそ %.4g・%.4g（英小文字 9〜10 塊の和は %.4g〜%.4g なので届かない）' % (
        float(Fr(T['S']) / KF), float(Fr(T['S'] + (2 << 32)) / KF), 9 * 0x61616161, 10 * 0x7a7a7a7a))

    print('== 既知の等価な鍵の σ（範囲の確かめ）')
    keys = []
    for fn in ('short_v_keys.txt', 'forged_v_keys.txt'):
        for line in open(fn, encoding='utf-8'):
            if line.startswith('#') or not line.strip(): continue
            keys.append(bytes.fromhex(line.split('|')[0].strip()))
    bad = 0; byn = {}
    for k in keys:
        n = (len(k) + 3) // 4; s = sigma_of(k); r = refined(n)
        byn.setdefault(n, []).append(s)
        if not (r[0] <= s <= r[1]): bad += 1
    print('  鍵 %d 本（塊 %s 個）: 範囲の外 %d 本。σ は %d〜%d' % (len(keys), '・'.join(str(n) for n in sorted(byn)), bad, min(min(v) for v in byn.values()), max(max(v) for v in byn.values())))

    print('== σ が範囲に入る鍵の数（log2。「なし」は 1 本もない。S と b0b1 から出した範囲。括弧は S だけから出した広い範囲）')
    Ls = list(range(29, 45))
    print('  %-18s' % '長さ' + ''.join('%7d' % L for L in Ls))
    for name, chars in ALPH:
        pm = pm_of(chars); row = []; roww = []
        for L in Ls:
            cnt, n = colcounts(L)
            t, _, _, _ = count_window(pm, L, refined(n)); tw, _, _, _ = count_window(pm, L, wide(n))
            row.append('   なし' if t == 0 else '%7.1f' % math.log2(t)); roww.append('   なし' if tw == 0 else '%7.1f' % math.log2(tw))
        print('  %-18s' % name + ''.join(row)); print('  %-18s' % '  （S だけ）' + ''.join(roww))

    print('== 英小文字だけ: 成り立つ長さでの、列ごとの和（S と b0b1 から出した範囲）')
    pm = pm_of(LOW)
    for L in Ls:
        cnt, n = colcounts(L)
        t, parts, cnt, n = count_window(pm, L, refined(n))
        if t == 0: continue
        print('  %d 文字（4 文字そろった塊 %d 個＋端数 %d 文字）: 2^%.1f 本。全部の小文字の列 26^%d = 2^%.1f のうち 2^%.1f 分の 1' % (L, L // 4, L % 4, math.log2(t), L, L * math.log2(26), L * math.log2(26) - math.log2(t)))
        for (B3, B2, B1), p in sorted(parts.items(), key=lambda x: -x[1])[:6]:
            print('     4 文字目の列の和 %d（%d 字、平均 %.1f = %s）、3 文字目の列 %d（%d 字、平均 %.1f = %s）、2 文字目の列 %d（%d 字、平均 %.1f = %s）: 割合 %.3f' % (
                B3, cnt[3], B3 / cnt[3], chr(round(B3 / cnt[3])), B2, cnt[2], B2 / cnt[2], chr(round(B2 / cnt[2])), B1, cnt[1], B1 / cnt[1], chr(round(B1 / cnt[1])), p / t))
    print('== 英小文字だけ・36〜38 文字: いちばん近い σ までの隔たり（S だけから出した範囲に対して）')
    for L in (36, 37, 38):
        cnt, n = colcounts(L); lo, hi = wide(n)
        # 範囲のすぐ下とすぐ上にある、とりうる σ
        below = max(((B3 << 24) + (B2 << 16) + (122 * cnt[1] << 8) + 122 * cnt[0]) for B3 in range(97 * cnt[3], 122 * cnt[3] + 1) for B2 in range(97 * cnt[2], 122 * cnt[2] + 1)
                    if (B3 << 24) + (B2 << 16) + (122 * cnt[1] << 8) + 122 * cnt[0] < lo)
        above = min(((B3 << 24) + (B2 << 16) + (97 * cnt[1] << 8) + 97 * cnt[0]) for B3 in range(97 * cnt[3], 122 * cnt[3] + 1) for B2 in range(97 * cnt[2], 122 * cnt[2] + 1)
                    if (B3 << 24) + (B2 << 16) + (97 * cnt[1] << 8) + 97 * cnt[0] > hi)
        print('  %d 文字: 範囲 %d〜%d。とりうる σ は、下が %d 以下（隔たり %d）、上が %d 以上（隔たり %d）' % (L, lo, hi, below, lo - below, above, above - hi))
