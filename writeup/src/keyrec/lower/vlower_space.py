#!/usr/bin/env python3
# 英小文字と空白だけの鍵（英文の形）で、σ が鍵 v の範囲に入るときの、列ごとの空白の数（docs/password.md §8.9）
#   python3 lower/vlower_space.py
# 列 j = 鍵の j+1, j+5, j+9, … 文字目（4 文字目の列は 4・8・12… 文字目）。文字は英語の頻度で独立、空白は割合 p で独立に入るとみた粗い型。
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from vlower import refined, colcounts, ways
from vlower_text import window_prob, letter_models

if __name__ == '__main__':
    uni, en, ro, words = letter_models()
    avg = sum(f * len(w) for w, f in words) / sum(f for w, f in words); p = 1 / (1 + avg)
    print('語の長さの平均 %.2f 文字、空白の割合 p = %.3f' % (avg, p))
    SMAX = 4
    for L in range(36, 41):
        cnt, n = colcounts(L); win = refined(n)
        Ws = []                                   # Ws[j][s] = 列 j に空白が s 個あるときの、和の分布（重みつき）
        for j in range(4):
            m = cnt[j]; row = []
            for s in range(SMAX + 1):
                w = ways(en, m - s) * math.comb(m, s) * p**s * (1 - p)**(m - s)
                row.append(np.concatenate([np.zeros(32 * s), w]))
            Ws.append(row)
        tot = 0.0; marg = [[0.0] * (SMAX + 1) for _ in range(4)]; best = []
        for s3 in range(SMAX + 1):
            for s2 in range(SMAX + 1):
                for s1 in range(SMAX + 1):
                    for s0 in range(SMAX + 1):
                        q = window_prob([Ws[0][s0], Ws[1][s1], Ws[2][s2], Ws[3][s3]], win)
                        if q == 0: continue
                        tot += q; best.append((q, (s0, s1, s2, s3)))
                        for j, s in enumerate((s0, s1, s2, s3)): marg[j][s] += q
        best.sort(reverse=True)
        print('%d 文字: σ が範囲に入る割合 2^%.1f' % (L, math.log2(tot)))
        for j, nm in ((3, '4 文字目の列'), (2, '3 文字目の列'), (1, '2 文字目の列'), (0, '1 文字目の列')):
            print('   %s（%d 字）の空白の数: ' % (nm, cnt[j]) + '、'.join('%d 個 %.0f%%' % (s, 100 * marg[j][s] / tot) for s in range(SMAX + 1) if marg[j][s] / tot >= 0.005))
        print('   多い組み合わせ（1・2・3・4 文字目の列の空白の数）: ' + '、'.join('%s %.0f%%' % ('-'.join(map(str, c)), 100 * q / tot) for q, c in best[:4]))
