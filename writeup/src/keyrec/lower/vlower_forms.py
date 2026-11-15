#!/usr/bin/env python3
# 参考: 鍵 k と同じ形（全角の語を「・」でつなぐ。IISI あり・なし）で、σ が鍵 v の範囲に入る割合（docs/password.md §8.9 の (4)）
#   python3 lower/vlower_forms.py [標本の数]
# 英文の形（lower/vlower_text.py）と同じ数え方: 列ごとの和の分布を標本から作り、列どうしは独立とみる。語は頻度に比例して選ぶ。
import sys, os, math, re
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from vlower import refined, colcounts
from vlower_text import window_prob

def ja_words(n2=11000, n3=5000):
    import wordfreq
    fd = wordfreq.get_frequency_dict('ja', 'large'); out = {2: [], 3: []}
    for w, f in sorted(fd.items(), key=lambda x: -x[1]):
        if len(w) in (2, 3) and (re.fullmatch('[一-鿿]+', w) or re.fullmatch('[ァ-ヶー]+', w)):
            try: b = w.encode('cp932')
            except UnicodeEncodeError: continue
            if len(b) == 2 * len(w) and len(out[len(w)]) < (n2 if len(w) == 2 else n3): out[len(w)].append((b, f))
    return out[2] + out[3]

def col_pmfs_bytes(sampler, L, N, rng):
    cnt, n = colcounts(L)
    W = [np.zeros(256 * cnt[j] + 1) for j in range(4)]
    for _ in range(N):
        s = sampler(L, rng); b = [0, 0, 0, 0]
        for i, ch in enumerate(s): b[i & 3] += ch
        for j in range(4): W[j][b[j]] += 1
    return [w / N for w in W]

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
    rng = np.random.default_rng(20261115)
    words = ja_words(); ws = [w for w, f in words]; cum = np.cumsum([f for w, f in words]); cum /= cum[-1]
    SEP = '・'.encode('cp932')
    def form(L, rng, iisi):
        while True:
            el = []; n = 4 + 2 if iisi else 0          # IISI と、その区切り 1 つ分
            while n < L + 2:
                w = ws[int(np.searchsorted(cum, rng.random()))]; el.append(w); n += len(w) + 2
            if n == L + 2:
                if iisi: el.insert(int(rng.integers(0, len(el) + 1)), b'IISI')
                return SEP.join(el)
    print('語は 2 文字の上位 11,000 と 3 文字の上位 5,000（漢字だけ・カタカナだけ。頻度に比例して選ぶ）。標本 %d 本ずつ' % N)
    for name, iisi in (('全角の語を「・」でつなぐ（IISI なし）', False), ('全角の語と IISI を「・」でつなぐ（IISI の位置は乱数）', True)):
        row = []
        for L in (34, 36, 38, 40):
            W = col_pmfs_bytes(lambda L, r: form(L, r, iisi), L, N, rng); p = window_prob(W, refined(colcounts(L)[1]))
            row.append('%d バイト %s' % (L, '標本では 0' if p == 0 else '2^%.1f' % math.log2(p)))
        print('  %s: ' % name + '、'.join(row))
    print('  例: ' + form(34, rng, False).decode('cp932') + ' / ' + form(38, rng, True).decode('cp932'))

def direct6(N=200000000, top=11000, seed=1):
    """確かめ: 2 文字の語 6 個を「・」でつないだ 34 バイトの鍵で、σ が範囲に入る割合を、列の独立を仮定せずに直接数える"""
    rng = np.random.default_rng(seed)
    w2 = [(b, f) for b, f in ja_words(top, 0) if len(b) == 4]
    cum = np.cumsum([f for b, f in w2]); cum /= cum[-1]
    contrib = np.zeros((6, len(w2)), dtype=np.int64)
    for k in range(6):
        for i, (b, f) in enumerate(w2):
            contrib[k, i] = sum(b[j] << (8 * ((6 * k + j) & 3)) for j in range(4))
    SEP = '・'.encode('cp932'); const = sum(SEP[j] << (8 * ((6 * k + 4 + j) & 3)) for k in range(5) for j in range(2))
    lo, hi = refined(9); hit = 0; done = 0
    while done < N:
        m = min(5000000, N - done)
        s = np.full(m, const, dtype=np.int64)
        for k in range(6): s += contrib[k][np.searchsorted(cum, rng.random(m))]
        hit += int(((s >= lo) & (s <= hi)).sum()); done += m
    return hit, N, len(w2)

if __name__ == '__main__' and len(sys.argv) > 2 and sys.argv[2] == 'direct':
    hit, n, nw = direct6(int(sys.argv[3]) if len(sys.argv) > 3 else 200000000)
    print('確かめ（直接数える）: 2 文字の語（上位 %d、頻度に比例）6 個を「・」でつないだ 34 バイト、%d 本のうち σ が範囲に入るもの %d 本 = 2^%.1f' % (nw, n, hit, math.log2(max(hit, 0.5) / n)))
