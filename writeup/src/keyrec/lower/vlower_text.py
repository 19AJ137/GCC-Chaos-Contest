#!/usr/bin/env python3
# 英小文字だけの鍵・英文の形の鍵が、鍵 v の S・b0b1 の条件（塊の値の和 σ の範囲）に合う割合の見積もり（docs/password.md §8.9）
#   python3 lower/vlower_text.py [標本の数（文の型）] [標本の数（S の割合）]
# 1) 文字が独立に出るとみた型（a〜z の一様、英語の文字の頻度、ローマ字の文字の頻度）
# 2) 英単語を頻度に比例して並べた文（空白なし / 空白あり / 語頭だけ大文字で空白なし）。列ごとの和の分布を標本から作り、列どうしは独立とみる
# 3) σ が範囲に入る英小文字だけの鍵のうち、S・b0b1・写像まで合う割合（鍵スケジュールを回して数える）
import sys, os, math, random, re
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from vlower import refined, colcounts, ways, count_window, pm_of, LOW, T
from gks import key_schedule

def window_prob(W, win):
    """W[j] = 列 j の和の分布。σ = Σ 256^j·B_j が win に入る確率（列どうしは独立とみる）"""
    lo, hi = win
    P0 = np.concatenate([[0.0], np.cumsum(W[0])])
    def r0(a, b):
        a = max(a, 0); b = min(b, len(W[0]) - 1)
        return P0[b + 1] - P0[a] if b >= a else 0.0
    tot = 0.0
    for B3 in np.nonzero(W[3])[0]:
        base3 = int(B3) << 24
        if base3 > hi: break
        b2lo = max(0, (lo - base3 - 256 * (len(W[1]) - 1) - (len(W[0]) - 1)) // 65536)
        for B2 in range(b2lo, min((hi - base3) // 65536, len(W[2]) - 1) + 1):
            if W[2][B2] == 0: continue
            for B1 in np.nonzero(W[1])[0]:
                base1 = base3 + (B2 << 16) + (int(B1) << 8)
                if base1 > hi: break
                v = r0(lo - base1, hi - base1)
                if v > 0: tot += W[3][B3] * W[2][B2] * W[1][B1] * v
    return tot

def letter_models():
    import wordfreq
    en = np.zeros(256); words = []
    fe = wordfreq.get_frequency_dict('en', 'large')
    for w, f in sorted(fe.items(), key=lambda x: -x[1])[:30000]:
        if re.fullmatch('[a-z]+', w):
            words.append((w, f))
            for ch in w: en[ord(ch)] += f
    en /= en.sum()
    # ローマ字: 日本語の頻度表の語に出るかな 1 字ずつを、ヘボン式に近い表で置き換えて数える（漢字は読みが分からないので数えない。粗い型）
    base = dict(zip('あいうえおかきくけこさしすせそたちつてとなにぬねのはひふへほまみむめもやゆよらりるれろわをんがぎぐげござじずぜぞだぢづでどばびぶべぼぱぴぷぺぽゃゅょっー',
                    'a i u e o ka ki ku ke ko sa shi su se so ta chi tsu te to na ni nu ne no ha hi fu he ho ma mi mu me mo ya yu yo ra ri ru re ro wa o n ga gi gu ge go za ji zu ze zo da ji zu de do ba bi bu be bo pa pi pu pe po ya yu yo t u'.split()))
    ro = np.zeros(256)
    for w, f in wordfreq.get_frequency_dict('ja', 'large').items():
        for ch in w:
            o = ord(ch)
            if 0x30A1 <= o <= 0x30F6: ch = chr(o - 0x60)
            for c in base.get(ch, ''): ro[ord(c)] += f
    ro /= ro.sum()
    uni = pm_of(LOW) / 26.0
    return uni, en, ro, words

def col_pmfs(sampler, L, N, rng):
    """文の型の標本から、列ごとの和の分布を作る"""
    cnt, n = colcounts(L)
    W = [np.zeros(127 * cnt[j] + 1) for j in range(4)]
    for _ in range(N):
        s = sampler(L, rng)
        b = [0, 0, 0, 0]
        for i, ch in enumerate(s): b[i & 3] += ch
        for j in range(4): W[j][b[j]] += 1
    return [w / N for w in W]

def make_samplers(words):
    ws = [w for w, f in words]; cum = np.cumsum([f for w, f in words]); cum /= cum[-1]
    def pick(rng): return ws[int(np.searchsorted(cum, rng.random()))]
    def text(L, rng, sep, cap):
        while True:
            out = []; n = 0
            while n < L:
                w = pick(rng)
                if cap: w = w[0].upper() + w[1:]
                if out and sep: n += 1
                out.append(w); n += len(w)
            if n == L: return (' ' if sep else '').join(out).encode()
    return (lambda L, rng: text(L, rng, False, False)), (lambda L, rng: text(L, rng, True, False)), (lambda L, rng: text(L, rng, False, True))

def sample_feasible(L, rng, tabs={}):
    """σ が範囲に入る英小文字だけの鍵を、一様に 1 本作る"""
    if L not in tabs:
        cnt, n = colcounts(L); pm = pm_of(LOW)
        t, parts, cnt, n = count_window(pm, L, refined(n))
        ks = list(parts); pw = np.cumsum([parts[k] for k in ks]); pw /= pw[-1]
        wt = {m: [ways(pm, k) for k in range(m + 1)] for m in set(cnt)}
        tabs[L] = (cnt, n, ks, pw, wt, refined(n))
    cnt, n, ks, pw, wt, (lo, hi) = tabs[L]
    B3, B2, B1 = ks[int(np.searchsorted(pw, rng.random()))]
    base = (B3 << 24) + (B2 << 16) + (B1 << 8)
    w0 = wt[cnt[0]][cnt[0]]; a = max(lo - base, 0); b = min(hi - base, len(w0) - 1)
    c0 = np.cumsum(w0[a:b + 1]); B0 = a + int(np.searchsorted(c0, rng.random() * c0[-1]))
    cols = []
    for j, B in enumerate((B0, B1, B2, B3)):
        m = cnt[j]; rem = B; col = []
        for k in range(m, 0, -1):
            w = np.array([wt[m][k - 1][rem - c] if 0 <= rem - c < len(wt[m][k - 1]) else 0.0 for c in LOW])
            c = LOW[int(np.searchsorted(np.cumsum(w), rng.random() * w.sum()))]
            col.append(c); rem -= c
        rng.shuffle(col); cols.append(col)
    return bytes(cols[i & 3][i >> 2] for i in range(L))

if __name__ == '__main__':
    N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
    N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
    rng = random.Random(20261115); nrng = np.random.default_rng(20261115)
    uni, en, ro, words = letter_models()
    asc = pm_of(range(32, 127)) / 95.0
    print('== 文字の平均（文字コード）: 一様 %.1f、英語 %.1f、ローマ字 %.1f。語の数（英語、a〜z だけ）%d' % (sum(i * uni[i] for i in range(256)), sum(i * en[i] for i in range(256)), sum(i * ro[i] for i in range(256)), len(words)))
    print('== 文字が独立に出る型: σ が範囲に入る割合（log2）と、2 文字目の列（10 字）の和が 1185・1186 になる割合')
    for L in (39, 40):
        cnt, n = colcounts(L); win = refined(n)
        for name, pm in (('a〜z の一様', uni), ('英語の文字の頻度', en), ('ローマ字の文字の頻度', ro)):
            W = [ways(pm, m) for m in cnt]
            p = window_prob(W, win); w1 = W[1]
            print('  %d 文字 %-12s: 2^%.1f。2 文字目の列だけで %.2e（%.0f 回に 1 回）、4 文字目の列（和 984）だけで %.2e' % (L, name, math.log2(p), w1[1185] + w1[1186], 1 / (w1[1185] + w1[1186]), W[3][984]))
    print('  参考（印字可能な ASCII の一様、36〜40 文字）: ' + '、'.join('%d 文字 2^%.1f' % (L, math.log2(window_prob([ways(asc, m) for m in colcounts(L)[0]], refined(colcounts(L)[1])))) for L in range(36, 41)))
    cat, spc, cam = make_samplers(words)
    print('== 英単語を頻度に比例して並べた文（標本 %d 本ずつ。列ごとの和の分布を標本から作り、列どうしは独立とみた見積もり）: σ が範囲に入る割合（log2）' % N1)
    for name, smp in (('空白なし・小文字', cat), ('空白あり・小文字', spc), ('空白なし・語頭だけ大文字', cam)):
        row = []
        for L in range(36, 41):
            W = col_pmfs(smp, L, N1, nrng); p = window_prob(W, refined(colcounts(L)[1]))
            row.append('%d 文字 %s' % (L, 'なし（標本では 0）' if p == 0 else '2^%.1f' % math.log2(p)))
        print('  %-14s: ' % name + '、'.join(row))
    print('  例: ' + ' / '.join(s(38, nrng).decode() for s in (cat, spc, cam)))
    print('== σ が範囲に入る英小文字だけの鍵（一様に作った標本 %d 本ずつ）のうち、S・b0b1・写像まで合う割合' % N2)
    for L in (39, 40):
        nS = nT = nM = 0; ex = []
        for i in range(N2):
            k = sample_feasible(L, rng)
            if i < 3: ex.append(k.decode())
            r = key_schedule(k)
            if r['S'] == T['S']:
                nS += 1
                if T['Tlo'] <= r['b0b1'] <= T['Thi']:
                    nT += 1
                    if r['m'] == 1: nM += 1
        cntw, _, _, _ = count_window(pm_of(LOW), L, refined(colcounts(L)[1]))
        est = math.log2(cntw) + math.log2(max(nM, 0.5) / N2) + math.log2(2 * 7.3e-12 * 5.8e-11)
        print('  %d 文字: S が合う %d 本（%.2f%%）、b0b1 も合う %d 本、写像 1 も合う %d 本（%.3f%%）。等価な鍵の数の見積もり 2^%.1f × %.2e × k780・k7b0 のセル 2^%.1f = 2^%.0f' % (
            L, nS, 100 * nS / N2, nT, nM, 100 * nM / N2, math.log2(cntw), nM / N2, math.log2(2 * 7.3e-12 * 5.8e-11), est))
        print('     σ が範囲に入る鍵の例（等価な鍵ではない）: ' + ' / '.join(ex))
