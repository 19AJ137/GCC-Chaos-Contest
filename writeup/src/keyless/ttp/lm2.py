"""lm2.py — バイト単位の n-gram 言語モデル（英文と Shift_JIS の日本語）。fulldec2.py 用。2026-11-12
  E 用（英文だけ）と A 用（英文＋日本語）の 2 つを作る。コーパス:
    英文（一般）: lm.py の corpus()（環境変数 LM_CORPUS）
    英文（IISI）: dom_all.txt（IISI の英語ページ GCC-over・prize-e・infoguard-e・chaosmail-e の本文を "\n=====\n" で区切って並べたもの）
    日本語（IISI）: domj/ の中のファイル（PRIZE.TXT・README.TXT・INFNOTE.TXT）と、IISI_PAGES の下の HTML（GCC-gaiyou.htm、siyourei.htm、infgd40.html）
    日本語（一般）: genj/ の中のファイル（雑誌の付録 CD-ROM の README など）
  IISI のページや配布物の文章は著作物なので、このリポジトリには入れていない。各自で用意して置く（Wayback Machine の iisi.co.jp、software/ の配布物）。
  環境変数 IISI_PAGES は、HTML を置いたフォルダ（中に webpage/GCC-gaiyou.htm などの並びがあるもの）"""
import os, glob, math, pickle, re, html, sys
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
import lm as L
ORDER = 5
A7 = [0x0d, 0x0a] + list(range(0x20, 0x7f))
A7IDX = {c: i for i, c in enumerate(A7)}
U = os.environ.get('IISI_PAGES', os.path.join(D, 'iisi')) + os.sep
def norm(t):
    t = t.replace(b'\r\n', b'\n').replace(b'\t', b' ').replace(b'\x1a', b'')
    t = bytes(c for c in t if c == 0x0a or c >= 0x20)
    return t.replace(b'\n', b'\r\n')
def strip_html(p):
    t = open(p, 'rb').read().decode('cp932', 'replace')
    t = re.sub(r'(?is)<(script|style).*?</\1>', ' ', t); t = re.sub(r'(?s)<[^>]+>', ' ', t); t = html.unescape(t)
    t = re.sub(r'[ \t]+', ' ', t); t = re.sub(r'\n\s*\n+', '\n', t)
    return t.encode('cp932', 'ignore')
def jp_domain():
    out = [norm(open(f, 'rb').read()) for f in sorted(glob.glob(os.path.join(D, 'domj', '*')))]
    for p in ('webpage/GCC-gaiyou.htm', 'webpage/wayback_1997/www.iisi.co.jp/reserch/siyourei.htm', 'webpage/wayback_dlsites/www.iisi.co.jp/vcs/software/win95/infgd40.html'):
        t = strip_html(U + p)
        i = t.find('ＧＣＣカオス暗号の概要'.encode('cp932'))
        if i > 0: t = t[i:]
        out.append(norm(t)); out.append(norm(t.replace(b'\n', b'')))
    return out
def jp_general(): return [norm(open(f, 'rb').read()) for f in sorted(glob.glob(os.path.join(D, 'genj', '*')))]
def en_domain(): return [t for t in open(os.path.join(D, 'dom_all.txt'), 'rb').read().split(b'\n=====\n')]
class LM2:
    def __init__(self, order=ORDER):
        self.order = order; self.c = [dict() for _ in range(order + 1)]; self.c[0][b''] = 0; self.v7 = {}; self.m1 = {}
    def add(self, t, w=1):
        n = len(t)
        for k in range(1, self.order + 1):
            d = self.c[k]
            for i in range(n - k + 1):
                g = t[i:i + k]; d[g] = d.get(g, 0) + w
        self.c[0][b''] += n * w
    def p1(self, ctx, ch):
        cb = bytes([ch]); bo = 1.0; cx = ctx; c = self.c
        while True:
            k = len(cx) + 1; num = c[k].get(cx + cb, 0)
            if num:
                den = c[k - 1].get(cx, 0) if k > 1 else c[0][b'']
                return math.log(bo * num / den)
            if not cx: return math.log(bo * 0.02 / c[0][b''])
            cx = cx[1:]; bo *= 0.4
    def logp1(self, ctx, ch):
        key = (ctx, ch); v = self.m1.get(key)
        if v is None:
            v = self.p1(ctx, ch)
            if len(self.m1) < 3000000: self.m1[key] = v
        return v
    def logp_vec7(self, ctx):
        v = self.v7.get(ctx)
        if v is None:
            v = [self.p1(ctx, ch) for ch in A7]
            if len(self.v7) < 300000: self.v7[ctx] = v
        return v
    def score(self, t):
        return sum(self.logp1(t[max(0, i - self.order + 1):i], t[i]) for i in range(len(t)))
def build():
    gen = L.corpus(); ed = en_domain(); jd = jp_domain(); jg = jp_general()
    print('英文（一般）%d B、英文（IISI）%d B、日本語（IISI）%d B、日本語（一般）%d B' % (sum(map(len, gen)), sum(map(len, ed)), sum(map(len, jd)), sum(map(len, jg))))
    le = LM2()
    for t in gen: le.add(t)
    for t in ed: le.add(t, 20)
    la = LM2()
    for t in gen: la.add(t)
    for t in ed: la.add(t, 5)
    for t in jd: la.add(t, 20)
    for t in jg: la.add(t, 4)
    pickle.dump((le.order, le.c), open(os.path.join(D, 'lm2_e.pkl'), 'wb'), protocol=4); pickle.dump((la.order, la.c), open(os.path.join(D, 'lm2_a.pkl'), 'wb'), protocol=4)
    return le, la
def load(name):
    m = LM2(); m.order, m.c = pickle.load(open(os.path.join(D, name), 'rb')); return m
if __name__ == '__main__':
    le, la = build()
    for s in ('カオス暗号は安全である。'.encode('cp932'), 'ＧＣＣカオス暗号方式'.encode('cp932'), b'You can not decrypt GCC7!', bytes([0x83, 0x99, 0xe1, 0x41, 0x92, 0x7f])):
        print('  A 用 %.2f ビット/バイト、E 用 %.2f  %r' % (-la.score(s) / len(s) / math.log(2), -le.score(s) / len(s) / math.log(2), s))
