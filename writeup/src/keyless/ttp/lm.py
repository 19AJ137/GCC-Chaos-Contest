"""lm.py — 文字 n-gram の言語モデル（英文。stupid backoff）。2 本の文章の XOR をほどく探索（ttp.py）用。2026-11-12
  コーパスは一般の英文（解読のときは、Debian 系の環境にある vim の説明書・ライセンス文と、Newton の Opticks を使った）。行末は CR LF に直す。
  別の英文を使うときは、環境変数 LM_CORPUS にファイルの glob パターンを指定する（複数は os.pathsep で区切る）。"""
import os, glob, math, pickle, random, re
ORDER = 5
ALPHA = [0x0d, 0x0a] + list(range(0x20, 0x7f))
AIDX = {c: i for i, c in enumerate(ALPHA)}
def corpus(maxbytes=6_000_000, seed=1):
    pat = os.environ.get('LM_CORPUS')
    if pat: files = sorted(f for p in pat.split(os.pathsep) for f in glob.glob(p))
    else: files = sorted(glob.glob('/usr/share/vim/vim91/doc/*.txt')) + sorted(glob.glob('/usr/share/common-licenses/*')) + ['/usr/local/go1.25.1/src/testdata/Isaac.Newton-Opticks.txt']
    rnd = random.Random(seed); out = []; tot = 0
    lic = [f for f in files if 'common-licenses' in f]; other = [f for f in files if 'common-licenses' not in f]
    rnd.shuffle(other)
    for f in lic + ['/usr/local/go1.25.1/src/testdata/Isaac.Newton-Opticks.txt'] + other:
        try: t = open(f, 'rb').read()
        except Exception: continue
        if os.path.islink(f): continue
        t = t.replace(b'\r\n', b'\n').replace(b'\t', b'    ')
        t = bytes(c for c in t if c == 0x0a or 0x20 <= c < 0x7f)
        t = re.sub(rb'[*|]([A-Za-z0-9_:.\-]+)[*|]', rb'\1', t)          # vim の説明書の印を外す
        t = re.sub(rb' {3,}', b'  ', t)
        t = t.replace(b'\n', b'\r\n')
        out.append(t); tot += len(t)
        if tot > maxbytes: break
    return out
class LM:
    def __init__(self, texts=None, order=ORDER):
        self.order = order; self.c = [dict() for _ in range(order + 1)]; self.cache = {}
        if texts:
            for t in texts: self.add(t)
    def add(self, t):
        n = len(t)
        for k in range(1, self.order + 1):
            d = self.c[k]
            for i in range(n - k + 1):
                g = t[i:i + k]; d[g] = d.get(g, 0) + 1
        self.c[0][b''] = self.c[0].get(b'', 0) + n
    def finalize(self):
        # 文脈の数え上げ: c[k-1] は「その文脈に続きがある回数」に近いが、末尾の差は無視する
        self.V = len(ALPHA)
    def logp_vec(self, ctx):
        """文脈 ctx（bytes、長さ order-1 まで）に続く各文字の log 確率（ALPHA の順）"""
        v = self.cache.get(ctx)
        if v is not None: return v
        out = []
        for ch in ALPHA:
            cb = bytes([ch]); bo = 1.0; p = None; cx = ctx
            while True:
                k = len(cx) + 1
                num = self.c[k].get(cx + cb, 0)
                if num:
                    den = self.c[k - 1].get(cx, 0) if k > 1 else self.c[0][b'']
                    p = bo * num / den; break
                if not cx: p = bo * 0.05 / self.c[0][b'']; break
                cx = cx[1:]; bo *= 0.4
            out.append(math.log(p))
        if len(self.cache) < 400000: self.cache[ctx] = out
        return out
    def score(self, t):
        s = 0.0
        for i in range(len(t)):
            ctx = t[max(0, i - self.order + 1):i]
            s += self.logp_vec(ctx)[AIDX[t[i]]]
        return s
def load(path='lm.pkl'):
    if os.path.exists(path):
        lm = LM(); lm.order, lm.c = pickle.load(open(path, 'rb')); lm.finalize(); return lm
    lm = LM(corpus()); lm.finalize(); pickle.dump((lm.order, lm.c), open(path, 'wb'), protocol=4); return lm
if __name__ == '__main__':
    import time; t0 = time.time(); lm = load(); print('n-gram の種類: ' + ', '.join('%d: %d' % (k, len(lm.c[k])) for k in range(1, lm.order + 1)), '%.0f 秒' % (time.time() - t0))
    for s in (b'the quick brown fox jumps over the lazy dog', b'xqzj kkvw pqqz mmnb', b'Please send the plaintext and the key to us.\r\nThe prize is one million yen.'):
        print('  %.2f ビット/文字  %r' % (-lm.score(s) / len(s) / math.log(2), s))
