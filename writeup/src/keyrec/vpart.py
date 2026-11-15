# 「鍵の一部が分かっている」形の仮説を調べる入口。分かっていない部分は全バイト値で解く（vres を呼ぶ）。
#   仮説 = バイト列の並び。分かっている部分は bytes、分かっていない部分はバイト数（整数）で書く。
#   例: [b'GoodLuckToIISIGoodLuckTo', 12]   先頭 24 バイトが既知、残り 12 バイトが未知
# 解ける形: 未知のバイトを含む塊（4 バイトそろったもの）が 3 個まで。位置は塊の境目に合っていなくてよい。
#   それ以外の未知のバイト（4 個目以降の塊の中・端数の塊）は外側で全部回す（1 バイトにつき 221 倍）。
#   目安（2 コア）: 塊 3 個で 0.01 秒（格子 16）・0.09 秒（格子 17）、＋1 バイトで 1〜3 秒・7〜20 秒、＋端数 2 バイトで 3 分・25 分。
#                  塊 2 個＋2 バイト 0.1〜0.4 秒、＋3 バイト 30 秒〜1 分。塊 1 個＋3 バイト 15 秒。塊 0 個＋4 バイト 10 分。
import sys, os, subprocess, struct, time, random, math, itertools
from gks import key_schedule, chunks_of, chunk_vals, map_of_i2, h2d, C01b
from maps import setlow, ksword
HERE = os.path.dirname(os.path.abspath(__file__))

# ---- 鍵 v の目標（../TARGETS.json と、docs/password.md §1.2 のセル）
_rep = h2d('3fe99b19be354000'); _b0 = h2d('3fcf8c8bdb300652')
TARGET_V = dict(name='v', S=0x3a9d3481, Tlo=0xab0f, Thi=0xabc5, m=1,
                c3=_b0 + (2.906e-11 - 2.915e-11) / 2, h3=(2.906e-11 + 2.915e-11) / 2,
                cells1=[(_rep + (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2), (1.0 - _rep - (5.457e-12 - 1.819e-12) / 2, (5.457e-12 + 1.819e-12) / 2)],
                rep=(_rep, _b0), warm=23, blk2=0xed)

def rep5(x, byte):
    q = bytearray(struct.pack('<d', x)); q[5] = byte; return struct.unpack('<d', bytes(q))[0]
def ks_map1(k780, k7b0, b, warm, n):
    y = setlow(2.0 * k780 - 1.0, b); p = setlow(k7b0 * C01b + 1.9, b); out = []
    for i in range(warm + n):
        y = setlow(y * y - p, b)
        if i >= warm: out.append(ksword((1.0 + y) / 2.0))
    return out
def full_v(k780, k7b0, b, warm=23, blk2=0xed): return ks_map1(k780, k7b0, b, warm, 256) + ks_map1(k780, rep5(k7b0, blk2), b, 0, 79)
_REFV = None
def keystream_ok_v(key):
    """鍵 v の最終判定: その鍵の鍵スケジュールから作った 335 語（1 ブロック目 256 語＋E の 2 ブロック目 79 語）が、復元済みの鍵ストリームと全部同じか"""
    global _REFV
    if _REFV is None: _REFV = full_v(_rep, _b0, 0xab6a)
    r = key_schedule(key)
    if r['m'] != 1 or r['S'] != TARGET_V['S']: return False
    return full_v(r['k'][0], r['k'][2], r['b0b1'], r['warm']) == _REFV

# ---- 仮説 → 塊の並び
# 1 割り当てあたりの秒（2 コア、実測）。格子で解く塊（U）の数ごと
COST = {16: {0: 0.25e-6, 1: 1.5e-6, 2: 3e-6, 3: 0.012}, 17: {0: 0.25e-6, 1: 1.5e-6, 2: 7.5e-6, 3: 0.090}}
NB = 0xFC - 0x20 + 1
def segs_to_known(segs):
    known = []
    for s in segs: known += list(s) if isinstance(s, (bytes, bytearray)) else [None] * s
    return known
def build_hyp(hid, segs, lat=16, budget=1e18):
    """仮説を塊の並びに直す。未知のバイトを含む塊（4 バイトそろったもの）のうち 3 個までを格子で解く塊（U）にし、
    残りの未知のバイトは外側で回す（M）。U にした塊の中の既知のバイトは、出てきた鍵を後で照らして確かめる。
    費用が最小になる選び方を取る。lat = 16 では、U が 2 個以上のとき、U の 4 バイト目が 0x33 以上でなければならない
    （既知の 4 バイト目が 0x20〜0x32 の塊を U にする選び方は lat = 17 でだけ使う）。"""
    known = segs_to_known(segs); L = len(known); nch = (L + 3) // 4
    info = []
    for i in range(nch):
        ch = known[4 * i:4 * i + 4]; info.append((sum(b is None for b in ch), len(ch) == 4, ch))
    touched = [i for i in range(nch) if info[i][0] and info[i][1]]
    tail_unk = info[-1][0] if not info[-1][1] else 0
    best = None
    import itertools
    for r in range(0, min(3, len(touched)) + 1):
        for us in itertools.combinations(touched, r):
            if r >= 2:
                bad = False
                for i in us:
                    b = info[i][2][3]
                    if b is not None and (b < 0x20 or (lat == 16 and b < 0x33)): bad = True
                if bad: continue
            E = sum(info[i][0] for i in touched if i not in us) + tail_unk
            na = NB ** E; cost = (na / 3.0 if tail_unk else na) * COST[lat][r]
            if best is None or cost < best[0]: best = (cost, us, E)
    if best is None: return dict(id=hid, feasible=False, L=L, known=known, segs=segs, nU=len(touched), E=tail_unk, cost=float('inf'), touched=len(touched), tail_unk=tail_unk)
    cost, us, E = best; lines = []
    for i in range(nch):
        nunk, full, ch = info[i]
        if i in us: lines.append('U')
        elif nunk == 0: lines.append('K %x %d' % (sum(b << (8 * j) for j, b in enumerate(ch)), len(ch)))
        else: lines.append('M %x %x %d' % (sum((b or 0) << (8 * j) for j, b in enumerate(ch)), sum(1 << j for j, b in enumerate(ch) if b is None), len(ch)))
    return dict(id=hid, lines=lines, nU=len(us), E=E, L=L, known=known, segs=segs, cost=cost, feasible=cost <= budget, touched=len(touched), tail_unk=tail_unk, lat=lat)
def matches(h, key):
    return len(key) == h['L'] and all(k is None or k == b for k, b in zip(h['known'], key))
_OKASC = set(b"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 .,-_'!?@&/:()+*#")
def textlike(key):
    """鍵が文字列として読めるか。0 = 読めない（Shift_JIS として不正、制御文字、外字）、1 = Shift_JIS の文字列としては正しい、
    2 = 全部がふつうに打てる文字（半角の印字可能な文字・半角カナ・全角の記号/英数/かな・JIS 第 1・第 2 水準の漢字）、
    3 = さらに絞った形（半角は英数字・空白とよく使う記号だけ、全角はよく使う記号・英数・かな・第 1 水準の漢字だけ。半角カナなし）"""
    i = 0; n = len(key); lvl = 3
    while i < n:
        b = key[i]
        if 0x20 <= b <= 0x7E:
            if b not in _OKASC: lvl = min(lvl, 2)
            i += 1; continue
        if 0xA1 <= b <= 0xDF: lvl = min(lvl, 2); i += 1; continue
        if (0x81 <= b <= 0x9F or 0xE0 <= b <= 0xFC) and i + 1 < n:
            t = key[i + 1]
            if not (0x40 <= t <= 0x7E or 0x80 <= t <= 0xFC): return 0
            try: ch = key[i:i + 2].decode('cp932')
            except UnicodeDecodeError: return 0
            if len(ch) != 1 or 0xE000 <= ord(ch) <= 0xF8FF: return 0
            if not (0x81 <= b <= 0x84 or 0x88 <= b <= 0x9F or 0xE0 <= b <= 0xEA): lvl = min(lvl, 1)
            elif not ((b == 0x81 and t <= 0x97) or b == 0x82 or (b == 0x83 and t <= 0x96) or (0x88 <= b <= 0x97 and not (b == 0x88 and t < 0x9F)) or (b == 0x98 and t <= 0x72)): lvl = min(lvl, 2)
            i += 2; continue
        return 0
    return lvl
def iter_logs(logdir='logs'):
    """vcamp・vone の記録（logs/vcamp_*.jsonl と、その .gz）を 1 行ずつ返す"""
    import gzip, json
    for fn in sorted(os.listdir(logdir)):
        if not fn.startswith('vcamp_'): continue
        if fn.endswith('.jsonl'): f = open(os.path.join(logdir, fn), encoding='utf-8')
        elif fn.endswith('.jsonl.gz'):
            if os.path.exists(os.path.join(logdir, fn[:-3])): continue          # 同じ名前の .jsonl があればそちらを読む
            f = gzip.open(os.path.join(logdir, fn), 'rt', encoding='utf-8')
        else: continue
        with f:
            for line in f:
                try: yield json.loads(line)
                except Exception: continue
def load_glog(logdir='logs'):
    """既知のバイトの並び → 記録（同じ形を 2 回解いたときは格子の大きいほう）"""
    glog = {}
    for r in iter_logs(logdir):
        if r['pat'] not in glog or r['lat'] >= glog[r['pat']]['lat']: glog[r['pat']] = r
    return glog
def show(h):
    out = ''
    for s in h['segs']: out += ('[%s]' % s.decode('cp932', 'replace')) if isinstance(s, (bytes, bytearray)) else '<%d>' % s
    return out

def run_vres(tgt, hyps, L=16, blo=0x20, bhi=0xFC, tag='vp', log=None, keep=False):
    """hyps: build_hyp の結果のリスト。戻り値: (鍵の一覧 {id: [bytes..]}, 統計 {id: dict}, 秒)"""
    nmax = max(len(h['lines']) for h in hyps); sjmax = 4 * (2 ** nmax - 1)
    path = os.path.join(HERE, 'tmp_%s.in' % tag)
    with open(path, 'w') as f:
        f.write('T %x %x %x %d %.17g %.17g %d %s %d' % (tgt['S'], tgt['Tlo'], tgt['Thi'], tgt['m'], tgt['c3'], tgt['h3'], len(tgt['cells1']),
                                                         ' '.join('%.17g %.17g' % c for c in tgt['cells1']), 1 if tgt.get('k788') else 0))
        if tgt.get('k788'): f.write(' %.17g %.17g %.17g %.17g' % (tgt['k788'][0], tgt['k788'][1], tgt['k7b8'][0], tgt['k7b8'][1]))
        f.write('\nO %d %d %d %d\n' % (L, blo, bhi, sjmax))
        for h in hyps: f.write('H %s %d\n%s\n' % (h['id'], len(h['lines']), '\n'.join(h['lines'])))
    t0 = time.time(); r = subprocess.run([os.path.join(HERE, 'vres'), path], capture_output=True, text=True); dt = time.time() - t0
    if r.returncode: raise RuntimeError('vres failed: %s' % r.stderr[-500:])
    keys = {}; stat = {}; raw = {}; byid = {h['id']: h for h in hyps}
    for line in r.stdout.split('\n'):
        p = line.split()
        if not p: continue
        if p[0] == 'KEY':
            kb = bytes.fromhex(p[3]); raw[p[1]] = raw.get(p[1], 0) + 1
            if matches(byid[p[1]], kb): keys.setdefault(p[1], set()).add(kb)
        elif p[0] == 'STAT': stat[p[1]] = dict(zip(p[2::2], [int(x) for x in p[3::2]]))
    for i, n in raw.items(): stat.setdefault(i, {})['rawkeys'] = n
    if not keep: os.remove(path)
    return keys, stat, dt

# ---- 対照: 乱数の鍵で目標を作る
def synth_target(key, rng, tw=183):
    r = key_schedule(key); off = rng.randrange(tw); h1 = (5.457e-12 + 1.819e-12) / 2; h3 = (2.906e-11 + 2.915e-11) / 2
    c1 = r['k'][0] + (rng.random() - 0.5) * 2 * h1 * 0.9; c3 = r['k'][2] + (rng.random() - 0.5) * 2 * h3 * 0.9
    return dict(name='synthetic', S=r['S'], Tlo=max(0, r['b0b1'] - off), Thi=min(0xFFFF, r['b0b1'] - off + tw - 1), m=r['m'], c3=c3, h3=h3,
                cells1=[(c1, h1), (1.0 - c1, h1)]), r

if __name__ == '__main__':
    print(__doc__ or 'vpart: see the comments at the top')
