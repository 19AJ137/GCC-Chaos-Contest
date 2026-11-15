# 鍵 v が「鍵 k と同じ形式」だという仮説を、語を一般の頻度辞書から取って調べる（11-14、3 回目）。
#   形: 2〜3 文字の語（全角）n 個と IISI を、区切り文字でつないだもの（IISI の位置は順不同）。
#   語の一覧: wordfreq（pip。日本語の頻度表 large）の、漢字だけ・カタカナだけの 2〜3 文字の語を頻度順に取る。
#   調べ方: kmitm（中間一致）。語の枠を 左 2・右 2（・外側 1）に分ける。出てきた候補は鍵スケジュール全体と鍵ストリームで判定する。
#   python3 vdict.py words                         語の一覧の大きさと、磁石・発光の順位
#   python3 vdict.py control n [本数]              対照（辞書の語で作った鍵を埋めて、出るかを見る）
#   python3 vdict.py run n N2 N3 [区切り…]         本番（鍵 v の目標）。N2・N3 = 2 文字・3 文字の語の数。区切りは 'dot' 'none' 'sp' など（既定は dot）
#   python3 vdict.py count n N2 N3 [区切り…]       形の数と見積もり
#   等価な形（塊を並べ替えると同じ鍵の集まりになる形。family_sig）は、1 つだけ回して、ほかの形には記録を写す（記録の same_as が写した元）。
#   環境変数 VDICT_NODEDUP=1 で、全部の形を別々に回す（11-14 の 4 回目の本番 125 形は、まとめる前の作りで全部を回した）。
import sys, os, json, time, random, subprocess, itertools
from gks import key_schedule, chunk_vals, map_of_i2
from vpart import TARGET_V, keystream_ok_v, synth_target
HERE = os.path.dirname(os.path.abspath(__file__)); TMP = os.path.join(HERE, 'tmp_dict'); os.makedirs(TMP, exist_ok=True)
IISI = b'IISI'
SEPS = {'dot': '・', 'none': '', 'sp': ' ', 'hy': '-', 'us': '_', 'pd': '.', 'cm': ',', 'sl': '/', 'pl': '+', 'am': '&', 'at': '@', 'cl': ':', 'eq': '=', 'as': '*',
        'zsp': '　', 'ten': '、', 'zcm': '，', 'zpd': '．', 'zcl': '：', 'zsl': '／', 'zhy': '－', 'zus': '＿', 'zpl': '＋', 'zam': '＆', 'zeq': '＝', 'zas': '＊', 'bar': 'ー', 'x': '×', 'arrow': '→', 'star': '☆'}

def cls(w):
    if all('一' <= c <= '鿿' or c == '々' for c in w): return 'kanji'
    if all('ァ' <= c <= 'ー' for c in w): return 'kata'
    return None
_WORDS = None
def load_words():
    """{2: [(語, 種類)…頻度順], 3: […]}。cp932 で 1 文字 2 バイトになるものだけ"""
    global _WORDS
    if _WORDS is None:
        import wordfreq
        d = wordfreq.get_frequency_dict('ja', 'large'); out = {2: [], 3: []}
        for w, f in sorted(d.items(), key=lambda kv: (-kv[1], kv[0])):
            if len(w) not in (2, 3) or cls(w) is None: continue
            try: b = w.encode('cp932')
            except UnicodeEncodeError: continue
            if len(b) != 2 * len(w): continue
            out[len(w)].append((w, cls(w)))
        _WORDS = out
    return _WORDS
def word_bytes(nchars, N):
    return [w.encode('cp932') for w, c in load_words()[nchars][:N]]

def make_struct(sizes, pos, sep):
    """sizes: 語ごとのバイト数、pos: IISI の位置（0〜n）、sep: 区切りのバイト列。戻り値: (ひな形, [(off, len)…語の枠])"""
    elems = [('w', s) for s in sizes]; elems.insert(pos, ('k', IISI)); tpl = b''; slots = []
    for i, (t, v) in enumerate(elems):
        if i: tpl += sep
        if t == 'w': slots.append((len(tpl), v)); tpl += b'\x00' * v
        else: tpl += v
    return tpl, slots
PLAN_OUTER4 = False
def chunk_slots(L, slots):
    nch = (L + 3) // 4
    return [[k for k, (o, n) in enumerate(slots) if o < 4 * i + 4 and o + n > 4 * i] for i in range(nch)]
def plan_groups(L, slots, lists_n, prefer_outer=None):
    """枠を 左・右（各 1〜2 個）・外側（0〜1 個）に分ける。1 つの塊が 2 つの組にまたがらない分け方のうち、手間が最小のもの。なければ None"""
    n = len(slots); cs = chunk_slots(L, slots); share = set()
    for s in cs:
        for a in s:
            for b in s:
                if a < b: share.add((a, b))
    best = None
    # 語 3 個は、外側 1 個＋左右 1 個ずつの分け方も使う（辞書の全語を入れるため）。枠 4 個で外側を使うのは、vdict3.py が PLAN_OUTER4 を立てたときだけ
    outers = [None] + list(range(n)) if (n == 3 or (n == 4 and PLAN_OUTER4)) else ([None] if n <= 4 else list(range(n)))
    for o in outers:
        rest = [k for k in range(n) if k != o]
        if o is not None and any(o in p for p in share): continue
        if len(rest) < 2 or len(rest) > 4: continue
        for r in range(1, len(rest)):
            for left in itertools.combinations(rest, r):
                right = tuple(k for k in rest if k not in left)
                if len(left) > 2 or len(right) > 2 or (left and right and left[0] > right[0]): continue
                if any((a in left) != (b in left) for a, b in share): continue
                nl = 1; nr = 1
                for k in left: nl *= lists_n[k]
                for k in right: nr *= lists_n[k]
                no = lists_n[o] if o is not None else 1
                if max(nl, nr) > 1.3e8: continue
                cost = no * min(nl, nr) * 2.5e-9 + no * nl * nr * 5.7e-8 * 0.1 * 2e-9 + (nl + nr) * 8e-8      # 2 コア・v2 の目安（和が合う割合を 0.1 とみる）
                if best is None or cost < best[0]: best = (cost, left, right, o)
    return best

def last_slot_filter(tpl, slots, words_by_slot, want=1):
    """最後の塊に枠が 1 つだけかかるとき、その枠の語を「写像が want になるもの」に絞る。戻り値: (語の一覧のリスト, 写像が決まっていて合わないなら False)"""
    L = len(tpl); cs = chunk_slots(L, slots); last = cs[-1]; base = 4 * (len(cs) - 1)
    if not last:
        c = int.from_bytes(tpl[base:], 'little'); return words_by_slot, map_of_i2(chunk_vals(c)[1][1]) == want
    if len(last) > 1: return words_by_slot, True
    k = last[0]; o, n = slots[k]; keep = []
    for w in words_by_slot[k]:
        key = bytearray(tpl); key[o:o + n] = w; c = int.from_bytes(bytes(key[base:]), 'little')
        if map_of_i2(chunk_vals(c)[1][1]) == want: keep.append(w)
    out = list(words_by_slot); out[k] = keep; return out, True

def run_kmitm(tgt, tpl, slots, words_by_slot, plan, jwin, tag='km', threads=None):
    cost, left, right, o = plan; spec = os.path.join(TMP, 'spec_%s.txt' % tag); lines = []
    c1, h1 = tgt['cells1'][0]; c1b, h1b = tgt['cells1'][1]
    lines.append('T %x %.17g %.17g %.17g %.17g 2e-12' % (tgt['S'], c1, h1, c1b, h1b)); lines.append('J %d %d' % jwin); lines.append('K %d %s' % (len(tpl), tpl.hex()))
    for k, (off, n) in enumerate(slots):
        p = os.path.join(TMP, 'wl_%s_%d.bin' % (tag, k))
        with open(p, 'wb') as f: f.write(b''.join(words_by_slot[k]))
        lines.append('W %d %d %d %s' % (k, n, len(words_by_slot[k]), p))
    for g, grp in ((0, left), (1, right), (2, () if o is None else (o,))):
        lines.append('G %d %d %s' % (g, len(grp), ' '.join('%d %d' % (slots[k][0], k) for k in grp)))
    lines.append('R')
    with open(spec, 'w') as f: f.write('\n'.join(lines) + '\n')
    env = dict(os.environ)
    if threads: env['OMP_NUM_THREADS'] = str(threads)
    t0 = time.time(); r = subprocess.run([os.path.join(HERE, os.environ.get('KMITM_EXE', 'kmitm')), spec], capture_output=True, text=True, env=env); dt = time.time() - t0
    if r.returncode: raise RuntimeError('kmitm failed: ' + r.stderr[-300:])
    cands = []; stat = {}
    for line in r.stdout.split('\n'):
        p = line.split()
        if not p: continue
        if p[0] == 'CAND': cands.append(bytes.fromhex(p[1]))
        elif p[0] == 'STAT': stat = dict(zip(p[1::2], p[2::2]))
    for k in range(len(slots)):
        try: os.remove(os.path.join(TMP, 'wl_%s_%d.bin' % (tag, k)))
        except OSError: pass
    return cands, stat, dt

def sched_ok(tgt, key):
    """鍵スケジュール全体（写像・S・b0b1 の範囲・k780・k7b0 のセル）"""
    r = key_schedule(key)
    if r['m'] != tgt['m'] or r['S'] != tgt['S'] or not (tgt['Tlo'] <= r['b0b1'] <= tgt['Thi']): return False
    if not any(abs(r['k'][0] - c) <= h for c, h in tgt['cells1']): return False
    return abs(r['k'][2] - tgt['c3']) <= tgt['h3']

def struct_list(n, seps, maxlen=40):
    for sizes in itertools.product((4, 6), repeat=n):
        for pos in range(n + 1):
            for sn in seps:
                sep = SEPS[sn].encode('cp932'); tpl, slots = make_struct(sizes, pos, sep)
                if len(tpl) > maxlen: continue
                yield sizes, pos, sn, tpl, slots
def show_struct(sizes, pos, sn):
    el = ['語%d' % (s // 2) for s in sizes]; el.insert(pos, 'IISI'); return (SEPS[sn] or '').join(el) if SEPS[sn] else ''.join(el)

def family_sig(tpl, slots):
    """塊を並べ替えると同じ鍵の集まりになる形どうしで等しくなる印（最後の塊は動かさない）。1 つの塊に語の枠が 2 つかかる形は None（まとめない）"""
    L = len(tpl); owner = {}; nch = (L + 3) // 4
    for k, (o, n) in enumerate(slots):
        for j in range(n): owner[o + j] = (k, j)
    per = {k: [] for k in range(len(slots))}; fixed = []; last = None
    for i in range(nch):
        d = []; ks = set()
        for b in range(4 * i, min(4 * i + 4, L)):
            if b in owner: ks.add(owner[b][0]); d.append(('w', owner[b][1]))
            else: d.append(('f', tpl[b]))
        if len(ks) > 1: return None
        d = tuple(d); k = next(iter(ks)) if ks else None
        if i == nch - 1: last = (d, k)
        elif k is None: fixed.append(d)
        else: per[k].append(d)
    sl = sorted(repr((slots[k][1], sorted(per[k]), last[0] if last[1] == k else None)) for k in range(len(slots)))
    return repr((L, sl, sorted(fixed), last[0] if last[1] is None else 'slot'))

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'words':
        W = load_words()
        for nc in (2, 3):
            ws = W[nc]; print('%d 文字の語: %d（漢字だけ %d、カタカナだけ %d）' % (nc, len(ws), sum(c == 'kanji' for w, c in ws), sum(c == 'kata' for w, c in ws)))
            pos = {w: i + 1 for i, (w, c) in enumerate(ws)}
            for t in ['磁石', '発光', '暗号', '解読', '懸賞', '情報', '科学', '国際', 'カオス', '研究所', '暗号化', '情報科学', '乱数', '写像', '秘密', '鍵']:
                if len(t) == nc: print('    %s: %s 位' % (t, pos.get(t)))
    elif mode in ('control',):
        n = int(sys.argv[2]); cnt = int(sys.argv[3]) if len(sys.argv) > 3 else 6; N2 = int(sys.argv[4]) if len(sys.argv) > 4 else 3000; N3 = int(sys.argv[5]) if len(sys.argv) > 5 else 1500
        seps = sys.argv[6:] or ['dot', 'none', 'sp', 'zsp']
        rng = random.Random(20261114 + n); W2 = word_bytes(2, N2); W3 = word_bytes(3, N3); ok = 0; done = 0; skipped = 0
        structs = [s for s in struct_list(n, seps)]; rng.shuffle(structs)
        for sizes, pos, sn, tpl, slots in structs:
            if done >= cnt: break
            wbs = [W2 if s == 4 else W3 for s in sizes]
            wbs, mok = last_slot_filter(tpl, slots, wbs)
            if not mok: continue
            plan = plan_groups(len(tpl), slots, [len(w) for w in wbs])
            if plan is None: skipped += 1; continue
            key = bytearray(tpl)
            for k, (o, ln) in enumerate(slots): key[o:o + ln] = rng.choice(wbs[k])
            key = bytes(key); r = key_schedule(key)
            if r['m'] != 1: continue                       # 最後の塊に枠が 2 つかかる形では、写像が 1 の鍵だけ対照にする
            tgt, _ = synth_target(key, rng); sj = sum(r['J']); jlo = max(0, sj - rng.randrange(241)); jwin = (jlo, jlo + 240)
            cands, stat, dt = run_kmitm(tgt, tpl, slots, wbs, plan, jwin, tag='ctl')
            good = [c for c in cands if sched_ok(tgt, c)]; hit = key in good; ok += hit; done += 1
            print('%s  長さ %d  %s  組 %s|%s 外 %s  nL %s nR %s 候補 %d 合格 %d  %.1f 秒  %s' % (show_struct(sizes, pos, sn), len(tpl), key.decode('cp932'), plan[1], plan[2], plan[3], stat.get('nL'), stat.get('nR'), len(cands), len(good), dt, '出た' if hit else '!!! 出ない')); sys.stdout.flush()
        print('対照 %d 本のうち %d 本で、埋めた鍵が出た（組に分けられず飛ばした形 %d）' % (done, ok, skipped))
    elif mode in ('run', 'count'):
        n = int(sys.argv[2]); N2 = int(sys.argv[3]); N3 = int(sys.argv[4]); seps = sys.argv[5:] or ['dot']
        W2 = word_bytes(2, N2); W3 = word_bytes(3, N3); log = os.path.join(HERE, 'logs', 'vdict_%d.jsonl' % n)
        done = {}                                           # 済みの形 → その記録
        if os.path.exists(log):
            for line in open(log, encoding='utf-8'):
                try: r = json.loads(line); done[(r['tpl'], r['N2'] if '語2' in r['form'] else 0, r['N3'] if '語3' in r['form'] else 0)] = r
                except Exception: pass
        todo = []; mapng = 0; nogroup = 0; total = 0
        for sizes, pos, sn, tpl, slots in struct_list(n, seps):
            total += 1
            wbs = [W2 if s == 4 else W3 for s in sizes]
            if any(not w for w in wbs): total -= 1; continue          # 語の一覧が空の枠がある形（N3 = 0 のときの 3 文字の語）は数えない
            wbs, mok = last_slot_filter(tpl, slots, wbs)
            if not mok: mapng += 1; continue
            plan = plan_groups(len(tpl), slots, [len(w) for w in wbs])
            if plan is None: nogroup += 1; continue
            todo.append((plan[0], sizes, pos, sn, tpl, slots, wbs, plan))
        tkey = lambda t: (t[4].hex(), N2 if 4 in t[1] else 0, N3 if 6 in t[1] else 0)
        todo.sort(key=lambda x: x[0]); new = [t for t in todo if tkey(t) not in done]
        # 塊を並べ替えると同じ鍵の集まりになる形（等価な形）は、1 つだけ回して記録を写す（11-14、4 回目）。VDICT_NODEDUP=1 で、全部を別々に回す
        fams = {}
        for t in todo:
            s = None if os.environ.get('VDICT_NODEDUP') else family_sig(t[4], t[5])
            fams.setdefault(s if s is not None else 'solo ' + t[4].hex(), []).append(t)
        work = []                                           # (見積もり, 回す形（None なら回さない）, 記録を写す形, 写す元の記録)
        for mem in fams.values():
            newm = [t for t in mem if tkey(t) not in done]
            if not newm: continue
            donem = [t for t in mem if tkey(t) in done]
            if donem: work.append((0.0, None, newm, done[tkey(donem[0])]))
            else: rep = min(newm, key=lambda t: t[0]); work.append((rep[0], rep, [t for t in newm if t is not rep], None))
        work.sort(key=lambda w: w[0])
        print('語 %d 個、2 文字の語 %d・3 文字の語 %d、区切り %s: 形 %d、写像で否定 %d、組に分けられない %d、回す %d（済み %d。等価な形をまとめると %d 組で、これから回すのは %d 組）、見積もり %.0f 秒' % (
            n, len(W2), len(W3), seps, total, mapng, nogroup, len(todo), len(todo) - len(new), len(fams), sum(w[1] is not None for w in work), sum(w[0] for w in work))); sys.stdout.flush()
        if mode == 'count':
            byL = {}
            for t in todo: byL.setdefault(len(t[4]), [0, 0.0]); byL[len(t[4])][0] += 1
            for w in work:
                if w[1] is not None: byL[len(w[1][4])][1] += w[0]
            print('  長さごと（形の数・これから回す分の見積もり秒）: ' + '  '.join('%d: %d・%.0f' % (L, v[0], v[1]) for L, v in sorted(byL.items()))); sys.exit(0)
        tlimit = float(os.environ.get('VDICT_TLIMIT', '540')); T0 = time.time(); nrun = 0; ncopy = 0; nhit = 0
        def run_one(f, t):
            cost, sizes, pos, sn, tpl, slots, wbs, plan = t
            cands, stat, dt = run_kmitm(TARGET_V, tpl, slots, wbs, plan, (90, 330), tag='run%d' % n)
            good = [c for c in cands if sched_ok(TARGET_V, c)]; hits = [c for c in good if keystream_ok_v(c)]
            rec = dict(n=n, form=show_struct(sizes, pos, sn), sep=sn, L=len(tpl), tpl=tpl.hex(), N2=N2, N3=N3, lists=[len(w) for w in wbs], groups=[plan[1], plan[2], plan[3]],
                       nL=int(stat.get('nL', 0)), nR=int(stat.get('nR', 0)), checks=int(stat.get('checks', 0)), cands=len(cands), sched=len(good), sec=round(dt, 1))
            if hits: rec['hit'] = [c.hex() for c in hits]
            f.write(json.dumps(rec, ensure_ascii=False) + '\n'); f.flush()
            print('  %-34s 長さ %d  語の数 %s  候補 %d 合格 %d  %.1f 秒（見積もり %.0f）%s' % (rec['form'], rec['L'], rec['lists'], len(cands), len(good), dt, cost, '  !!!! 当たり ' + ' '.join(c.decode('cp932', 'replace') for c in hits) if hits else '')); sys.stdout.flush()
            return rec, bool(hits)
        with open(log, 'a', encoding='utf-8') as f:
            for cost, rep, copies, src in work:
                if rep is not None:
                    if time.time() - T0 + cost * 0.6 > tlimit and nrun: break
                    src, h = run_one(f, rep); nrun += 1; nhit += h
                for t in copies:
                    cost2, sizes, pos, sn, tpl, slots, wbs, plan = t
                    if src.get('sched') or sorted(src['lists']) != sorted(len(w) for w in wbs):      # 鍵スケジュール全体に合う候補が出た組と、語の数がそろわない形は、写さずに回す
                        rec, h = run_one(f, t); nrun += 1; nhit += h; continue
                    rec = dict(n=n, form=show_struct(sizes, pos, sn), sep=sn, L=len(tpl), tpl=tpl.hex(), N2=N2, N3=N3, lists=[len(w) for w in wbs], groups=[plan[1], plan[2], plan[3]],
                               nL=0, nR=0, checks=src['checks'], cands=src['cands'], sched=src['sched'], sec=0.0, same_as=src['tpl'])
                    f.write(json.dumps(rec, ensure_ascii=False) + '\n'); f.flush(); ncopy += 1
                    print('  %-34s 長さ %d  語の数 %s  %s と等価（塊の並べ替え）。回さずに記録を写した' % (rec['form'], rec['L'], rec['lists'], src['form'])); sys.stdout.flush()
        print('今回 %d 形（回したもの %d、等価な形の記録を写したもの %d）、当たり %d、残り %d、%.0f 秒' % (nrun + ncopy, nrun, ncopy, nhit, len(new) - nrun - ncopy, time.time() - T0))
