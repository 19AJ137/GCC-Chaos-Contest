# 「語・語・語・固定の語」の形（語 3 個＋IISI・Gao・gao。順不同）を、辞書の全語で調べる（11-14、5 回目）。
#   語: wordfreq（日本語 large）の語のうち、全角だけでできたもの。長さは 1〜10 文字から選ぶ。
#       種類 kk = 漢字だけ・カタカナだけ（vdict.py と同じ）、all = ひらがなだけ・漢字かな交じりも入れる。
#   固定の語: IISI・Gao・gao など（半角）。位置は 0〜3（3 が末尾）。区切りは「・」のほか、半角の区切りも指定できる。
#   調べ方は vdict.py と同じ中間一致（kmitm）。語 3 個なので、外側 1 個＋左右 1 個ずつに分ける。
#   1 つの塊に語が 2 つかかる形（半角の区切りや、3 バイトの固定の語のあと）は、真ん中の語を「頭・中・尻」に分けて回す（split_mid）。
#     頭（尻）の数バイトは隣の語と同じ組に入れ、中は別の枠にする。分けた 3 つを自由に組み合わせるので、出た候補は「真ん中の語が辞書にあるか」で選び直す。
#   python3 vdict3.py count  種類 最短 最長 固定の語（, 区切り） 位置（, 区切り） 区切り（, 区切り）
#   python3 vdict3.py control 本数 （以下同じ）       対照（辞書の語で作った鍵を埋めて、出るかを見る）
#   python3 vdict3.py run    種類 最短 最長 固定の語 位置 区切り     本番（鍵 v の目標）。記録は logs/vdict3.jsonl
#     例: python3 vdict3.py run kk 2 8 IISI 0,1,2,3 dot      python3 vdict3.py run all 1 10 IISI 3 hdot,sp
import sys, os, json, time, random, itertools
import vdict
from vdict import TARGET_V, key_schedule, chunk_vals, map_of_i2, keystream_ok_v, synth_target, sched_ok, plan_groups, run_kmitm, family_sig, chunk_slots
vdict.PLAN_OUTER4 = True
HERE = vdict.HERE
SEPB = {'dot': '・'.encode('cp932'), 'zsp': '　'.encode('cp932'), 'none': b'', 'hdot': b'\xa5', 'sp': b' ', 'hy': b'-', 'us': b'_', 'pd': b'.', 'cm': b',', 'sl': b'/', 'pl': b'+',
        'am': b'&', 'at': b'@', 'cl': b':', 'sc': b';', 'eq': b'=', 'as': b'*', 'sh': b'#', 'ex': b'!', 'bar': b'|', 'tl': b'~'}
SEPT = {k: (v.decode('cp932') if k != 'hdot' else '･') for k, v in SEPB.items()}
LMIN, LMAX = 21, 40                                            # 20 バイト以下は、S と b0b1 の両方に合う塊の和がない（vklike.py sigma）。入力欄は 40 桁まで

def jclass(w):
    if all('一' <= c <= '鿿' or c == '々' for c in w): return 'kanji'
    if all('ァ' <= c <= 'ー' for c in w): return 'kata'
    if all('ぁ' <= c <= 'ゖ' or c == 'ー' for c in w): return 'hira'
    if all('一' <= c <= '鿿' or c == '々' or 'ぁ' <= c <= 'ゖ' or 'ァ' <= c <= 'ー' for c in w): return 'mixed'
    return None
_W = {}
def load(kind):
    """{文字数: [語のバイト列…頻度順]}。cp932 で 1 文字 2 バイトになるものだけ"""
    if kind not in _W:
        import wordfreq
        d = wordfreq.get_frequency_dict('ja', 'large'); out = {}
        ok = ('kanji', 'kata') if kind == 'kk' else ('kanji', 'kata', 'hira', 'mixed')
        for w, f in sorted(d.items(), key=lambda kv: (-kv[1], kv[0])):
            if jclass(w) not in ok: continue
            try: b = w.encode('cp932')
            except UnicodeEncodeError: continue
            if len(b) != 2 * len(w): continue
            out.setdefault(len(w), []).append(b)
        _W[kind] = out
    return _W[kind]

def make_struct(chars, pos, sep, fixed):
    """chars: 語ごとの文字数、pos: 固定の語の位置（0〜語の数）。戻り値: (ひな形, [(off, len)…語の枠])"""
    elems = [('w', 2 * c) for c in chars]; elems.insert(pos, ('k', fixed)); tpl = b''; slots = []
    for i, (t, v) in enumerate(elems):
        if i: tpl += sep
        if t == 'w': slots.append((len(tpl), v)); tpl += b'\x00' * v
        else: tpl += v
    return tpl, slots
def show(chars, pos, sepname, fixed):
    el = ['語%d' % c for c in chars]; el.insert(pos, fixed.decode()); return SEPT[sepname].join(el)

_LF = {}
def last_filter(tpl, slots, lists, want=1):
    """vdict.last_slot_filter と同じ（最後の塊に枠が 1 つだけかかるとき、写像が want になる語に絞る）。同じ条件の結果は使い回す"""
    L = len(tpl); cs = chunk_slots(L, slots); last = cs[-1]; base = 4 * (len(cs) - 1)
    if not last:
        c = int.from_bytes(tpl[base:], 'little'); return lists, map_of_i2(chunk_vals(c)[1][1]) == want
    if len(last) > 1: return lists, True
    k = last[0]; o, n = slots[k]
    fixedpart = bytes(tpl[base:]); key = (fixedpart, o - base, n, id(lists[k]))
    if key not in _LF:
        keep = []; lo = max(o, base); hi = min(o + n, L)
        for w in lists[k]:
            ch = bytearray(fixedpart); ch[lo - base:hi - base] = w[lo - o:hi - o]
            if map_of_i2(chunk_vals(int.from_bytes(bytes(ch), 'little'))[1][1]) == want: keep.append(w)
        _LF[key] = keep
    out = list(lists); out[k] = _LF[key]; return out, True

_SP = {}
def split_mid(tpl, slots, lists):
    """真ん中の語（枠 1）が隣の語と塊を共有するとき、真ん中の語を「頭の数バイト・中・尻の数バイト」に分ける。
       戻り値: None（共有がない）、'bad'（扱えない形）、または (新しい枠, 新しい一覧, (真ん中の語の位置, 長さ, 語の集合))"""
    L = len(tpl); cs = chunk_slots(L, slots)
    sh01 = [i for i, s in enumerate(cs) if 0 in s and 1 in s]; sh12 = [i for i, s in enumerate(cs) if 1 in s and 2 in s]; sh02 = [i for i, s in enumerate(cs) if 0 in s and 2 in s]
    if not sh01 and not sh12 and not sh02: return None
    if sh02 or len(sh01) > 1 or len(sh12) > 1: return 'bad'
    o, n = slots[1]; k1 = 0; k2 = 0
    if sh01: k1 = min(o + n, 4 * sh01[0] + 4) - o               # 真ん中の語が、前の語と共有する塊に持つバイト数（語の先頭側）
    if sh12: k2 = o + n - max(o, 4 * sh12[0])                   # 後ろの語と共有する塊に持つバイト数（語の末尾側）
    if k1 + k2 > n: return 'bad'
    key = (id(lists[1]), k1, k2)
    if key not in _SP:
        ws = lists[1]
        _SP[key] = (sorted(set(w[:k1] for w in ws)) if k1 else None, sorted(set(w[k1:n - k2] for w in ws)) if n - k1 - k2 > 0 else None, sorted(set(w[n - k2:] for w in ws)) if k2 else None, frozenset(ws))
    Y, M, U, wset = _SP[key]
    s = [slots[0]]; l = [lists[0]]
    if k1: s.append((o, k1)); l.append(Y)
    if M is not None: s.append((o + k1, n - k1 - k2)); l.append(M)
    if k2: s.append((o + n - k2, k2)); l.append(U)
    s.append(slots[2]); l.append(lists[2])
    return s, l, (o, n, wset)

def shapes(kind, cmin, cmax, fixeds, poss, seps):
    W = load(kind)
    for fixed in fixeds:
        for sepname in seps:
            sep = SEPB[sepname]
            for chars in itertools.product(range(cmin, cmax + 1), repeat=3):
                if any(c not in W for c in chars): continue
                L = sum(2 * c for c in chars) + len(fixed) + 3 * len(sep)
                if L < LMIN or L > LMAX: continue
                for pos in poss:
                    tpl, slots = make_struct(chars, pos, sep, fixed)
                    yield fixed, sepname, chars, pos, tpl, slots

def classify(kind, cmin, cmax, fixeds, poss, seps):
    """戻り値: (回せる形の一覧, 数の内訳)"""
    W = load(kind); todo = []; st = dict(total=0, mapng=0, nogroup=0, split=0, ex=[])
    for fixed, sepname, chars, pos, tpl, slots in shapes(kind, cmin, cmax, fixeds, poss, seps):
        st['total'] += 1
        lists, mok = last_filter(tpl, slots, [W[c] for c in chars])
        if not mok or any(not l for l in lists): st['mapng'] += 1; continue
        form = show(chars, pos, sepname, fixed); sp = split_mid(tpl, slots, lists)
        if sp == 'bad': st['nogroup'] += 1; st['ex'].append(form); continue
        rs, rl, post = (slots, lists, None) if sp is None else sp
        plan = plan_groups(len(tpl), rs, [len(l) for l in rl])
        if plan is None: st['nogroup'] += 1; st['ex'].append(form); continue
        sig = family_sig(tpl, rs)
        st['split'] += sp is not None
        todo.append(dict(cost=plan[0], fixed=fixed, sep=sepname, chars=chars, pos=pos, form=form, tpl=tpl, slots0=slots, lists0=lists, slots=rs, lists=rl, post=post, plan=plan,
                         sig=sig if sig is not None else 'solo ' + tpl.hex()))
    return todo, st

def why(tgt, key):
    r = key_schedule(key); out = []
    out.append('写像 %d%s' % (r['m'], '' if r['m'] == tgt['m'] else '（×）'))
    out.append('S %08x%s' % (r['S'], '' if r['S'] == tgt['S'] else '（×）'))
    out.append('b0b1 %04x%s' % (r['b0b1'], '' if tgt['Tlo'] <= r['b0b1'] <= tgt['Thi'] else '（×）'))
    d1 = min((abs(r['k'][0] - c) - h for c, h in tgt['cells1'])); out.append('k780 %s' % ('セルの中' if d1 <= 0 else 'セルの外 %.1e（×）' % d1))
    d3 = abs(r['k'][2] - tgt['c3']) - tgt['h3']; out.append('k7b0 %s' % ('セルの中' if d3 <= 0 else 'セルの外 %.1e（×）' % d3))
    return '、'.join(out)

def solve(t, tgt, jwin, tag):
    """1 形を回す。戻り値: (候補（真ん中の語を選び直す前の数）, 候補（選び直したあと）, 記録, 秒)"""
    cands, stat, dt = run_kmitm(tgt, t['tpl'], t['slots'], t['lists'], t['plan'], jwin, tag=tag)
    raw = len(cands)
    if t['post']: o, n, wset = t['post']; cands = [c for c in cands if bytes(c[o:o + n]) in wset]
    return raw, cands, stat, dt

if __name__ == '__main__':
    mode = sys.argv[1]; a = sys.argv[2:]
    cnt = int(a.pop(0)) if mode == 'control' else 0
    kind = a[0]; cmin = int(a[1]); cmax = int(a[2]); fixeds = [x.encode() for x in a[3].split(',')]; poss = [int(x) for x in a[4].split(',')]; seps = a[5].split(',')
    tag = '%s %d-%d %s pos %s %s' % (kind, cmin, cmax, a[3], a[4], a[5])
    todo, st = classify(kind, cmin, cmax, fixeds, poss, seps)
    fams = {}
    for t in todo: fams.setdefault(t['sig'], []).append(t)
    W = load(kind)
    print('%s: 語の数 %s' % (tag, ' '.join('%d文字:%d' % (c, len(W[c])) for c in range(cmin, cmax + 1) if c in W)))
    head = '  形 %d、写像で否定 %d、組に分けられない %d、回せる %d（うち、真ん中の語を分けて回す形 %d。等価な形をまとめると %d 組）' % (st['total'], st['mapng'], st['nogroup'], len(todo), st['split'], len(fams))
    if mode == 'count':
        print(head + '、見積もり %.0f 秒' % sum(min(t['cost'] for t in m) for m in fams.values()))
        byL = {}
        for t in todo: byL[len(t['tpl'])] = byL.get(len(t['tpl']), 0) + 1
        print('  長さごとの回せる形: ' + '  '.join('%d: %d' % kv for kv in sorted(byL.items())))
        if st['nogroup']: print('  組に分けられない形の例: ' + '、'.join(st['ex'][:6]))
    elif mode == 'control':
        rng = random.Random(20261114 + cnt + len(todo)); ok = 0; done = 0; cand = list(todo); rng.shuffle(cand); want_split = cnt // 2 if st['split'] else 0; nsp = 0
        for t in cand:
            if done >= cnt: break
            if t['post'] is None and want_split - nsp >= cnt - done: continue          # 分けて回す形を、対照の半分に入れる
            key = bytearray(t['tpl'])
            for k, (o, ln) in enumerate(t['slots0']): key[o:o + ln] = rng.choice(t['lists0'][k])
            key = bytes(key); r = key_schedule(key)
            if r['m'] != 1: continue
            tgt, _ = synth_target(key, rng); sj = sum(r['J']); jlo = max(0, sj - rng.randrange(241)); jwin = (jlo, jlo + 240)
            raw, cands, stat, dt = solve(t, tgt, jwin, 'ctl3')
            good = [c for c in cands if sched_ok(tgt, c)]; hit = key in good; ok += hit; done += 1; nsp += t['post'] is not None
            print('  %-28s 長さ %d  %s  %s 組 %s|%s 外 %s  語の数 %s  候補 %d→%d 合格 %d  %.1f 秒  %s' % (t['form'], len(key), key.decode('cp932', 'replace'), '分けた枠' if t['post'] else '', t['plan'][1], t['plan'][2], t['plan'][3],
                  [len(l) for l in t['lists']], raw, len(cands), len(good), dt, '出た' if hit else '!!! 出ない')); sys.stdout.flush()
        print('対照 %d 本（うち、真ん中の語を分けて回す形 %d 本）のうち %d 本で、埋めた鍵が出た' % (done, nsp, ok))
    elif mode == 'run':
        log = os.path.join(HERE, 'logs', 'vdict3.jsonl'); done = {}
        if os.path.exists(log):
            for line in open(log, encoding='utf-8'):
                try: r = json.loads(line); done[(r['tpl'], r['kind'], tuple(r['lists']))] = r
                except Exception: pass
        tkey = lambda t: (t['tpl'].hex(), kind, tuple(len(l) for l in t['lists']))
        work = []
        for mem in fams.values():
            newm = [t for t in mem if tkey(t) not in done]
            if not newm: continue
            donem = [t for t in mem if tkey(t) in done]
            if donem: work.append((0.0, None, newm, done[tkey(donem[0])]))
            else: rep = min(newm, key=lambda t: t['cost']); work.append((rep['cost'], rep, [t for t in newm if t is not rep], None))
        work.sort(key=lambda w: w[0])
        print(head + '。これから回すのは %d 組、見積もり %.0f 秒' % (sum(w[1] is not None for w in work), sum(w[0] for w in work))); sys.stdout.flush()
        tlimit = float(os.environ.get('VDICT_TLIMIT', '1e9')); T0 = time.time(); nrun = 0; ncopy = 0; nhit = 0; ncand = 0; nchk = 0; nleft = 0
        def rec_of(t, **kw):
            return dict(n=3, kind=kind, fixed=t['fixed'].decode(), sep=t['sep'], pos=t['pos'], chars=list(t['chars']), form=t['form'], L=len(t['tpl']), tpl=t['tpl'].hex(), split=t['post'] is not None,
                        lists=[len(l) for l in t['lists']], groups=[t['plan'][1], t['plan'][2], t['plan'][3]], **kw)
        def run_one(f, t):
            raw, cands, stat, dt = solve(t, TARGET_V, (90, 330), 'run3')
            good = [c for c in cands if sched_ok(TARGET_V, c)]; hits = [c for c in good if keystream_ok_v(c)]
            rec = rec_of(t, checks=int(stat.get('checks', 0)), cands_raw=raw, cands=len(cands), sched=len(good), sec=round(dt, 1))
            if cands: rec['cand_keys'] = [c.hex() for c in sorted(cands)[:40]]
            if hits: rec['hit'] = [c.hex() for c in hits]
            f.write(json.dumps(rec, ensure_ascii=False) + '\n'); f.flush()
            if cands or hits:
                print('  %-28s 長さ %d  語の数 %s  和が合った組 %d  候補 %d 合格 %d  %.1f 秒%s' % (rec['form'], rec['L'], rec['lists'], rec['checks'], len(cands), len(good), dt, '  !!!! 当たり ' + ' '.join(c.decode('cp932', 'replace') for c in hits) if hits else ''))
                for c in sorted(cands)[:40]: print('      候補 %s  [%s]  %s' % (c.hex(), c.decode('cp932', 'replace'), why(TARGET_V, c)))
                sys.stdout.flush()
            return rec, bool(hits)
        with open(log, 'a', encoding='utf-8') as f:
            for cost, rep, copies, src in work:
                if rep is not None:
                    if time.time() - T0 > tlimit and nrun: nleft += 1 + len(copies); continue
                    src, h = run_one(f, rep); nrun += 1; nhit += h; ncand += src['cands']; nchk += src['checks']
                    if nrun % 50 == 0: print('  … %d 組を回した（%.0f 秒）。和が合った組 %d、候補 %d' % (nrun, time.time() - T0, nchk, ncand)); sys.stdout.flush()
                for t in copies:
                    if src.get('sched') or sorted(src['lists']) != sorted(len(l) for l in t['lists']):
                        rec, h = run_one(f, t); nrun += 1; nhit += h; ncand += rec['cands']; nchk += rec['checks']; continue
                    rec = rec_of(t, checks=src['checks'], cands_raw=src.get('cands_raw', src['cands']), cands=src['cands'], sched=src['sched'], sec=0.0, same_as=src['tpl'])
                    f.write(json.dumps(rec, ensure_ascii=False) + '\n'); f.flush(); ncopy += 1
        print('今回 %d 形（回したもの %d、等価な形の記録を写したもの %d）: 和が合った組 %d（回した分）、候補 %d、当たり %d、残り %d 形、%.0f 秒' % (nrun + ncopy, nrun, ncopy, nchk, ncand, nhit, nleft, time.time() - T0))
