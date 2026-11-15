# 鍵 v の「一部が分かっている」形の仮説を、族ごとに作って vres で調べる（2026-11-14）。
#   python3 vcamp.py count [族 …]        族ごとの仮説の数と見積もり
#   python3 vcamp.py run 族 [制限秒]       その族を回す（記録は logs/vcamp_族.jsonl。途中で止めても続きから回せる）
# 仮説 = 既知のバイト列と未知の区間（バイト数）の並び。未知の区間は全バイト値で解く（vpart.build_hyp が解き方を決める）。
import sys, os, itertools, time, json
from gks import key_schedule, chunk_vals, map_of_i2
from vpart import build_hyp, run_vres, TARGET_V, keystream_ok_v, show, segs_to_known, textlike, iter_logs
E = lambda s: s.encode('cp932')
K = E('磁石・発光・IISI'); DOT = E('・')
KW = [E('磁石'), E('発光'), E('IISI')]
# 鍵 k と一緒に入りそうな語（社名・製品名・人名・平文の語・作者のほかの鍵）
EXTRA_A = ['GCC', 'GCC6', 'GCC7', 'IISI', 'Gao', 'GAO', 'gao', 'Chaos', 'CHAOS', 'chaos', 'Japan', 'JAPAN', 'JAPEN', '1996', '1997', 'Zhenyu', 'Yang', 'Changxing',
           'GaoZhenyu', 'ZhenyuGao', 'ChangxingYang', 'Dr.Gao', 'Ph.D.', 'Gao@IISI', 'gao1234', '1234', '123']
EXTRA_J = ['カオス', '暗号', 'カオス暗号', '懸賞', '懸賞金', '解読', '賞金', '挑戦', '高振宇', '国際情報科学研究所', 'ＧＣＣ', 'ＩＩＳＩ', 'インフォガード', '日本', '中国', '読み込み中']
EXTRA = [E(x) for x in EXTRA_A + EXTRA_J]
EXTRA_S = [E(x) for x in ['GCC', 'GCC6', 'GCC7', 'IISI', 'Gao', 'Chaos', '1996', '1997', 'カオス', '暗号', 'カオス暗号', '懸賞', '解読', '高振宇', 'ＧＣＣ', '日本']]   # 2 語足すとき用（絞ったもの）
PHRASES = ['国際情報科学研究所', '株式会社国際情報科学研究所', '(株)国際情報科学研究所', '（株）国際情報科学研究所', '㈱国際情報科学研究所',
           'Information Science Institute', "Gao's Chaos Cryptosystem", "Gao's Chaos cryptosystem", "GCC(Gao's Chaos Cryptosystem)", "GCC(Gao's Chaos cryptosystem)",
           "GCC6(Gao's Chaos Cryptosystem)", 'Gao Zhenyu, Ph.D.', 'Gao Zhenyu, Ph.D. JAPEN', 'Gao Zhenyu, Ph.D. JAPAN', 'Gao Zhenyu', 'Zhenyu Gao', 'I am Gao@IISI.',
           'GCCカオス暗号', 'ＧＣＣカオス暗号', 'カオスインフォガード', 'カオス・インフォガード', 'デジタル型カオス暗号', '共通キー型のデジタル型カオス暗号',
           'You can not decrypte GCC7!', "No, you cann't!", 'Because this is GCC6 only!', 'Do you think you can decrypte GCC?', 'Thank you for your help',
           'digital chaos cryptosystem', 'chaos cryptosystem', 'GCC7ならば絶対解けることができない！', 'これはGCC6過ぎではない。', '読み込み中・・・',
           'GoodLuckToIISIGoodLuckToChangxingYang', 'GoodLuckToChangxingYang', 'GoodLuckToIISI', 'International Information Science']
ENDS = ['', 'IISI', '・IISI', ' IISI', 'GCC', 'GCC6', 'GCC7', '1996', '1997']
GL_NAMES = ['IISI', 'GCC', 'GCC6', 'GCC7', 'Gao', 'GAO', 'You', 'Me', 'All', 'Everyone', 'Japan', 'China', 'Yang', 'Zhenyu', 'Changxing', 'ChangxingYang', 'YangChangxing',
            'GaoZhenyu', 'ZhenyuGao', 'Dr.Gao', 'Chaos', 'InfoGuard']

def join(slots, style):
    """slots: bytes（既知の語）か int（未知の区間のバイト数）。style: dot＝全部を「・」でつなぐ、dotabs＝既知どうしだけ「・」（未知の区間は区切りも含む）、
    none＝そのままつなぐ、space＝全部を空白でつなぐ"""
    sep = {'dot': DOT, 'dotabs': DOT, 'none': b'', 'space': b' '}[style]; out = []
    for i, s in enumerate(slots):
        if i:
            if style == 'dotabs' and (isinstance(slots[i - 1], int) or isinstance(s, int)): pass
            elif sep: out.append(sep)
        out.append(s)
    # 既知どうしをまとめる
    segs = []
    for s in out:
        if isinstance(s, int):
            if segs and isinstance(segs[-1], int): segs[-1] += s
            else: segs.append(s)
        else:
            if segs and not isinstance(segs[-1], int): segs[-1] = segs[-1] + s
            else: segs.append(s)
    return segs
def seglen(segs): return sum(s if isinstance(s, int) else len(s) for s in segs)
def unk(segs): return sum(s for s in segs if isinstance(s, int))

def fam_A0():
    """鍵 k をそのまま含み、ほかは全部未知: X + k + Y"""
    for u1 in range(0, 15):
        for u2 in range(0, 15 - u1):
            if u1 + u2 == 0: continue
            yield [x for x in (u1, K, u2) if x]
def fam_A1():
    """鍵 k の途中に未知の区間が 1 つ入る: k[:c] + X + k[c:]"""
    for c in range(1, 16):
        for u in range(1, 15): yield [K[:c], u, K[c:]]
def fam_A2():
    """鍵 k と既知の語 1 個と未知の区間 1 つを並べる（順番と区切りは全通り）"""
    for e in EXTRA:
        for perm in itertools.permutations([K, e, 0]):
            for style in ('dot', 'dotabs', 'none', 'space'):
                for u in range(1, 14):
                    segs = join([u if x == 0 else x for x in perm], style)
                    if 25 <= seglen(segs) <= 40: yield segs
def fam_A3():
    """鍵 k と既知の語 2 個と未知の区間 1 つを並べる"""
    for e1, e2 in itertools.permutations(EXTRA_S, 2):
        for perm in itertools.permutations([K, e1, e2, 0]):
            if perm.index(e1) > perm.index(e2): continue        # e1, e2 の順は外側の permutations で回している
            for style in ('dot', 'dotabs', 'none'):
                for u in range(1, 14):
                    segs = join([u if x == 0 else x for x in perm], style)
                    if 25 <= seglen(segs) <= 40: yield segs
def fam_B0():
    """鍵 k の 3 語（磁石・発光・IISI）と未知の区間 1 つを並べる（順番と区切りは全通り）"""
    for perm in itertools.permutations(KW + [0]):
        for style in ('dot', 'dotabs', 'none'):
            for u in range(1, 15):
                segs = join([u if x == 0 else x for x in perm], style)
                if 17 <= seglen(segs) <= 40: yield segs
def fam_B1():
    """鍵 k の 3 語と、未知の区間 2 つ（2・4・6 バイト＝全角 1〜3 文字）を並べる"""
    for u1, u2 in itertools.product((2, 4, 6), repeat=2):
        for perm in itertools.permutations(KW + [-1, -2]):
            if perm.index(-1) > perm.index(-2): continue
            for style in ('dot', 'none'):
                segs = join([u1 if x == -1 else u2 if x == -2 else x for x in perm], style)
                if 17 <= seglen(segs) <= 40: yield segs
                # 半角 1 文字（未知）を IISI の前か末尾に足した形（長さを 4 で割って 1 余らせるための 1 文字）
                for pos in ('before_iisi', 'end', 'start'):
                    sl = [u1 if x == -1 else u2 if x == -2 else x for x in perm]; s2 = join(sl, style)
                    if pos == 'end': s2 = s2 + [1]
                    elif pos == 'start': s2 = [1] + s2
                    else:
                        out = []
                        for s in s2:
                            if not isinstance(s, int) and KW[2] in s:
                                i = s.index(KW[2]); out += [x for x in (s[:i], 1, s[i:]) if x]
                            else: out.append(s)
                        s2 = out
                    s2 = join(s2, 'none')
                    if 17 <= seglen(s2) <= 40: yield s2
def fam_B2():
    """鍵 k の 3 語と既知の語 1 個と未知の区間 1 つを並べる"""
    for e in EXTRA:
        if e == KW[2]: continue
        for perm in itertools.permutations(KW + [e, 0]):
            for style in ('dot', 'dotabs', 'none'):
                for u in range(1, 14):
                    segs = join([u if x == 0 else x for x in perm], style)
                    if 25 <= seglen(segs) <= 40: yield segs
def fam_C():
    """4.02a の固定鍵（GoodLuckToIISIGoodLuckToChangxingYang）の言い回し"""
    for G0 in ('GoodLuckTo', 'goodluckto', 'GOODLUCKTO', 'Good Luck To ', 'Good luck to ', 'GoodLuck', 'Good Luck '):
        G = E(G0); sp = b' ' if G0.endswith(' ') else b''
        for a0 in GL_NAMES:
            a = E(a0.lower() if G0 == 'goodluckto' else a0.upper() if G0 == 'GOODLUCKTO' else a0)
            for u in range(1, 15):
                for segs in ([G + a + sp + G, u], [G, u, sp + G + a], [G + a + sp, u], [u, sp + G + a], [G + a, u, sp + G + a], [G + a + sp + G + a + sp, u]):
                    if 25 <= seglen(segs) <= 40: yield segs
        for u1 in range(1, 10):
            for u2 in range(1, 15 - u1):
                segs = [G, u1, sp + G, u2]
                if 25 <= seglen(segs) <= 40: yield segs
        for u in range(1, 15):
            for segs in ([G, u], [u, G], [G + G, u], [G, u, G]):
                if 25 <= seglen(segs) <= 40: yield segs
    F = E('GoodLuckToIISIGoodLuckToChangxingYang')
    for u in range(1, 4):
        yield [F, u]; yield [u, F]
    for c in (10, 14, 24):
        for u in range(1, 4): yield [F[:c], u, F[c:]]
def fam_D():
    """平文・社名・製品名の語句（12 バイト以上）＋未知の区間＋（あれば）末尾の語"""
    for p0 in PHRASES:
        p = E(p0)
        for q0 in ENDS:
            q = E(q0)
            for u in range(1, 15):
                cands = [[p, u, q], [q, u, p], [u, p, q], [p, q, u]] if q else [[p, u], [u, p]]
                for segs in cands:
                    if 25 <= seglen(segs) <= 40: yield segs
        for u1 in range(1, 14):          # 語句の前後に未知の区間
            for u2 in range(1, 15 - u1):
                segs = [u1, p, u2]
                if 25 <= seglen(segs) <= 40: yield segs
def fam_S():
    """重い個別の仮説（未知 14 バイト: 塊 3 個＋端数 2 バイト）"""
    yield [K, 14]
    yield [E('GoodLuckToIISIGoodLuckTo'), 14]
def fam_E():
    """既知の語 2〜3 個（鍵 k の語を含まないものも）を「・」か、区切りなしでつないだもの＋未知の区間"""
    toks = EXTRA + KW[:2]
    for m in (2, 3):
        for seq in itertools.permutations(toks if m == 2 else [E(x) for x in ['GCC', 'GCC6', 'GCC7', 'IISI', 'Gao', 'Chaos', '1996', 'カオス', '暗号', 'カオス暗号', '懸賞', '解読', '高振宇', '磁石', '発光']], m):
            for style in ('dot', 'none'):
                base = join(list(seq), style)[0]
                if len(base) < 12: continue
                for u in range(1, 14):
                    for segs in ([base, u], [u, base], [base + DOT, u] if style == 'dot' else None, [u, DOT + base] if style == 'dot' else None):
                        if segs and 29 <= seglen(segs) <= 40: yield segs
FAMS = dict(A0=fam_A0, A1=fam_A1, A2=fam_A2, A3=fam_A3, B0=fam_B0, B1=fam_B1, B2=fam_B2, C=fam_C, D=fam_D, E=fam_E, S=fam_S)

def last_chunk_map(known):
    """最後の塊が全部既知なら、その写像の番号。未知が混じれば None"""
    L = len(known); tail = known[4 * ((L - 1) // 4):]
    if any(b is None for b in tail): return None
    c = sum(b << (8 * j) for j, b in enumerate(tail)); (f, i, s, t) = chunk_vals(c)
    return map_of_i2(i[1])
def gen_segs(fam):
    """族の仮説（重複を除く）。戻り値: (仮説のリスト（segs）, 集計)。最後の塊が既知で写像が合わないもの・未知が無いものは除く"""
    seen = set(); out = []; st = dict(made=0, uniq=0, map_ng=0, nozero=0)
    for segs in FAMS[fam]():
        st['made'] += 1
        known = tuple(segs_to_known(segs))
        if known in seen: continue
        seen.add(known); st['uniq'] += 1
        if all(b is not None for b in known): st['nozero'] += 1; continue
        m = last_chunk_map(known)
        if m is not None and m != TARGET_V['m']: st['map_ng'] += 1; continue
        out.append(segs)
    return out, st
# 族ごとの方針。inner＝鍵 k だけの族（重いものまで回す）、core＝語を 1 個足した族と言い回し・語句の族、big＝語を 2 個足した族など（数が多いので軽いものだけ）
CLASS = dict(A0='inner', A1='inner', B0='inner', A2='core', B1='core', C='core', D='core', A3='big', B2='big', E='big', S='special')
def choose(cls, h16, h17):
    """回す格子を決める。戻り値: 16 / 17 / None（見送り）"""
    if cls == 'inner': return 17 if h17['cost'] <= 130 else 16 if h16['cost'] <= 130 else None
    if cls == 'core': return 17 if h17['cost'] <= 1.0 else 16 if h16['cost'] <= 6 else None
    if cls == 'big':
        if h17['nU'] <= 2 and h17['cost'] <= 1.0: return 17
        if h16['nU'] == 3 and h16['E'] == 0 and h16['cost'] <= 0.05: return 16
        if h17['nU'] == 3 and h17['E'] == 0 and h17['cost'] <= 0.3: return 17
        return None
    if cls == 'special':                  # 個別の重い形。環境変数 VCAMP_SPECIAL_LAT=17 で格子 17 でも回せる
        lat = int(os.environ.get('VCAMP_SPECIAL_LAT', '16')); return lat if (h17 if lat == 17 else h16)['cost'] <= 6000 else None
def plan(fam, segs_list):
    """仮説ごとに格子と解き方を決める。戻り値: (回す仮説 [(lat, h)], 見送り [h])"""
    run = []; skip = []
    for i, segs in enumerate(segs_list):
        h17 = build_hyp('%s_%06d' % (fam, i), segs, lat=17); h16 = build_hyp('%s_%06d' % (fam, i), segs, lat=16)
        lat = choose(CLASS[fam], h16, h17)
        if lat == 17: run.append((17, h17))
        elif lat == 16: run.append((16, h16))
        else: skip.append(h16 if h16['cost'] <= h17['cost'] else h17)
    return run, skip
def pat(h): return ''.join('..' if b is None else '%02x' % b for b in h['known'])
def load_done(logdir='logs'):
    done = {}
    for r in iter_logs(logdir): done[r['pat']] = max(done.get(r['pat'], 0), r['lat'])
    return done

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'count':
        tot = 0.0
        for fam in (sys.argv[2:] or list(FAMS)):
            segs, st = gen_segs(fam); run, skip = plan(fam, segs)
            by = {}
            for lat, h in run: by[(lat, h['nU'], h['E'])] = by.get((lat, h['nU'], h['E']), 0) + 1
            Ls = sorted(set(h['L'] for lat, h in run)); est = sum(h['cost'] for lat, h in run); tot += est
            print('%-3s 作った %7d  重複を除く %7d  写像で否定 %6d  回す %7d（見積もり %6.0f 秒）  見送り %6d  長さ %s〜%s' % (
                fam, st['made'], st['uniq'], st['map_ng'], len(run), est, len(skip), Ls[0] if Ls else '-', Ls[-1] if Ls else '-'))
            print('      格子,U,E:数  ' + ' '.join('%d,%d,%d:%d' % (k[0], k[1], k[2], v) for k, v in sorted(by.items())))
            sys.stdout.flush()
        print('見積もりの合計 %.0f 秒' % tot)
    elif mode == 'run':
        fam = sys.argv[2]; tlimit = float(sys.argv[3]) if len(sys.argv) > 3 else 540.0
        log = 'logs/vcamp_%s.jsonl' % fam
        segs, st = gen_segs(fam); run, skip = plan(fam, segs)
        done = load_done()
        todo = [(lat, h) for lat, h in run if done.get(pat(h), 0) < lat]
        todo.sort(key=lambda x: (x[0], x[1]['cost']))
        print('%s: 回す %d、済み %d、今回 %d（見積もり %.0f 秒）、見送り %d、写像で否定 %d' % (fam, len(run), len(run) - len(todo), len(todo), sum(h['cost'] for lat, h in todo), len(skip), st['map_ng'])); sys.stdout.flush()
        T0 = time.time(); pos = 0; nhit = 0; ncand = 0
        with open(log, 'a') as f:
            while pos < len(todo) and time.time() - T0 < tlimit:
                lat = todo[pos][0]; batch = []; est = 0.0
                while pos < len(todo) and todo[pos][0] == lat and (not batch or (est + todo[pos][1]['cost'] <= 40.0 and len(batch) < 3000)): batch.append(todo[pos][1]); est += todo[pos][1]['cost']; pos += 1
                keys, stat, dt = run_vres(TARGET_V, batch, L=lat, tag='camp_%s' % fam)
                for h in batch:
                    ks = sorted(keys.get(h['id'], set())); good = [k.hex() for k in ks if keystream_ok_v(k)]
                    sst = stat.get(h['id'], {})
                    rec = dict(fam=fam, lat=lat, pat=pat(h), known=show(h), L=h['L'], nU=h['nU'], E=h['E'], n=sst.get('assignments'), raw=sst.get('rawkeys', 0))
                    if ks: rec['cand'] = [k.hex() for k in ks]; ncand += 1
                    if good:
                        rec['hit'] = good; rec['text'] = max(textlike(bytes.fromhex(g)) for g in good)
                    f.write(json.dumps(rec, ensure_ascii=False) + '\n')
                    if good:
                        nhit += 1; tl = [g for g in good if textlike(bytes.fromhex(g))]
                        if tl: print('!!!! 鍵ストリームが合い、文字列として読める鍵', show(h), [(g, bytes.fromhex(g).decode('cp932', 'replace')) for g in tl])
                        else: print('  鍵ストリームが合う鍵（文字列としては読めない＝等価な鍵） %d 本  %s  例 %s' % (len(good), show(h), good[0]))
                f.flush()
                print('  %s 格子 %d  %d / %d  %.0f 秒（この回 %.1f 秒、見積もり %.1f 秒）' % (fam, lat, pos, len(todo), time.time() - T0, dt, est)); sys.stdout.flush()
        print('%s 終わり: 今回 %d 件、当たり %d 件、候補だけ %d 件、残り %d 件、%.0f 秒' % (fam, pos, nhit, ncand, len(todo) - pos, time.time() - T0))
