#!/usr/bin/env python3
# 鍵 v が「全角の語を『・』でつないだ 34 バイト（IISI なし）」だという仮説を、辞書の語で調べる（docs/password.md §8.10。11-15）。
#   形: 語 k 個（k = 2〜5）。文字数の合計は 18 − k（2 バイト × 文字数 ＋ 区切り 2 バイト ×（k − 1）= 34）。語 6 個は k6.py。
#   語の一覧: vdict3.load（wordfreq の日本語 large。kk = 漢字だけ・カタカナだけ、all = ひらがな・漢字かな交じりも）。長さごとに上位 N 語で切れる。
#   調べ方は vdict.py と同じ中間一致（kmitm）。出てきた候補は、鍵スケジュール全体と鍵ストリーム 335 語で判定する。
#   python3 v34.py count   k 種類 最短 最長 [長さ:語数,…]     形の数と見積もり       例: python3 v34.py count 4 kk 2 8 2:11000,3:11000,4:11000
#   python3 v34.py run     k 種類 最短 最長 [長さ:語数,…]     本番（鍵 v の目標）。記録は logs/v34.jsonl。途中で止めても続きから回せる
#   python3 v34.py control k 種類 最短 最長 [長さ:語数,…] [本数]  対照（辞書の語で作った鍵を埋めて、出るかを見る）
#   環境変数 V34_TLIMIT（秒。既定 520）を超えそうになったら、次の形に進まずに終わる。
import sys, os, json, time, random, itertools
import vdict, vdict3
from vdict import TARGET_V, keystream_ok_v, synth_target, sched_ok, plan_groups, run_kmitm, family_sig, last_slot_filter
from gks import key_schedule
HERE = vdict.HERE; SEP = '・'.encode('cp932'); TOTAL = 34; LOG = os.path.join(HERE, 'logs', 'v34.jsonl')

def compositions(total, k, lo, hi):
    if k == 1:
        if lo <= total <= hi: yield (total,)
        return
    for a in range(lo, min(hi, total - lo * (k - 1)) + 1):
        for rest in compositions(total - a, k - 1, lo, hi): yield (a,) + rest

def make(chars):
    tpl = b''; slots = []
    for i, c in enumerate(chars):
        if i: tpl += SEP
        slots.append((len(tpl), 2 * c)); tpl += b'\x00' * (2 * c)
    assert len(tpl) == TOTAL
    return tpl, slots

def shapes(k, kind, lo, hi, caps):
    W = vdict3.load(kind)
    for chars in compositions(18 - k, k, lo, hi):
        lists = [W.get(c, [])[:caps.get(c, 10**9)] for c in chars]
        if any(not l for l in lists): continue
        tpl, slots = make(chars)
        yield chars, tpl, slots, lists

def fill(tpl, slots, words):
    key = bytearray(tpl)
    for (o, n), w in zip(slots, words): key[o:o + n] = w
    return bytes(key)

def solve(tgt, t, tag, jwin=(90, 330)):
    chars, tpl, slots, lists, plan = t
    cands, stat, dt = run_kmitm(tgt, tpl, slots, lists, plan, jwin, tag=tag)
    good = [c for c in cands if sched_ok(tgt, c)]
    return cands, good, stat, dt

if __name__ == '__main__':
    mode = sys.argv[1]; k = int(sys.argv[2]); kind = sys.argv[3]; lo = int(sys.argv[4]); hi = int(sys.argv[5])
    caps = {int(a.split(':')[0]): int(a.split(':')[1]) for a in sys.argv[6].split(',')} if len(sys.argv) > 6 and ':' in sys.argv[6] else {}
    todo = []; mapng = 0; nogroup = 0; total = 0
    for chars, tpl, slots, lists in shapes(k, kind, lo, hi, caps):
        total += 1
        lists2, mok = last_slot_filter(tpl, slots, lists)
        if not mok: mapng += 1; continue
        plan = plan_groups(len(tpl), slots, [len(l) for l in lists2])
        if plan is None: nogroup += 1; continue
        todo.append((chars, tpl, slots, lists2, plan))
    if mode == 'control':
        n = int(sys.argv[7]) if len(sys.argv) > 7 else 6; rng = random.Random(20261115 + k); ok = 0; tried = 0
        cheap = sorted(todo, key=lambda t: t[4][0])[:max(1, len(todo) * 2 // 3)]
        while tried < n:
            t = rng.choice(cheap); chars, tpl, slots, lists, plan = t
            key = fill(tpl, slots, [rng.choice(l) for l in lists])
            if key_schedule(key)['m'] != 1: continue
            tried += 1; tgt, _ = synth_target(key, rng); J = sum(key_schedule(key)['J']); jlo = max(0, J - rng.randrange(241))
            cands, good, stat, dt = solve(tgt, t, 'c34_%d' % k, jwin=(jlo, jlo + 240))
            hit = key in good; ok += hit
            print('  対照 %d: 語の文字数 %s  鍵 %s  候補 %d 合格 %d  %s  %.1f 秒' % (tried, '・'.join(map(str, chars)), key.decode('cp932'), len(cands), len(good), '出た' if hit else '!!!! 出ない', dt)); sys.stdout.flush()
        print('対照 %d 本中 %d 本で、埋めた鍵が出た' % (n, ok)); sys.exit(0)
    done = {}
    if os.path.exists(LOG):
        for line in open(LOG, encoding='utf-8'):
            try: r = json.loads(line); done[(r['tpl'], tuple(r['lists']), r['kind'])] = r
            except Exception: pass
    tkey = lambda t: (t[1].hex(), tuple(len(l) for l in t[3]), kind)
    fams = {}
    for t in todo:
        s = family_sig(t[1], t[2])
        fams.setdefault((s, tuple(sorted(len(l) for l in t[3]))) if s is not None else ('solo', t[1].hex()), []).append(t)
    work = []
    for mem in fams.values():
        newm = [t for t in mem if tkey(t) not in done]
        if not newm: continue
        donem = [t for t in mem if tkey(t) in done]
        if donem: work.append((0.0, None, newm, done[tkey(donem[0])]))
        else: rep = min(newm, key=lambda t: t[4][0]); work.append((rep[4][0], rep, [t for t in newm if t is not rep], None))
    work.sort(key=lambda w: w[0])
    print('語 %d 個（%s、%d〜%d 文字、上限 %s）: 形 %d、写像で否定 %d、組に分けられない %d、回せる形 %d（済み %d）。等価な形をまとめると %d 組、これから回すのは %d 組、見積もり %.0f 秒' % (
        k, kind, lo, hi, caps or 'なし', total, mapng, nogroup, len(todo), sum(tkey(t) in done for t in todo), len(fams), sum(w[1] is not None for w in work), sum(w[0] for w in work))); sys.stdout.flush()
    if mode == 'count': sys.exit(0)
    tlimit = float(os.environ.get('V34_TLIMIT', '520')); T0 = time.time(); nrun = ncopy = nhit = 0; stopped = False
    def rec_of(t, cands, good, stat, dt, hits, same=None):
        chars, tpl, slots, lists, plan = t
        r = dict(k=k, kind=kind, chars=list(chars), L=len(tpl), tpl=tpl.hex(), lists=[len(l) for l in lists], groups=[plan[1], plan[2], plan[3]],
                 checks=int(stat.get('checks', 0)), cands=cands, sched=good, sec=round(dt, 1))
        if hits: r['hit'] = [c.hex() for c in hits]
        if same: r['same_as'] = same
        return r
    def run_one(f, t):
        cands, good, stat, dt = solve(TARGET_V, t, 'r34_%d' % k)
        hits = [c for c in good if keystream_ok_v(c)]
        r = rec_of(t, len(cands), len(good), stat, dt, hits); f.write(json.dumps(r, ensure_ascii=False) + '\n'); f.flush()
        print('  %-14s 語の数 %s  和が合った組 %d  候補 %d 合格 %d  %.1f 秒%s' % ('・'.join(map(str, t[0])), r['lists'], r['checks'], len(cands), len(good), dt, '  !!!! 当たり ' + ' '.join(c.decode('cp932', 'replace') for c in hits) if hits else '')); sys.stdout.flush()
        return r, bool(hits)
    with open(LOG, 'a', encoding='utf-8') as f:
        for cost, rep, copies, src in work:
            if rep is not None:
                if time.time() - T0 + cost * 0.6 > tlimit and nrun: stopped = True; break
                src, h = run_one(f, rep); nrun += 1; nhit += h
            for t in copies:
                if src.get('sched'): r, h = run_one(f, t); nrun += 1; nhit += h; continue
                r = rec_of(t, src['cands'], src['sched'], dict(checks=src['checks']), 0.0, None, same=src['tpl']); f.write(json.dumps(r, ensure_ascii=False) + '\n'); f.flush(); ncopy += 1
    print('今回 %d 形（回したもの %d、等価な形の記録を写したもの %d）、当たり %d、%.0f 秒%s' % (nrun + ncopy, nrun, ncopy, nhit, time.time() - T0, '。まだ残りがある（もう一度実行すると続きから回る）' if stopped else ''))
