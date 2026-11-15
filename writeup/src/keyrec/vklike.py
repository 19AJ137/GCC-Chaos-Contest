# 鍵 v が「鍵 k と同じ形式」だという仮説（11-14、3 回目）:
#   2〜3 文字の語（全角。中身は未知＝文字の総当り）と、区切り文字と、IISI を、順不同に並べたもの。
#   語の中身は辞書を使わず、未知のバイト（0x20〜0xFC の全部）として入口（vres）に解かせる。
#   python3 vklike.py count [語の数…]      形の数と、入口に載るかどうか（見積もり）
#   python3 vklike.py run   語の数 [予算秒]  載る形を解く（記録は logs/vcamp_KL<語の数>.jsonl）
#   python3 vklike.py sigma                 塊の数ごとに、S と b0b1 の両方に合う「塊の値の和」の候補があるか
import sys, json, time, itertools
from vpart import build_hyp, run_vres, TARGET_V, keystream_ok_v, show, textlike, segs_to_known
from vcamp import last_chunk_map, pat, load_done

IISI = b'IISI'
SEP_TEXT = ['・', '', ' ', '-', '_', '.', ',', '/', '+', '&', '@', ':', ';', '=', '*', '#', '!', '|', '~',
            '　', '、', '。', '，', '．', '：', '／', '－', '―', '‐', '＿', '＋', '＆', '＠', '＝', '＊', '～', 'ー', '×', '☆', '★', '→', '＃']
SEPS = [s.encode('cp932') for s in SEP_TEXT] + [b'\xa5']          # 最後は半角カナの中点
WORD_BYTES = (4, 6)                                                 # 全角 2 文字・3 文字

def forms(nw, seps=SEPS):
    """語 nw 個と IISI を並べた形。区切りは、隣り合う要素の間に同じものを入れる。戻り値: (segs, 説明)"""
    for sizes in itertools.product(WORD_BYTES, repeat=nw):
        for pos in range(nw + 1):                                   # IISI の位置
            for sep in seps:
                elems = [s for s in sizes]; elems.insert(pos, IISI)
                segs = []
                for i, e in enumerate(elems):
                    if i and sep: segs.append(sep)
                    segs.append(e)
                # 既知どうしをつなぐ
                out = []
                for s in segs:
                    if out and isinstance(s, bytes) and isinstance(out[-1], bytes): out[-1] = out[-1] + s
                    else: out.append(s)
                yield out, (sizes, pos, sep)

def classify(nw, budget):
    run = []; skip = []; mapng = 0; seen = set(); total = 0
    for segs, info in forms(nw):
        known = tuple(segs_to_known(segs))
        if known in seen: continue
        seen.add(known); total += 1
        m = last_chunk_map(known)
        if m is not None and m != TARGET_V['m']: mapng += 1; continue
        hid = 'KL%d_%05d' % (nw, len(run) + len(skip))
        h17 = build_hyp(hid, segs, lat=17); h16 = build_hyp(hid, segs, lat=16)
        if h17['cost'] <= budget: run.append((17, h17))
        elif h16['cost'] <= budget: run.append((16, h16))
        else: skip.append(h16 if h16['cost'] <= h17['cost'] else h17)
    return total, mapng, run, skip

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'sigma':
        from stfilter import sigma_candidates, TV
        for n in range(1, 12):
            c = sigma_candidates(TV['S'], TV['Tlo'], TV['Thi'], 0, 0, n, n, cmin=0x20000000)
            print('塊 %2d 個（%2d〜%2d バイト）: 和の候補 %4d 個' % (n, 4 * n - 3, 4 * n, len(set(sg for sg, j3, w in c))))
    elif mode == 'count':
        budget = 600.0
        for nw in [int(x) for x in sys.argv[2:]] or [2, 3, 4]:
            total, mapng, run, skip = classify(nw, budget)
            byL = {}
            for lat, h in run: byL.setdefault(h['L'], [0, 0, 0.0]); byL[h['L']][0] += 1; byL[h['L']][2] += h['cost']
            for h in skip: byL.setdefault(h['L'], [0, 0, 0.0]); byL[h['L']][1] += 1
            unk = sorted(set(sum(x is None for x in h['known']) for lat, h in run)); unks = sorted(set(sum(x is None for x in h['known']) for h in skip))
            print('語 %d 個: 形 %d、写像で否定 %d、入口に載る %d（見積もり %.0f 秒。未知 %s バイト）、載らない %d（未知 %s バイト）' % (
                nw, total, mapng, len(run), sum(h['cost'] for lat, h in run), unk, len(skip), unks))
            print('   長さごと（載る/載らない）: ' + '  '.join('%d: %d/%d' % (L, v[0], v[1]) for L, v in sorted(byL.items())))
            if skip: print('   載らない形のうち、いちばん軽いものの見積もり: %.3g 秒  例 %s' % (min(h['cost'] for h in skip), show(min(skip, key=lambda h: h['cost']))))
    elif mode == 'run':
        nw = int(sys.argv[2]); budget = float(sys.argv[3]) if len(sys.argv) > 3 else 600.0; tlimit = float(sys.argv[4]) if len(sys.argv) > 4 else 540.0
        total, mapng, run, skip = classify(nw, budget)
        fam = 'KL%d' % nw; log = 'logs/vcamp_%s.jsonl' % fam; done = load_done()
        todo = [(lat, h) for lat, h in run if done.get(pat(h), 0) < lat]; todo.sort(key=lambda x: (x[0], x[1]['cost']))
        print('%s: 形 %d、写像で否定 %d、載る %d（済み %d）、今回 %d（見積もり %.0f 秒）、載らない %d' % (fam, total, mapng, len(run), len(run) - len(todo), len(todo), sum(h['cost'] for lat, h in todo), len(skip))); sys.stdout.flush()
        T0 = time.time(); pos = 0; nhit = 0; ncand = 0
        with open(log, 'a') as f:
            while pos < len(todo) and time.time() - T0 < tlimit:
                lat = todo[pos][0]; batch = []; est = 0.0
                while pos < len(todo) and todo[pos][0] == lat and (not batch or (est + todo[pos][1]['cost'] <= 40.0 and len(batch) < 3000)): batch.append(todo[pos][1]); est += todo[pos][1]['cost']; pos += 1
                keys, stat, dt = run_vres(TARGET_V, batch, L=lat, tag='camp_%s' % fam)
                for h in batch:
                    ks = sorted(keys.get(h['id'], set())); good = [k.hex() for k in ks if keystream_ok_v(k)]; sst = stat.get(h['id'], {})
                    rec = dict(fam=fam, lat=lat, pat=pat(h), known=show(h), L=h['L'], nU=h['nU'], E=h['E'], n=sst.get('assignments'), raw=sst.get('rawkeys', 0))
                    if ks: rec['cand'] = [k.hex() for k in ks]; ncand += 1
                    if good: rec['hit'] = good; rec['text'] = max(textlike(bytes.fromhex(g)) for g in good); nhit += 1; print('  鍵ストリームが合う鍵 %d 本  %s  区分 %d' % (len(good), show(h), rec['text']))
                    f.write(json.dumps(rec, ensure_ascii=False) + '\n')
                f.flush(); print('  %s 格子 %d  %d / %d  %.0f 秒（この回 %.1f 秒、見積もり %.1f 秒）' % (fam, lat, pos, len(todo), time.time() - T0, dt, est)); sys.stdout.flush()
        print('%s 終わり: 今回 %d 件、当たり %d 件、候補だけ %d 件、残り %d 件、%.0f 秒' % (fam, pos, nhit, ncand, len(todo) - pos, time.time() - T0))
