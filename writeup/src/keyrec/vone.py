# 個別の仮説を、格子を指定して解く。記録は logs/vcamp_C.jsonl などに足す（vcamp と同じ形式）。
#   python3 vone.py 族 格子 '既知:未知の並び' …    並びは「文字列」と「<バイト数>」を | で区切る。例: 'GoodLuckToIISIGoodLuckTo|<11>'
import sys, json, time
from vpart import build_hyp, run_vres, TARGET_V, keystream_ok_v, show, textlike
from vcamp import pat
fam = sys.argv[1]; lat = int(sys.argv[2])
for spec in sys.argv[3:]:
    segs = [int(x[1:-1]) if x.startswith('<') and x.endswith('>') else x.encode('cp932') for x in spec.split('|')]
    h = build_hyp('one', segs, lat=lat)
    print('%s  格子 %d  解き方 U%d E%d  見積もり %.0f 秒' % (show(h), lat, h['nU'], h['E'], h['cost'])); sys.stdout.flush()
    keys, stat, dt = run_vres(TARGET_V, [h], L=lat, tag='one_%s' % fam)
    ks = sorted(keys.get('one', set())); good = [k.hex() for k in ks if keystream_ok_v(k)]; sst = stat.get('one', {})
    rec = dict(fam=fam, lat=lat, pat=pat(h), known=show(h), L=h['L'], nU=h['nU'], E=h['E'], n=sst.get('assignments'), raw=sst.get('rawkeys', 0))
    if ks: rec['cand'] = [k.hex() for k in ks]
    if good: rec['hit'] = good; rec['text'] = max(textlike(bytes.fromhex(g)) for g in good)
    with open('logs/vcamp_%s.jsonl' % fam, 'a') as f: f.write(json.dumps(rec, ensure_ascii=False) + '\n')
    lv = [textlike(bytes.fromhex(g)) for g in good]
    print('   %.0f 秒  鍵ストリームが合う鍵 %d 本（区分 %s）' % (dt, len(good), {l: lv.count(l) for l in sorted(set(lv))})); sys.stdout.flush()
    for g in good:
        if textlike(bytes.fromhex(g)) >= 1: print('     区分 %d  %s  %s' % (textlike(bytes.fromhex(g)), g, bytes.fromhex(g).decode('cp932', 'replace')))
