# 打てる文字だけでできた短い等価鍵を探す。固定する塊を「4 バイトだけで打てる文字になるもの」から選び、解けた 3 塊と合わせて、
# 塊の並べ替え（最後の塊以外は入れ替えても同じ内部状態）のどれかが Shift_JIS の文字列として打てるものを拾う。
#   python3 short/vshort_typ.py <長さ> [制限秒] [乱数の種] [打てる鍵がこの数出たら止める]
# 記録: logs/vshort_typ.jsonl（まとめた logs/vshort_typ.jsonl.gz があれば、済みの仮説としてそれも読む）
import sys, os, json, time, math, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from vpart import build_hyp, run_vres, TARGET_V, keystream_ok_v, textlike
from gks import key_schedule, chunk_vals, map_of_i2
from vshort import load_hi, read_log, SIGMA, A_U, B_U
b4 = lambda c: c.to_bytes(4, 'little')
D = B_U - A_U
def nice(key):
    """区分（vpart.textlike）と、空白を含まないか"""
    return textlike(key), (0x20 not in key)
def best_order(chunks, last):
    best = None
    for p in itertools.permutations(chunks):
        k = b''.join(b4(c) for c in p) + last; t = nice(k)
        if best is None or t > best[0]: best = (t, k)
    return best

if __name__ == '__main__':
    L = int(sys.argv[1]); tlimit = float(sys.argv[2]) if len(sys.argv) > 2 else 540.0; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    want = int(sys.argv[4]) if len(sys.argv) > 4 else 1; lat = 16
    nfull = L // 4; nb = L % 4; nk = nfull - 3
    hi = load_hi(); rng = random.Random(seed * 7919 + L)
    typ = [x for x in hi if min(x[1]) >= 0.95 and textlike(b4(x[0])) >= 2]
    typ1 = [x for x in typ if x[2] == 1]
    if nb == 1: lasts = [b'Q', b'j']
    elif nb == 2: lasts = [bytes.fromhex('4b46')]
    elif nb == 3: lasts = [c.to_bytes(3, 'little') for c, f, m in hi if (1 << 16) <= c < (1 << 24) and m == 1 and textlike(c.to_bytes(3, 'little')) >= 2][:8]
    else: lasts = None
    # 解く 3 塊の和の位置 x（0〜3）。21〜23 バイトは 1.58 以上にしかならない。打てる割合と密度の積が大きい範囲に限る
    xlo, xhi = (1.55, 1.95) if nfull == 5 else (1.0, 2.0)
    log = os.path.join(os.path.dirname(HERE), 'logs', 'vshort_typ.jsonl'); done = set(); nsets0 = 0
    for r in read_log(log):
        if r['L'] == L: done.update(r.get('known', [])); nsets0 += len(r.get('sets', []))
    print('長さ %d: 固定に使う打てる塊 %d 個（写像 1 は %d 個）、端数 %s、済みの仮説 %d（鍵の組 %d）' % (L, len(typ), len(typ1), [x.hex() for x in lasts] if lasts else '-', len(done), nsets0), flush=True)
    T0 = time.time(); nh = 0; nsets = 0; found = []
    ordered = None
    if nfull == 5 and nb == 1:
        # 21 バイト: 固定する 2 塊と端数の組を、当たりやすさの見積もりの順に使う（short/vtyp_pool.py）
        from vtyp_pool import ranked_pool
        ordered, _ = ranked_pool(hi, done, log=lambda *a: print(*a, flush=True)); print('  使える組 %d 通り（済みを除く。当たりやすい順）' % len(ordered), flush=True)
    opos = 0
    while time.time() - T0 < tlimit and len(found) < want:
        hyps = []; kn = []; tries = 0
        while ordered is not None and len(hyps) < 2000 and opos < len(ordered):
            known = ordered[opos]; opos += 1; done.add(known.hex()); h = build_hyp('h%06d' % len(hyps), [12, known], lat=lat)
            if not h['feasible'] or h['nU'] != 3: continue
            hyps.append(h); kn.append(known.hex())
        if ordered is not None and not hyps: print('組を使い切った'); break
        while ordered is None and len(hyps) < (400 if nfull >= 6 and nb else 2000 if nfull == 5 else 1000):
            tries += 1
            if nb: last = rng.choice(lasts); ks = [rng.choice(typ)[0] for _ in range(nk)]; lastv = int.from_bytes(last, 'little')
            else: ks = [rng.choice(typ)[0] for _ in range(nk - 1)]; lastv = rng.choice(typ1)[0]; last = b4(lastv)
            if len(set(ks)) != len(ks) or (not nb and lastv in ks): continue
            su = SIGMA - sum(ks) - lastv; x = (su - 3 * A_U) / D
            if not (xlo <= x <= xhi): continue
            known = b''.join(b4(c) for c in sorted(ks)) + last
            if known.hex() in done: continue
            done.add(known.hex()); h = build_hyp('h%06d' % len(hyps), [12, known], lat=lat)
            if not h['feasible'] or h['nU'] != 3: continue
            hyps.append(h); kn.append(known.hex())
        keys, stat, dt = run_vres(TARGET_V, hyps, L=lat, tag='vshorttyp%d' % L)
        nh += len(hyps); sets = {}
        for hid, ks_ in keys.items():
            for k in ks_:
                if not keystream_ok_v(k): continue
                full = [int.from_bytes(k[i:i + 4], 'little') for i in range(0, 4 * nfull, 4)]; last = k[4 * nfull:]
                if not nb: last = b4(full[-1]); full = full[:-1]
                sets[(tuple(sorted(full)), last)] = 1
        rec = []
        for (full, last) in sets:
            (t, nosp), k = best_order(list(full), last); nsets += 1
            assert keystream_ok_v(k)
            rec.append(dict(key=k.hex(), lv=t, nosp=nosp))
            if t >= 2:
                print('   *** 打てる鍵（区分 %d%s）: %s  %s' % (t, '' if nosp else '、空白あり', k.hex(), k.decode('cp932')), flush=True)
                if nosp: found.append(k)
        with open(log, 'a') as f: f.write(json.dumps(dict(L=L, n=len(hyps), sec=round(dt, 2), sets=rec, known=kn)) + '\n')
        print('  仮説 %d 個 %.1f 秒、鍵の組 %d（累計 %d、うち区分 1 以上 %d）、打てる鍵 累計 %d、経過 %.0f 秒' % (len(hyps), dt, len(sets), nsets, sum(1 for r in rec if r['lv'] >= 1), len(found), time.time() - T0), flush=True)
    print('計: 仮説 %d 個、%.0f 秒、鍵の組 %d、打てる鍵 %d 本' % (nh, time.time() - T0, nsets, len(found)))
