#!/usr/bin/env python3
# 鍵 v が「2 文字の語 6 個を『・』でつないだ 34 バイト」だという仮説を、辞書の語で調べる（docs/password.md §8.10。11-15）。
#   鍵の塊: [語1] [・ 語2の1字目] [語2の2字目 ・] [語3] [・ 語4の1字目] [語4の2字目 ・] [語5] [・ 語6の1字目] ＋ 端数 [語6の2字目]
#   語 1・3・5 は塊 1 個、語 2・4 は塊 2 個、語 6 は塊 1 個＋端数。最後の塊（端数）で写像が決まる。最後以外の塊は、どう並べても同じ内部状態になる。
#   調べ方: k6band（3 個＋3 個の中間一致。帯ごとに貯めて引く）。貯める側 = 語 1・3・5 の組（順不同）、探す側 = 語 2・4（順不同）と語 6。
#   語の一覧: vdict3.load('kk') の 2 文字の語（漢字だけ・カタカナだけ）を頻度順に上位 V 語。
#   python3 k6.py est V                       数の見積もり
#   python3 k6.py control V [本数]            対照（上位 V 語で作った鍵を埋めて、出るかを見る）
#   python3 k6.py run V [秒]                  本番（鍵 v の目標）。続きから回せる。途中の状態は logs/k6_run_V.json、終わったら logs/k6.jsonl に 1 行
import sys, os, time, random, subprocess, json, math
import numpy as np
import vdict, vdict3
from vdict import TARGET_V, keystream_ok_v, synth_target, sched_ok
from gks import chunk_vals, key_schedule, map_of_i2
HERE = vdict.HERE; TMP = vdict.TMP; SEP = '・'.encode('cp932'); LOG = os.path.join(HERE, 'logs', 'k6.jsonl')
DT = np.dtype([('s', '<u8'), ('t', '<u4'), ('pad', '<u4'), ('P', '<f8', (4,))])
NT, DB, NR = int(os.environ.get("K6_NT", 60000)), int(os.environ.get("K6_DB", 19)), 2048          # 帯の組数、表のビット数、区間の数
JWIN = (115, 300)                     # 本番の ΣJ の窓。S と b0b1 の範囲から決まる範囲は 121〜297（k6jwin.py）。少し広げてある

def unit(chunks):
    s = t = 0; P = [1.0] * 4
    for c in chunks:
        f, i, s1, t1 = chunk_vals(c); s += s1; t += t1
        for j in range(4): P[j] *= (1.0 + f[j])
    return s, t & 0xFFFFFFFF, 0, P

def tables(words, want=1):
    """語ごとの単位: W = 塊 1 個（語そのまま）、M = 塊 2 個（・＋1 字目、2 字目＋・）、Z = 塊 1 個＋端数（写像が want になる語だけ）"""
    W = np.zeros(len(words), DT); M = np.zeros(len(words), DT); Zl = []; zi = []
    for n, w in enumerate(words):
        a, b = w[:2], w[2:]; tail = int.from_bytes(b, 'little')
        W[n] = unit([int.from_bytes(w, 'little')])
        M[n] = unit([int.from_bytes(SEP + a, 'little'), int.from_bytes(b + SEP, 'little')])
        if map_of_i2(chunk_vals(tail)[1][1]) == want: Zl.append(unit([int.from_bytes(SEP + a, 'little'), tail])); zi.append(n)
    Z = np.zeros(len(Zl), DT)
    for n, u in enumerate(Zl): Z[n] = u
    return W, M, Z, zi

def run(tgt, words, tabs, jwin, tag='k6', tlimit=0, todo=None, nt=NT, db=DB, nr=NR):
    W, M, Z, zi = tabs; paths = []
    for nm, arr in (('W', W), ('M', M), ('Z', Z)):
        p = os.path.join(TMP, '%s_%s.bin' % (tag, nm)); arr.tofile(p); paths.append((nm, len(arr), p))
    (c1, h1), (c1b, h1b) = tgt['cells1']
    spec = os.path.join(TMP, 'spec_%s.txt' % tag)
    with open(spec, 'w') as f:
        f.write('T %x %.17g %.17g %.17g %.17g %.17g %.17g %d %d 2e-12\nJ %d %d\n' % (tgt['S'], c1, h1, c1b, h1b, tgt['c3'], tgt['h3'], tgt['Tlo'], tgt['Thi'], jwin[0], jwin[1]))
        for nm, n, p in paths: f.write('%s %d %s\n' % (nm, n, p))
        f.write('N %d %d %d\nL %g\n' % (nt, db, nr, tlimit))
        for a, b in (todo or []): f.write('R %d %d\n' % (a, b))
    t0 = time.time(); r = subprocess.run([os.path.join(HERE, 'k6band'), spec], capture_output=True, text=True); dt = time.time() - t0
    if r.returncode: raise RuntimeError('k6band failed: ' + r.stderr[-400:])
    keys = []; stat = {}; done = []; nrange = None
    for line in r.stdout.split('\n'):
        p = line.split()
        if not p: continue
        if p[0] == 'CAND':
            i, j, k, a, b, z = map(int, p[1:7]); keys.append(SEP.join([words[i], words[a], words[j], words[b], words[k], words[zi[z]]]))
        elif p[0] == 'DONE': done.append((int(p[1]), int(p[2]), int(p[3]), int(p[4]), float(p[5])))
        elif p[0] == 'NRANGE': nrange = int(p[1])
        elif p[0] == 'STAT': stat = dict(zip(p[1::2], p[2::2]))
    for nm, n, p in paths: os.remove(p)
    return keys, stat, dt, done, nrange, r.stderr

def same_class(k, key):
    """塊の並べ替えで同じになる鍵か（語 1・3・5 の入れ替え、語 2・4 の入れ替え）"""
    return sorted([k[0:4], k[12:16], k[24:28]]) == sorted([key[0:4], key[12:16], key[24:28]]) and k[28:] == key[28:] and sorted([k[6:10], k[18:22]]) == sorted([key[6:10], key[18:22]])

def intervals(flags):
    out = []; a = None
    for i, d in enumerate(flags + [True]):
        if not d and a is None: a = i
        if d and a is not None: out.append((a, i - 1)); a = None
    return out

if __name__ == '__main__':
    mode = sys.argv[1]; V = int(sys.argv[2]); words = vdict3.load('kk')[2][:V]
    t0 = time.time(); tabs = tables(words); W, M, Z, zi = tabs
    nA = V * (V + 1) * (V + 2) // 6; nB = V * (V + 1) // 2 * len(Z)
    print('2 文字の語 上位 %d（写像 1 になる末尾の語 %d）。貯める側 %.3g 組、探す側 %.3g 組、語の並び %.3g 通り（表の作成 %.1f 秒）' % (V, len(Z), nA, nB, float(V)**5 * len(Z), time.time() - t0)); sys.stdout.flush()
    if mode == 'est': sys.exit(0)
    if mode == 'control':
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 3
        rng = random.Random(20261115 + V); ok = 0
        for t in range(n):
            ws = [rng.choice(words) for _ in range(5)] + [words[rng.choice(zi)]]; key = SEP.join(ws); r = key_schedule(key); assert r['m'] == 1 and len(key) == 34
            tgt, _ = synth_target(key, rng); sj = sum(r['J']); jlo = max(0, sj - rng.randrange(241))
            keys, stat, dt, done, nrange, err = run(tgt, words, tabs, (jlo, jlo + 240), tag='k6c')
            good = [k for k in keys if sched_ok(tgt, k)]; hit = any(same_class(k, key) for k in good); ok += hit
            print('  対照 %d: 鍵 %s  候補 %d 合格 %d  %s  %.1f 秒  区間 %s/%s 帯 %s A %s B %s 引いた回数 %s 粗い一致 %s k7b0 %s k780 %s S %s' % (t + 1, key.decode('cp932'), len(keys), len(good), '出た' if hit else '!!!! 出ない', dt,
                  stat.get('done'), stat.get('ranges'), stat.get('bands'), stat.get('A'), stat.get('B'), stat.get('probe'), stat.get('hit'), stat.get('l3'), stat.get('c1'), stat.get('S'))); sys.stdout.flush()
        print('対照 %d 本中 %d 本で、埋めた鍵（塊の並べ替えを含む）が出た' % (n, ok)); sys.exit(0)
    if mode == 'run':
        tlimit = float(sys.argv[3]) if len(sys.argv) > 3 else 500.0
        sp = os.path.join(HERE, 'logs', 'k6_run_%d.json' % V)
        st = json.load(open(sp)) if os.path.exists(sp) else dict(V=V, nZ=len(Z), nA=nA, nB=nB, jwin=list(JWIN), nt=NT, db=DB, nr=NR, nrange=None, done=None, A=0, B=0, probe=0, sec=0.0, runs=0, cands=[], sched=[], hit=[], finished=False)
        if st['finished']: print('上位 %d 語は終わっている: 候補 %d、鍵スケジュール全体に合うもの %d、鍵ストリームまで合うもの %d' % (V, len(st['cands']), len(st['sched']), len(st['hit']))); sys.exit(0)
        todo = intervals(st['done']) if st['done'] else None
        keys, stat, dt, done, nrange, err = run(TARGET_V, words, tabs, JWIN, tlimit=tlimit, todo=todo)
        if st['nrange'] is None: st['nrange'] = nrange; st['done'] = [False] * nrange
        assert st['nrange'] == nrange
        for q, a, b, p, s in done: st['done'][q] = True; st['A'] += a; st['B'] += b; st['probe'] += p
        st['sec'] += dt; st['runs'] += 1
        for k in keys:
            if k.hex() in st['cands']: continue
            st['cands'].append(k.hex())
            if sched_ok(TARGET_V, k):
                st['sched'].append(k.hex())
                if keystream_ok_v(k): st['hit'].append(k.hex()); print('  !!!! 当たり', k.decode('cp932', 'replace'), k.hex())
        nd = sum(st['done']); st['finished'] = nd == nrange
        json.dump(st, open(sp, 'w'))
        print('上位 %d 語: 区間 %d/%d 済み（今回 %d、%.0f 秒。通算 %.0f 秒）。A %.4g B %.4g 引いた回数 %.4g。候補 %d、鍵スケジュール全体に合うもの %d、鍵ストリームまで合うもの %d。今回: %s' % (
            V, nd, nrange, len(done), dt, st['sec'], st['A'], st['B'], st['probe'], len(st['cands']), len(st['sched']), len(st['hit']), ' '.join('%s=%s' % kv for kv in stat.items())))
        if st['finished']:
            rec = {k: st[k] for k in ('V', 'nZ', 'nA', 'nB', 'jwin', 'nt', 'db', 'nr', 'nrange', 'A', 'B', 'probe', 'runs')}; rec.update(sec=round(st['sec'], 1), cands=len(st['cands']), sched=len(st['sched']), hit=st['hit'])
            with open(LOG, 'a', encoding='utf-8') as f: f.write(json.dumps(rec, ensure_ascii=False) + '\n')
            print('上位 %d 語は全部終わった。' % V)
