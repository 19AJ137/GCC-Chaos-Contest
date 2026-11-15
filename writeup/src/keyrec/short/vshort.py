# 短い等価鍵を作る: 小数部が 4 つとも大きい塊（short/hi90.txt）を固定し、残りの 3 塊を入口（vres）で解く。
#   先に塊の一覧を作る: ./short/hiscan 0.90 > short/hi90.txt（2 コアで 4 秒、36 MB）
#   python3 short/vshort.py <長さ 21..28> [仮説の数 上限] [制限秒] [格子 16|17] [乱数の種] [鍵の組がこの数出たら止める]
# 記録: logs/vshort.jsonl（1 回の束ごとに 1 行: 長さ・仮説の数・秒・出た鍵・解いた仮説の既知の部分）。まとめた logs/vshort.jsonl.gz があれば、済みの仮説としてそれも読む
import sys, os, json, time, math, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from vpart import build_hyp, run_vres, TARGET_V, keystream_ok_v, textlike
from gks import key_schedule, chunk_vals, map_of_i2
SIGMA = 16582419900                      # σ（塊の値の和）の候補の中ほど
A_U, B_U = 0x3243F6A8, 0xFFFFFFFF        # 格子で解く塊の値の範囲（格子 16）

def read_log(log):
    """記録（.jsonl と、まとめた .jsonl.gz があればそれも）を 1 行ずつ返す"""
    import gzip
    for fn, op in ((log + '.gz', gzip.open), (log, open)):
        if os.path.exists(fn):
            with op(fn, 'rt', encoding='utf-8') as f:
                for line in f:
                    if line.strip(): yield json.loads(line)
def load_hi():
    hi = []
    for line in open(os.path.join(HERE, 'hi90.txt')):
        p = line.split(); hi.append((int(p[0], 16), tuple(float(x) for x in p[1:5]), int(p[5])))
    hi.sort(key=lambda x: -sum(math.log2(1 + v) for v in x[1]))
    return hi
def last_candidates(nb, hi):
    """端数の塊（nb バイト）で写像 1 になるもののうち、(1+f) の積が大きいもの"""
    if nb == 3: return [c for c, f, m in hi if (1 << 16) <= c < (1 << 24) and m == 1][:24]
    out = []
    for c in range(1 << (8 * (nb - 1)) if nb > 1 else 0, 1 << (8 * nb)):
        f, i, s, t = chunk_vals(c)
        if map_of_i2(i[1]) == 1: out.append((-sum(math.log2(1 + v) for v in f), c))
    out.sort(); return [c for _, c in out[:(2 if nb == 1 else 6)]]
def gen(L, hi, rng, ntop=1500):
    """長さ L の仮説を 1 つ作る。戻り値: (既知のバイト列, 未知 12 バイトを前に置いた segs)"""
    nfull = L // 4; nb = L % 4; nk = nfull - 3
    full = [x for x in hi[:ntop * 4] if x[0] >= 0x01000000][:ntop]
    while True:
        if nb:
            last = rng.choice(gen.last[nb]); ks = [rng.choice(full)[0] for _ in range(nk)]; lastb = last.to_bytes(nb, 'little')
        else:
            m1 = [x for x in full if x[2] == 1]; ks = [rng.choice(full)[0] for _ in range(nk - 1)]; last = rng.choice(m1)[0]; lastb = last.to_bytes(4, 'little')
        if len(set(ks)) != len(ks) or (not nb and last in ks): continue
        su = SIGMA - sum(ks) - last
        # 未知 3 塊の和が、取りうる範囲の中ほど（端から 15% 以上内側）にあるものだけ
        lo, hi_ = 3 * A_U, 3 * B_U
        if not (lo + 0.15 * (hi_ - lo) <= su <= hi_ - 0.15 * (hi_ - lo)): continue
        known = b''.join(c.to_bytes(4, 'little') for c in sorted(ks)) + lastb
        return known
gen.last = {}

if __name__ == '__main__':
    L = int(sys.argv[1]); nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 2000; tlimit = float(sys.argv[3]) if len(sys.argv) > 3 else 300.0
    lat = int(sys.argv[4]) if len(sys.argv) > 4 else 16; seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    want = int(sys.argv[6]) if len(sys.argv) > 6 else 1          # 鍵の組がこの数だけ出たら止める
    hi = load_hi(); rng = random.Random(seed * 1000 + L)
    for nb in (1, 2, 3): gen.last[nb] = last_candidates(nb, hi)
    log = os.path.join(os.path.dirname(HERE), 'logs', 'vshort.jsonl'); os.makedirs(os.path.dirname(log), exist_ok=True)
    done = set()
    for r in read_log(log):
        if r['L'] == L and r['lat'] == lat: done.update(r.get('known', []))
    T0 = time.time(); nh = 0; sets = set(); allkeys = []; batch_n = 400 if L >= 25 else 2000
    print('長さ %d、格子 %d、済みの仮説 %d' % (L, lat, len(done)), flush=True)
    while nh < nmax and time.time() - T0 < tlimit and len(sets) < want:
        hyps = []; kn = []
        while len(hyps) < min(batch_n, nmax - nh):
            known = gen(L, hi, rng)
            if known.hex() in done: continue
            done.add(known.hex()); h = build_hyp('h%06d' % len(hyps), [12, known], lat=lat)
            if not h['feasible'] or h['nU'] != 3: continue
            hyps.append(h); kn.append(known.hex())
        keys, stat, dt = run_vres(TARGET_V, hyps, L=lat, tag='vshort%d' % L)
        nh += len(hyps); good = []
        for hid, ks in keys.items():
            for k in sorted(ks):
                r = key_schedule(k); ok = keystream_ok_v(k)
                if ok:
                    good.append(k.hex()); allkeys.append(k)
                    sets.add(tuple(sorted(int.from_bytes(k[i:i + 4], 'little') for i in range(0, 12, 4))) + (k[12:],))
        with open(log, 'a') as f: f.write(json.dumps(dict(L=L, lat=lat, n=len(hyps), sec=round(dt, 2), keys=good, known=kn)) + '\n')
        print('  仮説 %d 個 %.1f 秒（1 個 %.4f 秒）、鍵ストリームが合う鍵 %d 本（累計 %d 組）、経過 %.0f 秒' % (len(hyps), dt, dt / len(hyps), len(good), len(sets), time.time() - T0), flush=True)
    print('計: 仮説 %d 個、%.0f 秒、鍵の組 %d、鍵 %d 本' % (nh, time.time() - T0, len(sets), len(allkeys)))
    for k in allkeys[:12]:
        r = key_schedule(k); print('   %s  J %s ΣJ %d b0b1 %04x 区分 %d' % (k.hex(), r['J'], sum(r['J']), r['b0b1'], textlike(k)))
