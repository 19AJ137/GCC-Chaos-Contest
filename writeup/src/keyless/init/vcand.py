"""vcand.py — 鍵 v（写像 1）の初期値 y_init の候補を全部出す。2026-11-12
  語 1 の状態 y_1（lo ごとに確定）→ 語 0 の状態 y_0 の候補（符号と下位ビット）→ K'_0 ごとに IV_v = S_v → 空回し warm = max(11, S_v mod 79)
  → y_0 から warm+1 歩さかのぼった格子点で |y| < 1 のものが y_init の候補。
  出力: vcand.pkl = [(lo, S_v, warm, 符号, idx の配列)]"""
import sys, os, pickle, time, importlib.util
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
import back1 as B1
spec = importlib.util.spec_from_file_location('state0', os.path.join(D, 'state0.py')); V = importlib.util.module_from_spec(spec); spec.loader.exec_module(V)
A = open(os.path.join(CT, 'A.bin'), 'rb').read()[:1024]; Aw = V.words(A)
K, F = V.stream_k(); KP, FP = V.stream_v()
M0 = int.from_bytes(b'Is t', 'big')
p0 = V.u2d(V.V_P); y0v = V.u2d(V.V_Y0); y1 = V.setlow(y0v * y0v - p0, V.V_LO)
def le(iv): return int.from_bytes(iv.to_bytes(4, 'big'), 'little')
def stream_from(y, p, lo, n):
    ks = []; fl = []
    for i in range(n):
        w, f = V.word((1.0 + y) * 0.5); ks.append(w); fl.append(f); y = V.setlow(y * y - p, lo)
    return ks, fl
LO_A, LO_B = 0xab0f, 0xabcf
if __name__ == '__main__':
    t0 = time.time(); out = []; tot = 0; over = []; stat = {}
    los = range(LO_A, LO_B + 1)
    # 検算: どの lo でも語 1 以降の鍵ストリームが同じか
    for lo in (LO_A, V.V_LO, LO_B):
        p = V.u2d((V.d2u(p0) & ~0xFFFF) | lo); y = V.u2d((V.d2u(y1) & ~0xFFFF) | lo); ks, fl = stream_from(y, p, lo, 255)
        print('lo = %04x: 語 1〜255 の鍵ストリームが verify.py と %s' % (lo, '一致' if ks == KP[1:] and fl == FP[1:] else '不一致'))
    for lo in los:
        p = V.u2d((V.d2u(p0) & ~0xFFFF) | lo); y1l = V.u2d((V.d2u(y1) & ~0xFFFF) | lo)
        sg1, k1 = B1.tkey(y1l); r = B1.back_n([(sg1, k1, k1)], p, lo, 1)
        for sg, i1, i2 in zip(*r):
            idx = np.arange(i1, i2 + 1, dtype=np.int64); y = sg * ((idx.astype(np.uint64) << np.uint64(16)) | np.uint64(lo)).view(np.float64)
            x = (1.0 + y) * 0.5; b = x.view(np.uint8).reshape(-1, 8)
            w = (b[:, 4].astype(np.uint32) << 24) | (b[:, 2].astype(np.uint32) << 16) | (b[:, 3].astype(np.uint32) << 8) | b[:, 5].astype(np.uint32)
            assert (b[:, 6] <= 0xE6).all()                      # 語 0 は鍵 v で連鎖あり
            cut = np.nonzero(np.diff(w.astype(np.int64)) != 0)[0] + 1; st = np.concatenate([[0], cut]); en = np.concatenate([cut, [len(w)]])
            for a, e in zip(st, en):
                kw = int(w[a]); iv = Aw[0] ^ M0 ^ K[0] ^ kw; Sv = le(iv); warm = max(11, Sv % 79)
                res = B1.back_n([(int(sg), int(idx[a]), int(idx[e - 1]))], p, lo, warm + 1)
                key = (int(sg), kw)
                if res is None: over.append((lo, kw, warm)); continue
                s2, j1, j2 = res
                if len(s2) == 0: continue
                # |y| < 1 だけ（y_init = 2·k780 − 1、0 <= k780 < 1）。idx が 1.0 の idx 未満
                one = V.d2u(1.0) >> 16; j2 = np.minimum(j2, one - 1); ok = j1 <= j2; s2, j1, j2 = s2[ok], j1[ok], j2[ok]
                # y = −1 ちょうど（k780 = 0）は除いてよい（積の小数部が 0 になるのは事実上ない）
                n = int((j2 - j1 + 1).sum())
                if n == 0: continue
                if n > 20000000: over.append((lo, kw, warm, n)); continue
                pts = np.concatenate([np.arange(a1, a2 + 1, dtype=np.int64) * s for s, a1, a2 in zip(s2, j1, j2)]) if len(s2) < 200000 else None
                if pts is None:
                    rep = (j2 - j1 + 1); base = np.repeat(j1, rep); off = np.arange(n) - np.repeat(np.cumsum(rep) - rep, rep); pts = (base + off) * np.repeat(s2, rep)
                out.append((lo, Sv, warm, kw, pts)); tot += n; stat[key] = stat.get(key, 0) + n
        if (lo - LO_A) % 16 == 15: print('  lo %04x まで: 候補 %d 点、%.0f 秒' % (lo, tot, time.time() - t0), flush=True)
    print('y_init の候補: %d 点（lo × 語 0 の候補の組 %d 個）。上限を超えて飛ばした組 %d: %s、%.0f 秒' % (tot, len(out), len(over), over[:10], time.time() - t0))
    print('語 0 の鍵ストリームの候補ごとの点の数:')
    for (sg, kw), n in sorted(stat.items()): print('   符号 %+d K\'_0 = %08x: %d' % (sg, kw, n))
    pickle.dump(out, open(os.path.join(D, 'vcand.pkl'), 'wb'))
