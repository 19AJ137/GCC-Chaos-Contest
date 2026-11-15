"""mkstate.py — 2 つの鍵の「鍵スケジュールの出力」（DLL が鍵から作る内部の値の一式）を、倍精度の具体的な値として 1 組ずつ作る。2026-11-12
  1 ブロック目と 2 ブロック目の定数（上位 48 ビット）の両方を再現する k7b0・k7b8 を、近似値のまわりで探す。
  出力: state.json（verify.py に埋め込む値）"""
import sys, os, pickle, struct, json, math
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
import eblock2 as E2M, eblock2x as X2
V = E2M.V; gcc = E2M.gcc; u2d, d2u, setlow = V.u2d, V.d2u, V.setlow
LOK = V.K_LO; LOV = V.V_LO
ex = pickle.load(open(os.path.join(D, 'eblock2_exact.pkl'), 'rb')); s2 = [s for s in ex['sols'] if s['lo'] == LOK and s['u0'] > 0][0]
found = pickle.load(open(os.path.join(D, 'eblock2_found.pkl'), 'rb')); kp, vp, p48c = found['kp'], found['vp'], found['p48c']
BK = (E2M.RK0 >> 16) & 0xFF; BV = (E2M.RV0 >> 16) & 0xFF
def up(x): return d2u(x) >> 16
def kconst(k7b0, k7b8):
    P = k7b0 * gcc.C_0_4 + gcc.C_0_9; Q = (k7b8 * gcc.C_M014 - P * gcc.C_0_8) + gcc.C_1_72; A = Q * P + 1.0; C = 1.0 - Q
    return up(Q), up(A), up(C)
T1 = (V.K_Q >> 16, V.K_A >> 16, V.K_C >> 16); T2 = (up(s2['Q']), up(s2['A']), up(s2['C']))
# 鍵 k: k7b0・k7b8 を格子で探す
good = []; c0 = kp['k7b0']; c8 = kp['k7b8']
for i in range(-400, 401):
    k7b0 = c0 + i * 5e-13
    for j in range(-400, 401):
        k7b8 = c8 + j * 5e-13
        if kconst(k7b0, k7b8) != T1: continue
        if kconst(E2M.repl5(k7b0, BK)[0], k7b8) != T2: continue
        good.append((i, j, k7b0, k7b8))
print('鍵 k: 両方のブロックの定数を再現する (k7b0, k7b8) の格子点 %d 個' % len(good))
ci = sum(g[0] for g in good) / len(good); cj = sum(g[1] for g in good) / len(good); g = min(good, key=lambda g: (g[0] - ci) ** 2 + (g[1] - cj) ** 2)
k7b0k, k7b8k = g[2], g[3]; print('  採用: k7b0 = %.17g（%016x）、k7b8 = %.17g（%016x）。範囲: k7b0 %.3e 幅、k7b8 %.3e 幅' % (k7b0k, d2u(k7b0k), k7b8k, d2u(k7b8k), (max(x[0] for x in good) - min(x[0] for x in good)) * 5e-13, (max(x[1] for x in good) - min(x[1] for x in good)) * 5e-13))
# 鍵 k の初期値（k780・k788）: kinit.pkl の lo = LOK の値
ki = [x for x in pickle.load(open(os.path.join(D, 'kinit.pkl'), 'rb')) if x[0] == LOK][0]; k780k, k788k = ki[1], ki[2]
# 鍵 v: k7b0
goodv = []; cv = vp['k7b0']
for i in range(-4000, 4001):
    k7b0 = cv + i * 1e-13
    if up(k7b0 * gcc.C_TENTH2 + gcc.C_1_9) != (V.V_P >> 16): continue
    if up(E2M.repl5(k7b0, BV)[0] * gcc.C_TENTH2 + gcc.C_1_9) != p48c + 1: continue
    goodv.append(k7b0)
print('鍵 v: 両方のブロックの p を再現する k7b0 の格子点 %d 個（幅 %.3e）' % (len(goodv), goodv[-1] - goodv[0])); k7b0v = goodv[len(goodv) // 2]
# 鍵 v の初期値: y_init（lo = LOV の候補のうち小さい方の idx）→ k780 = (1 + y)/2
vc = [c for c in pickle.load(open(os.path.join(D, 'vcand.pkl'), 'rb')) if c[0] == LOV][0]; idxs = sorted(set(abs(int(i)) for i in vc[4])); Sv = vc[1]
yc = u2d((idxs[0] << 16) | 0x8000); k780v = (1.0 + yc) * 0.5; assert up(k780v * 2.0 - 1.0) == idxs[0], (up(k780v * 2.0 - 1.0), idxs)
print('  y_init の idx の候補 %s（y ≈ ±%.15f）、S_v = %08x、k780 = %.17g（もう一方の符号なら %.17g）、k7b0 = %.17g' % (idxs, yc, Sv, k780v, (1.0 - yc) * 0.5, k7b0v))
Sk = int.from_bytes(V.K_IV.to_bytes(4, 'big'), 'little')
st = dict(k=dict(mode=2, k780=d2u(k780k), k788=d2u(k788k), k7b0=d2u(k7b0k), k7b8=d2u(k7b8k), S=Sk, lo=LOK), v=dict(mode=1, k780=d2u(k780v), k7b0=d2u(k7b0v), S=Sv, lo=LOV))
json.dump(st, open(os.path.join(D, 'state.json'), 'w'), indent=1)
print('K = dict(mode=2, k780=0x%016x, k788=0x%016x, k7b0=0x%016x, k7b8=0x%016x, S=0x%08x, lo=0x%04x)' % (d2u(k780k), d2u(k788k), d2u(k7b0k), d2u(k7b8k), Sk, LOK))
print('V = dict(mode=1, k780=0x%016x, k7b0=0x%016x, S=0x%08x, lo=0x%04x)' % (d2u(k780v), d2u(k7b0v), Sv, LOV))
print('値: k: k780 = %.15f k788 = %.15f k7b0 = %.15f k7b8 = %.15f; v: k780 = %.15f k7b0 = %.15f' % (k780k, k788k, k7b0k, k7b8k, k780v, k7b0v))
