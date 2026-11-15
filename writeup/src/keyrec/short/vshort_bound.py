# 鍵 v と等価な鍵の長さの下限。S と b0b1 の 2 条件だけを使い、塊の値は 0〜2^32-1 の全部（端数の塊・小さい塊を含む）を対象にする。
#   python3 short/vshort_bound.py
# 考え方:
#   S + w·2^32 = Σi1 + R + ΣJ、 T ≡ Σi3 + Q (mod 2^16)。 R = Σ ftol(i3/i2)、Q = Σ ftol(i1/i4)。
#   K_PI·σ = Σi1 + F1 + E1、 K_E·σ = Σi3 + F3 + E3 （F = 小数部の和、E = clr と丸めの切り捨ての和）
#   → T16 ≡ ρ·(S_w − Y) + Z (mod 65536)、 ρ = K_E/K_PI、 Y = R + ΣJ − F1 − E1、 Z = Q − F3 − E3
#   Y と Z の取りうる範囲から、T16 が [Tlo, Thi] に入りうるかを判定する。入りえなければ、その塊の数の鍵は存在しない。
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fractions import Fraction as Fr
from gks import K_PI, K_E, chunk_vals, map_of_i2
S, Tlo, Thi = 0x3a9d3481, 0xab0f, 0xabc5
KP, KE = Fr(K_PI), Fr(K_E); RHO = KE / KP; M16 = 65536

# ---- 小さい塊（例外）: ftol(i3/i2)・ftol(i1/i4) が 3・2 にならない c
def rq(c):
    f, i, s, t = chunk_vals(c); return (s - int(i[0])) & 0xFFFFFFFF, (t - int(i[2])) & 0xFFFFFFFF
exo = [c for c in range(0, 1 << 16) if rq(c) != (3, 2)]
CX = max(exo)
print('3・2 にならない塊の値: %d 個、最大 %d。r の範囲 %d〜%d、q の範囲 %d〜%d' % (len(exo), CX, min(rq(c)[0] for c in exo), max(rq(c)[0] for c in exo), min(rq(c)[1] for c in exo), max(rq(c)[1] for c in exo)))
# 解析的な確認: c > 141 なら 3 < i3/i2 < 4、c > 47 なら 2 < i1/i4 < 3（本文）。数値でも 2^24 まで全部と、その上は乱数で見る
import random
bad = sum(1 for c in range(CX + 1, 1 << 20) if rq(c) != (3, 2)); rnd = random.Random(1)
bad += sum(1 for _ in range(2000000) if rq(rnd.randrange(1 << 20, 1 << 32)) != (3, 2))
print('  %d より大きい塊で 3・2 にならないもの（2^20 まで全部＋乱数 200 万個）: %d' % (CX, bad))
RX = max(rq(c)[0] for c in exo); QN = min(rq(c)[1] for c in exo); QX = max(rq(c)[1] for c in exo)
EPS = Fr(1, 1000)                     # 切り捨ての和（n·2^-13 未満）と丸めの分の余裕

def feasible(n, e, jmax=None):
    """塊 n 個、うち小さい塊（c <= CX）が e 個のとき、S と b0b1 を両立できる (w, Y の範囲) を返す"""
    out = []
    smax = (n - e) * 0xFFFFFFFF + e * CX; smin = (n - e) * (CX + 1)
    Rlo, Rhi = 3 * (n - e), 3 * (n - e) + RX * e
    Qlo, Qhi = 2 * (n - e) + QN * e, 2 * (n - e) + QX * e
    Jhi = 4 * (2 ** n - 1) if jmax is None else jmax
    Ylo, Yhi = Fr(Rlo) - n - EPS, Fr(Rhi + Jhi) + EPS          # Y = R + ΣJ − F1 − E1
    Zlo, Zhi = Fr(Qlo) - n - EPS, Fr(Qhi) + EPS                # Z = Q − F3 − E3
    for w in range(0, n + 1):
        Sw = S + (w << 32)
        # σ = (S_w − Y)/K_PI が [smin, smax] に入る Y だけ
        ylo = max(Ylo, Sw - KP * smax); yhi = min(Yhi, Sw - KP * smin)
        if ylo > yhi: continue
        # T16 ≡ ρ·(S_w − Y) + Z。値の範囲 [a, b] が、mod 65536 で [Tlo, Thi] と重なるか
        a = RHO * (Sw - yhi) + Zlo; b = RHO * (Sw - ylo) + Zhi
        m0 = math.floor((a - Thi) / M16)
        for m in range(m0, m0 + 3):
            lo, hi = Tlo + m * M16, Thi + m * M16
            if hi >= a and lo <= b:
                # 重なる部分に対応する Y の範囲（Z の幅も含めた外側の見積もり）
                y1 = Sw - (min(b, hi) - Zlo) / RHO; y2 = Sw - (max(a, lo) - Zhi) / RHO
                out.append((w, float(max(ylo, y1)), float(min(yhi, y2))))
    return out

print()
print('塊の数 n、小さい塊の数 e ごとに、S と b0b1 を両立できるか（Y = R + ΣJ − Σf1 の範囲。なければ「不可能」）')
for n in range(1, 8):
    for e in range(0, n + 1):
        r = feasible(n, e)
        Yhi = 3 * (n - e) + RX * e + 4 * (2 ** n - 1)
        print(' n=%d e=%d: Y の上限 %4d  %s' % (n, e, Yhi, '不可能' if not r else '  '.join('w=%d: Y %.1f〜%.1f' % x for x in r)))
