# 「鍵 v が、この塊を（塊の境目に合った位置で）含む」という仮定が、S と b0b1 だけで否定できるかを見る。
# 残りの未知の塊（4 バイトそろったもの）nU 個の和 sigma_U の候補があるかどうか。候補が 0 なら、その仮定は成り立たない。
import sys
from gks import chunk_vals, chunks_of, map_of_i2, K_PI, K_E, h2d
from sums import MP, SHP, ME, SHE, fl_mul
TV = dict(S=0x3a9d3481, Tlo=0xab0f, Thi=0xabc5, m=1)

def known_sums(chunks):
    """既知の塊の S・T への寄与"""
    s = t = 0
    for c in chunks:
        (f, i, ss, tt) = chunk_vals(c); s += ss; t += tt
    return s & 0xFFFFFFFF, t & 0xFFFFFFFF

def sigma_candidates(S, Tlo, Thi, SK, TK, nU, ntot, cmin=0x20000000, sjmax=None):
    """未知の塊 nU 個の和の候補（sigma, j3, wrap）。S = SK + sum_U(i1 + 3) + SJ、T = TK + sum_U(i3 + 2)"""
    sjmax = 4 * (2**ntot - 1) if sjmax is None else sjmax
    Sp = (S - SK) & 0xFFFFFFFF; out = []
    for w in range(nU + 1):
        Sw = Sp + (w << 32)
        lo = int((Sw - 3 * nU - sjmax) / K_PI) - 3; hi = int((Sw - 3 * nU - 4 + nU) / K_PI) + 3
        lo = max(lo, nU * cmin); hi = min(hi, nU * 0xFFFFFFFF)
        if hi < lo: continue
        vlo = fl_mul(ME, SHE, lo); vhi = fl_mul(ME, SHE, hi)
        for j3 in range(nU):
            off = (2 * nU + TK - j3) & 0xFFFF
            x = (vlo + off) & 0xFFFF
            v = vlo + (max(0, Tlo - x) if x <= Thi else (65536 - x + Tlo))
            while v <= vhi:
                t16 = (v + off) & 0xFFFF
                if not (Tlo <= t16 <= Thi): v += 65536 - (Thi - Tlo + 1); continue
                s0 = -((-(v << SHE)) // ME)                # 最小の sigma で floor(K_E sigma) = v
                sg = s0
                while fl_mul(ME, SHE, sg) == v:
                    if lo <= sg <= hi: out.append((sg, j3, w))
                    sg += 1
                v += 1
    return out

if __name__ == '__main__':
    k = '磁石・発光・IISI'.encode('cp932')
    tests = [('鍵 k の 4 塊（磁石／・発／光・／IISI）', chunks_of(k)),
             ('鍵 k の先頭 3 塊（磁石／・発／光・）', chunks_of(k)[:3]),
             ('GoodLuckToIISIGoodLuckTo（6 塊）', chunks_of(b'GoodLuckToIISIGoodLuckTo')),
             ("Gao's Chaos Cryptosystem（6 塊）", chunks_of(b"Gao's Chaos Cryptosystem")),
             ('equivalent KEY: GCCv（等価な鍵の先頭 5 塊。陽性対照）', chunks_of(b'equivalent KEY: GCCv'))]
    for name, ch in tests:
        SK, TK = known_sums(ch)
        row = []
        for nU in range(1, 8):
            ntot = len(ch) + nU
            c = sigma_candidates(TV['S'], TV['Tlo'], TV['Thi'], SK, TK, nU, ntot)
            row.append('%d:%d' % (nU, len(set(x[0] for x in c))))
        print('%-52s 既知 %d 塊 | 未知の塊の数:sigma の候補  %s' % (name, len(ch), '  '.join(row)))
