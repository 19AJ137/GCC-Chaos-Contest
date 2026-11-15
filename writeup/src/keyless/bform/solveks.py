"""solveks.py — 平文 M と連鎖フラグが分かっている範囲で、A・B の式から鍵ストリームの語 K_n（外側の鍵 k）・K'_n（内側の鍵 v）を解く（鍵なし。2026-11-12）
  A = E_k(E_v(M))。B は指定の形（例 Dv,Dk = D_k(D_v(M))、Ek,Ev = E_v(E_k(M))）。GF(2) の連立 1 次方程式で、一意に決まる語を挙げる。
  python3 -I solveks.py head132_fixed.json 形 [N]     （例: 形 = Ek,Ev。出力は ks_Ek_Ev.pkl）"""
import sys, os, pickle
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
from bform import apply, W, name
def solve(form, N, flags, Aw, Bw, Mw):
    Asym = apply(('Ev', 'Ek'), N, flags); Bsym = apply(form, N, flags)
    idx = {}; eqs = []
    def bit(s):
        if s not in idx: idx[s] = len(idx)
        return 1 << idx[s]
    for S, V in ((Asym, Aw), (Bsym, Bw)):
        for n in range(N):
            mask = 0; val = V[n]
            for s in S[n]:
                if s[0] == 'M': val ^= Mw[s[1]]
                else: mask ^= bit(s)
            eqs.append((mask, val))
    # 掃き出し（RREF）
    piv = {}   # 最下位ビット -> (mask, val)
    incons = 0
    for mask, val in eqs:
        for p, (m2, v2) in piv.items():
            if mask >> p & 1: mask ^= m2; val ^= v2
        if mask == 0:
            incons += (val != 0); continue
        p = (mask & -mask).bit_length() - 1
        for q in list(piv):
            m2, v2 = piv[q]
            if m2 >> p & 1: piv[q] = (m2 ^ mask, v2 ^ val)
        piv[p] = (mask, val)
    inv = {v: k for k, v in idx.items()}; det = {}
    for p, (m, v) in piv.items():
        if m == 1 << p: det[inv[p]] = v
    return det, len(idx), len(piv), incons
def main():
    from headdata import load_head
    d = load_head(sys.argv[1]); form = tuple(sys.argv[2].split(',')); N = int(sys.argv[3]) if len(sys.argv) > 3 else 132
    A = open(os.path.join(CT, 'A.bin'), 'rb').read(); B = open(os.path.join(CT, 'B.bin'), 'rb').read()
    Aw = [int.from_bytes(W(A, n), 'big') for n in range(N)]; Bw = [int.from_bytes(W(B, n), 'big') for n in range(N)]; Mw = [int.from_bytes(W(d['MA'], n), 'big') for n in range(N)]
    det, nunk, rank, incons = solve(form, N, d['flags'], Aw, Bw, Mw)
    print('B = %s: 未知数 %d 語、階数 %d、矛盾する式 %d、一意に決まった語 %d 個' % (name(form), nunk, rank, incons, len(det)))
    for key, lab in (('K', '外側の鍵 k の鍵ストリーム K_n'), ('V', '内側の鍵 v の鍵ストリーム K\'_n')):
        ns = sorted(s[1] for s in det if s[0] == key)
        print('  %s: %d 語。語の番号: %s' % (lab, len(ns), ns))
    for s in det:
        if s[0].startswith('IV'): print('  %s = %08x' % (s[0], det[s]))
    pickle.dump(dict(det=det, form=form, N=N), open(os.path.join(D, 'ks_%s.pkl' % '_'.join(form)), 'wb'))
    return det
if __name__ == '__main__': main()
