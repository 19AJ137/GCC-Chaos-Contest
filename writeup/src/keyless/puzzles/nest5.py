"""nest5.py — 行ごとに「どの前の語と XOR を取ると文章どうしの XOR らしくなるか」を、読みを決めずに探す（鍵なし、2026-11-12）
  X_n ⊕ (前の 10 語以内から 0〜3 語) で、最上位ビットが 4 列とも 0 になる組をすべて挙げ、ビット 5 が 0 のバイト数 → 0x20 未満のバイト数で並べる。
  読み A の組（nest4.py）が何番目に来るかを見る。
  python3 -I nest5.py [lo hi]"""
import sys, os, itertools, collections
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
from nest import W, X, msb, pr
from nest4 import XW, decode, same
def score(v): return (sum(1 for b in v if not b & 0x20), sum(1 for b in v if (b & 0x7f) < 0x20))
def readA_set(R):
    """読み A の MX を X の語の集合で表す"""
    out = {}; prevo = None; previ = None
    for r in R:
        n = r['n']; s = {n}
        if r['fo']:
            if prevo is None: prevo = n; continue
            s ^= {prevo}; prevo = n
        px = frozenset(s)
        if r['fi']:
            if previ is None: previ = px; continue
            s = set(px) ^ set(previ); previ = px
        out[n] = frozenset(s) - {n}
    return out
def main():
    lo, hi, oc, ic = (int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv) > 4 else (55, 118, 3, 2)
    R = decode(lo, hi, oc, ic); SA = readA_set(R); rk = collections.Counter(); byc = collections.defaultdict(list)
    for r in R:
        n = r['n']
        if n not in SA or n < lo + 6: continue
        cands = []
        prev = list(range(max(lo, n - 10), n))
        for k in range(0, 4):
            for T in itertools.combinations(prev, k):
                v = X(XW[n], *[XW[t] for t in T])
                if any(b & 0x80 for b in v): continue
                cands.append((score(v), -k, T, v))
        cands.sort(reverse=True)
        a = tuple(sorted(SA[n])); va = X(XW[n], *[XW[t] for t in a]); sa = score(va)
        rank = sum(1 for c in cands if c[0] > sa) + 1; ties = sum(1 for c in cands if c[0] == sa)
        cls = (r['fo'], r['fi']); byc[cls].append((rank, len(cands)))
        best = cands[0]
        print('  語 %3d 型(外%d内%d) X=%s: 読み A = ⊕%-16s %s %-4s 点 %s、順位 %d（同点 %d）/ %d。最良 ⊕%-16s %s %-4s 点 %s%s' % (
            n, r['fo'], r['fi'], msb(XW[n]), list(a), va.hex(), pr(va), sa, rank, ties, len(cands), list(best[2]), best[3].hex(), pr(best[3]), best[0], '   A=B' if n in same else ''))
    print()
    for cls in sorted(byc):
        v = byc[cls]; print('  型(外%d内%d): %d 行。読み A の組が 1 位（同点を含む）%d、3 位以内 %d、候補の数の平均 %.1f、順位の中央値 %d' % (cls[0], cls[1], len(v), sum(1 for r, c in v if r == 1), sum(1 for r, c in v if r <= 3), sum(c for r, c in v) / len(v), sorted(r for r, c in v)[len(v) // 2]))
if __name__ == '__main__': main()
