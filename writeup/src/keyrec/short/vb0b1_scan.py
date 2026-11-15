# b0b1 の範囲（ab0f〜abc5）が、k780 のセルの中の位置や、隣の区切りによって変わらないかを確かめる。
#   写像 1 は毎回の値の下位 16 ビットを b0b1 で置き換えるので、初期値 y0 = 2·k780 − 1 は 2^-37 刻みの区切りでしか効かない。
#   k780 を区切り 1 つ分（3.638e-12）ずつずらした 17 通り × b0b1 の全 65,536 通りで、鍵ストリーム 335 語が鍵 v のものと同じになる組を全部拾う。
#   python3 short/vb0b1_scan.py
import sys, os, struct
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
from gks import h2d, C01b
from vpart import full_v, rep5, _rep, _b0
pk = struct.Struct('<d'); uq = struct.Struct('<Q')
def setlow(x, b): return pk.unpack(uq.pack((uq.unpack(pk.pack(x))[0] & ~0xFFFF) | b))[0]
def word(x):
    z = pk.pack(x); return ((z[4] << 24) | (z[2] << 16) | (z[3] << 8) | z[5], z[6] > 0xE6)
REF = full_v(_rep, _b0, 0xab6a); R1 = REF[:256]; R2 = REF[256:]
def same(k780, k7b0, b, warm=23):
    y = setlow(2.0 * k780 - 1.0, b); p = setlow(k7b0 * C01b + 1.9, b)
    for i in range(warm): y = setlow(y * y - p, b)
    for i in range(256):
        y = setlow(y * y - p, b)
        if word((1.0 + y) / 2.0) != R1[i]: return False
    y = setlow(2.0 * k780 - 1.0, b); p = setlow(rep5(k7b0, 0xed) * C01b + 1.9, b)
    for i in range(79):
        y = setlow(y * y - p, b)
        if word((1.0 + y) / 2.0) != R2[i]: return False
    return True
Q = 3.638e-12
print('k780 を区切り q 個分ずらしたときに、335 語が同じになる b0b1')
allb = set()
for q in range(-8, 9):
    ok = [b for b in range(65536) if same(_rep + q * Q, _b0, b)]
    allb.update(ok)
    runs = []
    for b in ok:
        if runs and b == runs[-1][1] + 1: runs[-1][1] = b
        else: runs.append([b, b])
    print('  q = %+d: %5d 通り  %s' % (q, len(ok), '  '.join('%04x〜%04x' % (a, c) for a, c in runs) or '-'), flush=True)
print('全部を合わせた b0b1: %d 通り、%04x〜%04x' % (len(allb), min(allb), max(allb)))
# k7b0 の側: セルの中の 5 点と、セルの外側（定数の区切りの隣）で同じことを見る（b0b1 は ab00〜abff）
for name, d in (('セルの下端近く', -2.8e-11), ('中ほど', 0.0), ('セルの上端近く', 2.8e-11), ('セルの外（下）', -6e-11), ('セルの外（上）', 6e-11), ('定数の区切りの隣（下）', -1.455e-10), ('定数の区切りの隣（上）', 1.455e-10)):
    ok = [b for b in range(0xa800, 0xb000) if same(_rep, _b0 + d, b)]
    print('  k7b0 %+.3g（%s）: %d 通り %s' % (d, name, len(ok), ('%04x〜%04x' % (ok[0], ok[-1])) if ok else '-'))
