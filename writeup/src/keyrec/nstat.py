# 鍵 v の塊の数（長さ）の見立て: S と b0b1 の範囲から、Q = S - floor(K_PI * (塊の値の和)) の範囲が決まる。
# Q = (積の整数部 4 つの和 SJ) + 3 * (塊の数) - (f1 の和の整数部)。SJ は塊の数とともに増える（1 塊あたり平均 1.47 倍）ので、
# 塊の数ごとに「乱数の鍵で Q がこの範囲に入る割合」を出せば、どの長さがありそうかが分かる。
import sys, random, math
from gks import key_schedule, chunks_of, chunk_vals, K_PI, K_E
from sums import MP, SHP, ME, SHE, fl_mul
from stfilter import sigma_candidates, TV
def qwin(n, tail=0):
    c = sigma_candidates(TV['S'], TV['Tlo'], TV['Thi'], 0, 0, n, n, cmin=0x20000000)
    qs = sorted(set(((TV['S'] - fl_mul(MP, SHP, sg)) & 0xFFFFFFFF) for sg, j3, w in c))
    return qs
def gen(kind, L, rng):
    if kind == 'ascii': return bytes(rng.randrange(0x20, 0x7F) for _ in range(L))
    if kind == 'lower': return bytes(rng.randrange(0x61, 0x7B) for _ in range(L))
    if kind == 'alnum': return bytes(rng.choice(b'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789') for _ in range(L))
    if kind == 'kanji':            # 全角（第 1 水準の漢字）だけ
        out = b''
        while len(out) < L - 1: out += bytes([rng.randrange(0x88, 0x98), rng.choice([x for x in range(0x40, 0xFD) if x != 0x7F])])
        return out[:L] if len(out) >= L else out + b'I' * (L - len(out))
    if kind == 'klike':            # 鍵 k の形: 漢字 2 文字の語を「・」でつなぎ、最後に IISI
        out = b''
        while len(out) + 4 + 6 <= L: out += bytes([rng.randrange(0x88, 0x98), rng.choice([x for x in range(0x40, 0xFD) if x != 0x7F]), rng.randrange(0x88, 0x98), rng.choice([x for x in range(0x40, 0xFD) if x != 0x7F])]) + b'\x81\x45'
        pad = L - len(out) - 4
        return out + bytes(rng.randrange(0x30, 0x7B) for _ in range(pad)) + b'IISI'
    if kind == 'kata':             # カタカナ（0x83 0x40〜0x96）
        out = b''
        while len(out) < L - 1: out += bytes([0x83, rng.choice([x for x in range(0x40, 0x97) if x != 0x7F])])
        return out[:L] if len(out) >= L else out + b'I' * (L - len(out))
if __name__ == '__main__':
    rng = random.Random(20261114); NS = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
    print('鍵 v の Q の範囲（塊の数ごと。sigma の候補から）')
    W = {}
    for n in range(4, 13):
        qs = qwin(n); W[n] = (qs[0], qs[-1]) if qs else None
        print('  塊 %2d 個: Q = %s  （候補 %d 個）' % (n, '%d〜%d' % (qs[0], qs[-1]) if qs else 'なし', len(qs)))
    print('乱数の鍵で、Q が鍵 v の範囲に入る割合（%d 本ずつ）' % NS)
    print('  長さ  塊   ' + '  '.join('%-7s' % k for k in ['ascii', 'lower', 'alnum', 'kanji', 'klike', 'kata']))
    for L in range(16, 49):
        n = (L + 3) // 4; row = []
        for kind in ['ascii', 'lower', 'alnum', 'kanji', 'klike', 'kata']:
            hit = tot = 0
            for _ in range(NS):
                key = gen(kind, L, rng)
                if key is None or len(key) != L: continue
                r = key_schedule(key); sig = sum(chunks_of(key)); tot += 1
                # 鍵 v と同じ条件: sigma が、この鍵の S・b0b1 から見て「鍵 v の S・b0b1 の関係」と同じ位置にあるか
                Q = (r['S'] - fl_mul(MP, SHP, sig)) & 0xFFFFFFFF
                # 鍵 v の側: T16 の範囲から sigma の位置が決まり、Q の範囲が決まる。乱数の鍵では Q だけを見る（T16 は一様）
                w = W.get(n)
                if w and w[0] <= Q <= w[1]: hit += 1
            row.append('%5.1f%%' % (100.0 * hit / tot) if tot else '   -   ')
        print('  %3d  %3d   %s' % (L, n, '  '.join('%-7s' % x for x in row)))
