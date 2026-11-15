"""state0.py — 途中段階の値（1 ブロック目の「語 0 の時点の写像の状態」から鍵ストリームを作る版）。init/ のスクリプトが読み込む。
  最終版は ../../verify.py（鍵スケジュールの出力から始め、E の 2 ブロック目まで復号する）。
  注意: この版の V_Y0（鍵 v の語 0 の状態）と V_IV は、符号を決める前の仮の組（y_0 が正、IV_v = b411b957）。語 0 の鍵ストリームと IV の XOR は正しい値と同じなので、
        先頭 1,024 バイトの復号結果は最終版と同じになる。正しい組は y_0 が負、K'_0 = 8c8c921c、IV_v = 81349d3a（init/vcand.py が決める）。
  以下は 2026-11-12 の最初の版の verify.py の中身。"""
import sys, os, struct, hashlib
K_LO = 0x7ed8; K_Q = 0x3feae59489bb7ed8; K_A = 0x3ffcef6b106b7ed8; K_C = 0x3fc469add9107ed8; K_U0 = 0xbffbe15ee8017ed8; K_V0 = 0x3ff6c7cf07a97ed8; K_IV = 0x4d0478a6
V_LO = 0xab6e; V_P = 0x3ffecb5b5923ab6e; V_Y0 = 0x3f9c6e6daa6fab6e; V_IV = 0xb411b957
def u2d(u): return struct.unpack('<d', struct.pack('<Q', u))[0]
def d2u(x): return struct.unpack('<Q', struct.pack('<d', x))[0]
def setlow(x, lo): return u2d((d2u(x) & ~0xFFFF) | lo)
def word(x):
    b = struct.pack('<d', x); return (b[4] << 24) | (b[2] << 16) | (b[3] << 8) | b[5], (0 if b[6] > 0xE6 else 1)
def stream_k(n=256):
    Q, A, C, u, v = u2d(K_Q), u2d(K_A), u2d(K_C), u2d(K_U0), u2d(K_V0); ks = []; fl = []
    for i in range(n):
        if i:
            a = Q * u; b = a * u; nu = A * u - a * v; nv = C * v + b; u = setlow(nu, K_LO); v = setlow(nv, K_LO)
        w, f = word(v); ks.append(w); fl.append(f)
    return ks, fl
def stream_v(n=256):
    p, y = u2d(V_P), u2d(V_Y0); ks = []; fl = []
    for i in range(n):
        if i: y = setlow(y * y - p, V_LO)
        w, f = word((1.0 + y) * 0.5); ks.append(w); fl.append(f)
    return ks, fl
def enc(P, ks, fl, iv):
    out = []; prev = None
    for n, x in enumerate(P):
        c = x ^ ks[n]
        if fl[n]: c ^= (iv if prev is None else out[prev]); prev = n
        out.append(c)
    return out
def dec(Cw, ks, fl, iv):
    out = []; prev = None
    for n, c in enumerate(Cw):
        x = c ^ ks[n]
        if fl[n]: x ^= (iv if prev is None else Cw[prev]); prev = n
        out.append(x)
    return out
def words(d): return [int.from_bytes(d[4 * i:4 * i + 4], 'big') for i in range(len(d) // 4)]
def tobytes(ws): return b''.join(w.to_bytes(4, 'big') for w in ws)
def main():
    here = os.path.dirname(os.path.abspath(__file__)); cdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, '..', '..', '..', 'ciphertext')
    A = open(os.path.join(cdir, 'A.bin'), 'rb').read()[:1024]; B = open(os.path.join(cdir, 'B.bin'), 'rb').read()[:1024]; E = open(os.path.join(cdir, 'E.bin'), 'rb').read()[:1024]
    K, F = stream_k(); KP, FP = stream_v()
    MA = dec(dec(words(A), K, F, K_IV), KP, FP, V_IV)           # A: 外側 k を外し、内側 v を外す
    ME = dec(dec(words(E), K, F, K_IV), KP, FP, V_IV)
    B2 = enc(enc(MA, K, F, K_IV), KP, FP, V_IV)                  # B: 同じ平文を k で暗号化してから v で暗号化
    Bw = words(B); same = sum(1 for a, b in zip(B2, Bw) if a == b)
    MB = dec(dec(Bw, KP, FP, V_IV), K, F, K_IV)                  # B を直接復号（外側 v、内側 k）
    ma, me, mb = tobytes(MA), tobytes(ME), tobytes(MB)
    print('連鎖なしの語: 鍵 k %d 語、鍵 v %d 語、両方 %d 語' % (F.count(0), FP.count(0), sum(1 for i in range(256) if not F[i] and not FP[i])))
    print('A = B の語: %s' % [i for i in range(256) if A[4 * i:4 * i + 4] == B[4 * i:4 * i + 4]])
    print('両方とも連鎖なしの語: %s' % [i for i in range(256) if not F[i] and not FP[i]])
    print('A の平文から作り直した B と実物の B: %d / 256 語が一致（%s）' % (same, '全部一致' if same == 256 else '不一致あり'))
    print('B を直接復号した平文と A の平文: %s' % ('同じ' if mb == ma else '違う'))
    print('A・B の平文 SHA-256: %s' % hashlib.sha256(ma).hexdigest()); print('E の平文（先頭 1,024 B）SHA-256: %s' % hashlib.sha256(me).hexdigest())
    for nm, d, fn in (('A・B の平文（1,024 バイト、Shift_JIS）', ma, 'plain_A.txt'), ('E の平文（先頭 1,024 バイト、Shift_JIS）', me, 'plain_E_first1024.txt')):
        print('\n===== %s =====' % nm); print(d.decode('cp932', 'replace').replace('\r\n', '\n'))
        open(os.path.join(here, fn), 'wb').write(d)
if __name__ == '__main__': main()
