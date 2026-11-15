"""verify.py — 懸賞の暗号文 3 本を、復元した「鍵スケジュールの出力」だけで全部復号し、B を作り直して照合する（単独で動く確認用）

  構造:  A（A.bin、1,024 B）= E_k(E_v(M_A))      … 鍵 v で暗号化してから、鍵 k で暗号化
         E（E.bin、1,338 B）= E_k(E_v(M_E))      … 同じ 2 つの鍵、同じ順序
         B（B.bin、1,024 B）= E_v(E_k(M_A))      … A と同じ平文を、鍵の順序を逆にして暗号化
  どの層も GCC6 の File_cryption の通常の暗号化（モード 1）。

  鍵（パスワード）は使わない。下の K・V は、GCC6 の DLL が鍵から作る内部の値（鍵スケジュールの出力）を、鍵なしの解析で復元したもの:
         mode   … 写像の番号（k は 2、v は 1）
         k780, k788 … 写像の初期値（空回しの前）。写像 2 は u = k780、v = k788。写像 1 は y = 2·k780 − 1
         k7b0, k7b8 … 写像の定数のもと。写像 1 は p = 1.9 + 0.1·k7b0。写像 2 は P = 0.9 + 0.4·k7b0、Q = 1.72 − 0.8·P − 0.14·k7b8、A = Q·P + 1、C = 1 − Q
         S      … 連鎖の初期値（IV）で、空回しの回数（S mod 79、最小 11）も決める
         lo     … 状態と定数の下位 16 ビットに毎回書き込まれる値（b0・b1）
  このスクリプトは、DLL の処理（空回し、鍵ストリーム、連鎖、1,024 バイトごとの作り直し、端数バイト）を、この値から始めてそのままなぞる。
  鍵そのものが分かっていなくても、この値があれば、この 2 つの鍵で暗号化したファイルは長さによらず復号できる。

  値の幅について: 下の数は、条件を満たす値の 1 組である。上位のビットは決まっているが、次の範囲は暗号文からは決まらない（復号の結果は同じ）。
         lo: k は 0x7ec1〜0x7eec、v は 0xab0f〜0xabc5 のどれでもよい。k780・k788・k7b0・k7b8 は、およそ 10^-11 の幅。
         鍵 v の k780 は 0.80018317365…と 0.19981682634…（y の符号が逆）のどちらでもよい。鍵 v の k788・k7b8 は使われないので分からない。

  使い方:  python3 -I -B verify.py [-o 出力フォルダ] [暗号文のフォルダ]     （既定は ../ciphertext）
  出力: 2 本の平文（画面。-o を付けると A.txt・E.txt にも書く）、B の作り直しの一致数、SHA-256、../plaintext/ との照合
"""
import sys, os, struct, hashlib
K = dict(mode=2, k780=0x3fe500287933f358, k788=0x3fe4f462e3890aef, k7b0=0x3fc3ca0d1ca9ebf4, k7b8=0x3fe924aef33bd88d, S=0xa678044d, lo=0x7ed8)
V = dict(mode=1, k780=0x3fe99b19be354000, k7b0=0x3fcf8c8bdb300652, S=0x3a9d3481, lo=0xab6e)
def u2d(u): return struct.unpack('<d', struct.pack('<Q', u))[0]
def d2u(x): return struct.unpack('<Q', struct.pack('<d', x))[0]
def setlow(x, lo): return u2d((d2u(x) & ~0xFFFF) | lo)
# DLL の定数（.rdata のバイト列のまま。0.1 は真値と少し違う）
C_TENTH2 = u2d(0x3fb99999999999a0); C_1_9 = u2d(0x3ffe666666666666); C_0_4 = u2d(0x3fd999999999999a); C_0_9 = u2d(0x3feccccccccccccd)
C_M014 = u2d(0xbfc1eb851eb851ec); C_0_8 = u2d(0x3fe999999999999a); C_1_72 = u2d(0x3ffb851eb851eb85)
class Stream:
    """鍵スケジュールの出力から始める GCC6 の鍵ストリーム"""
    def __init__(self, st):
        self.st = st; self.lo = st['lo']; self.setup(u2d(st['k7b0']))
        self.warm = max(11, st['S'] % 79)
        for _ in range(self.warm): self.step()
    def setup(self, k7b0):
        st, lo = self.st, self.lo
        if st['mode'] == 1:
            self.y = setlow(u2d(st['k780']) * 2.0 - 1.0, lo); self.p = setlow(k7b0 * C_TENTH2 + C_1_9, lo)
        else:
            self.u = u2d(st['k780']); self.v = u2d(st['k788'])
            P = k7b0 * C_0_4 + C_0_9; Q = (u2d(st['k7b8']) * C_M014 - P * C_0_8) + C_1_72
            self.Q = setlow(Q, lo); self.A = setlow(Q * P + 1.0, lo); self.C = setlow(1.0 - Q, lo)
    def step(self):
        """写像を 1 歩進めて、出力の倍精度の値を返す"""
        if self.st['mode'] == 1:
            self.y = setlow(self.y * self.y - self.p, self.lo); return (1.0 + self.y) * 0.5
        a = self.Q * self.u; b = a * self.u; nu = self.A * self.u - a * self.v; nv = self.C * self.v + b
        self.u = setlow(nu, self.lo); self.v = setlow(nv, self.lo); return self.v
    def word(self):
        """鍵ストリームの 4 バイトと、連鎖なしかどうか"""
        b = struct.pack('<d', self.step()); return bytes((b[4], b[2], b[3], b[5])), b[6] > 0xE6
    def reseed(self, byte):
        """1,024 バイトごとの作り直し: k7b0 のバイト 5 を連鎖レジスタの 2 バイト目に置き換えて、定数と初期値を作り直す（空回しなし）"""
        b = bytearray(struct.pack('<Q', self.st['k7b0'])); b[5] = byte; self.setup(struct.unpack('<d', bytes(b))[0])
def file_cryption(st, data, encrypt, trace=None):
    """File_cryption と同じ処理（鍵の代わりに鍵スケジュールの出力 st を使う）"""
    g = Stream(st); chain = st['S'].to_bytes(4, 'little'); bl = chain[1]; out = bytearray()
    for off in range(0, len(data), 1024):
        blk = data[off:off + 1024]; nw = len(blk) // 4
        for i in range(nw):
            ks, nochain = g.word(); x = blk[4 * i:4 * i + 4]; o = bytes(a ^ b for a, b in zip(x, ks))
            if trace is not None: trace.append((ks, nochain))
            if not nochain:
                o = bytes(a ^ b for a, b in zip(o, chain)); chain = o if encrypt else x
            out += o
        if len(blk) % 4:
            ks, _ = g.word()
            for j in range(len(blk) % 4):
                d = blk[4 * nw + j]; c = ks[j] ^ d ^ bl; out.append(c); bl = c if encrypt else d
        g.reseed(chain[1])
    return bytes(out)
def main():
    here = os.path.dirname(os.path.abspath(__file__)); args = sys.argv[1:]; odir = None
    if args[:1] == ['-o']: odir = args[1]; args = args[2:]
    cdir = args[0] if args else os.path.join(here, '..', 'ciphertext'); pdir = os.path.join(here, '..', 'plaintext')
    A = open(os.path.join(cdir, 'A.bin'), 'rb').read(); B = open(os.path.join(cdir, 'B.bin'), 'rb').read(); E = open(os.path.join(cdir, 'E.bin'), 'rb').read()
    tk = []; tv = []
    MA = file_cryption(V, file_cryption(K, A, False, tk), False, tv)        # A: 外側 k を外し、内側 v を外す
    ME = file_cryption(V, file_cryption(K, E, False), False)                # E: 同じ（1,338 バイト。2 ブロック目と端数 2 バイトを含む）
    B2 = file_cryption(V, file_cryption(K, MA, True), True)                 # B: 同じ平文を k で暗号化してから v で暗号化
    MB = file_cryption(K, file_cryption(V, B, False), False)                # B を直接復号（外側 v、内側 k）
    A2 = file_cryption(K, file_cryption(V, MA, True), True); E2 = file_cryption(K, file_cryption(V, ME, True), True)
    same = sum(1 for i in range(0, 1024, 4) if B2[i:i + 4] == B[i:i + 4])
    fk = [0 if nc else 1 for _, nc in tk]; fv = [0 if nc else 1 for _, nc in tv]
    print('空回しの回数: 鍵 k %d 回、鍵 v %d 回（S mod 79）' % (max(11, K['S'] % 79), max(11, V['S'] % 79)))
    print('連鎖なしの語（先頭 1,024 バイト）: 鍵 k %d 語、鍵 v %d 語、両方 %d 語' % (fk.count(0), fv.count(0), sum(1 for i in range(256) if not fk[i] and not fv[i])))
    print('A = B の語: %s' % [i for i in range(256) if A[4 * i:4 * i + 4] == B[4 * i:4 * i + 4]])
    print('両方とも連鎖なしの語: %s' % [i for i in range(256) if not fk[i] and not fv[i]])
    print('A の平文から作り直した B と実物の B: %d / 256 語が一致（%s）' % (same, '全部一致' if B2 == B else '不一致あり'))
    print('B を直接復号した平文と A の平文: %s' % ('同じ' if MB == MA else '違う'))
    print('平文を暗号化し直すと元の暗号文に戻るか: A %s、E %s' % ('戻る' if A2 == A else '戻らない', '戻る' if E2 == E else '戻らない'))
    print('A・B の平文（%d バイト）SHA-256: %s' % (len(MA), hashlib.sha256(MA).hexdigest()))
    print('E の平文（%d バイト）SHA-256: %s' % (len(ME), hashlib.sha256(ME).hexdigest()))
    for nm, d, fn in (('A・B の平文（1,024 バイト、Shift_JIS）', MA, 'A.txt'), ('E の平文（1,338 バイト、Shift_JIS）', ME, 'E.txt')):
        print('\n===== %s =====' % nm); print(d.decode('cp932', 'replace').replace('\r\n', '\n'))
        ref = os.path.join(pdir, fn)
        if os.path.exists(ref): print('（../plaintext/%s と%s）' % (fn, '一致' if open(ref, 'rb').read() == d else '不一致'))
        if odir: os.makedirs(odir, exist_ok=True); open(os.path.join(odir, fn), 'wb').write(d)
if __name__ == '__main__': main()
