"""eblock2.py — E（E.bin）の 2 ブロック目（バイト 1024〜1337）を、鍵なしで復号する。2026-11-12
  File_cryption は 1,024 バイトごとに写像を作り直す: k7b0 のバイト 5 を「その層の連鎖レジスタの 2 バイト目」に置き換えて定数を作り直し、
  状態は鍵から作った初期値（空回しなし）に戻す。
  鍵 v（写像 1）: 初期値 y_init は vcand.py が出した候補（lo ごとに 2 通りの大きさ）。定数 p' は p から計算（上位 48 ビットに ±1 の幅）。
  鍵 k（写像 2）: 初期値は使わない。平文の予想（"Fea" の続き "tures of GCC ..."）と v の鍵ストリームから k の鍵ストリームの語を出し、
                  連続する語が写像 2 の式（定数は Q'・A'・C' の近似値）に合うかで v の候補と k の語の上位ビットを同時に決める。
                  そのあと、状態と定数を倍精度で厳密に決めて、残りの語を回す。
  python3 -I -B eblock2.py [-selftest]"""
import sys, os, pickle, struct, time, itertools, importlib.util
import numpy as np
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
spec = importlib.util.spec_from_file_location('state0', os.path.join(D, 'state0.py')); V = importlib.util.module_from_spec(spec); spec.loader.exec_module(V)
sys.path.insert(0, os.path.join(D, '..', '..', 'gcc6')); import gcc
u2d, d2u, setlow = V.u2d, V.d2u, V.setlow
EFULL = open(os.path.join(CT, 'E.bin'), 'rb').read()
K, F = V.stream_k(); KP, FP = V.stream_v()
Ew1 = V.words(EFULL[:1024]); Cp1 = V.dec(Ew1, K, F, V.K_IV); ME1 = V.dec(Cp1, KP, FP, V.V_IV)
lk = max(n for n in range(256) if F[n]); lv = max(n for n in range(256) if FP[n])
RK0 = Ew1[lk]; RV0 = Cp1[lv]                     # ブロック 1 の終わりの連鎖レジスタ（k 層は暗号文、v 層は内側の暗号文）
E2 = EFULL[1024:]; NW2 = len(E2) // 4; Ew2 = V.words(E2[:4 * NW2]); TAIL = E2[4 * NW2:]
PRED = (b'tures of GCC Chaos Encryption are as follows:\r\n1)  Suitable for the Age of Multimedia.\r\n2)  The Ultimate in Security\r\n3)  Super High Speed\r\n'
        b'4)  Suitability to Various Operating Systems\r\n5)  Strength against Noise\r\n6)  Suitability for Use on Networks\r\n7)  Versatility\r\n')
Mw = V.words(PRED[:len(PRED) // 4 * 4])
def repl5(x, byte):
    b = bytearray(struct.pack('<d', x)); old = b[5]; b[5] = byte; return struct.unpack('<d', bytes(b))[0], old
# ---- 鍵 k のブロック 2 の定数（近似値と、上位 48 ビットの候補）
def k_params2(byte):
    Q48, A48, C48 = V.K_Q >> 16, V.K_A >> 16, V.K_C >> 16
    Cm = u2d((C48 << 16) | 0x8000); Qm = 1.0 - Cm; Am = u2d((A48 << 16) | 0x8000); Pm = (Am - 1.0) / Qm
    k7b0 = (Pm - gcc.C_0_9) / gcc.C_0_4; k7b8 = (Qm + Pm * gcc.C_0_8 - gcc.C_1_72) / gcc.C_M014
    k7b0n, old = repl5(k7b0, byte)
    P2 = k7b0n * gcc.C_0_4 + gcc.C_0_9; Q2 = (k7b8 * gcc.C_M014 - P2 * gcc.C_0_8) + gcc.C_1_72; A2 = Q2 * P2 + 1.0; C2 = 1.0 - Q2
    return dict(k7b0=k7b0, k7b8=k7b8, old=old, k7b0n=k7b0n, P=P2, Q=Q2, A=A2, C=C2)
def v_params2(byte):
    p48 = V.V_P >> 16; pm = u2d((p48 << 16) | 0x8000); k7b0 = (pm - gcc.C_1_9) / gcc.C_TENTH2
    k7b0n, old = repl5(k7b0, byte); p2 = k7b0n * gcc.C_TENTH2 + gcc.C_1_9
    return dict(k7b0=k7b0, old=old, k7b0n=k7b0n, p=p2)
def vstream2(yinit, p2, lo, n):
    """ブロック 2 の v 層: 初期値 yinit から。返り値 (語, 連鎖フラグ) を n 語"""
    y = yinit; ks = []; fl = []
    for _ in range(n):
        y = setlow(y * y - p2, lo); w, f = V.word((1.0 + y) * 0.5); ks.append(w); fl.append(f)
    return ks, fl
def inner_ct(M, ks, fl, reg):
    """v 層の暗号化（平文 M → 内側の暗号文 C'）。reg = 連鎖レジスタの初期値"""
    out = []
    for m, k, f in zip(M, ks, fl):
        c = m ^ k
        if f: c ^= reg; reg = c
        out.append(c)
    return out, reg
# ---- 仮説（出力 v の上位 16 ビット = b7・b6）の表
HYP = [(0x40, 0x00 | n) for n in range(8)] + [(0x3F, (e << 4) | n) for e in range(0xF, 0x3, -1) for n in range(16)]
HB7 = np.array([h[0] for h in HYP], dtype=np.uint64); HB6 = np.array([h[1] for h in HYP], dtype=np.uint64); HF = (HB6 <= 0xE6).astype(np.uint64) | (HB7 == 0x40).astype(np.uint64)
HF = np.where(HB7 == 0x40, 1, (HB6 <= 0xE6).astype(np.int64)).astype(np.uint64)
def vfrom(Kw, b7, b6, lo):
    """鍵ストリームの語（配列）と b7・b6 から出力の倍精度の値。語のバイト並びは (b4, b2, b3, b5)"""
    Kw = Kw.astype(np.uint64); b4 = (Kw >> np.uint64(24)) & np.uint64(0xFF); b2 = (Kw >> np.uint64(16)) & np.uint64(0xFF); b3 = (Kw >> np.uint64(8)) & np.uint64(0xFF); b5 = Kw & np.uint64(0xFF)
    bits = (b7 << np.uint64(56)) | (b6 << np.uint64(48)) | (b5 << np.uint64(40)) | (b4 << np.uint64(32)) | (b3 << np.uint64(24)) | (b2 << np.uint64(16)) | np.uint64(lo)
    return bits.view(np.float64)
def m24_of_word(Kw):
    Kw = Kw.astype(np.int64); b4 = (Kw >> 24) & 0xFF; b3 = (Kw >> 8) & 0xFF; b5 = Kw & 0xFF; return (b5 << 16) | (b4 << 8) | b3
def m24_of_val(x):
    u = x.view(np.uint64); return ((u >> np.uint64(24)) & np.uint64(0xFFFFFF)).astype(np.int64), ((u >> np.uint64(48)) & np.uint64(0xFF)), (u >> np.uint64(56))
def consist(X, Ecw, R0, Q, A, C, lo, nchk=3):
    """X[n] = K_n xor f_n·(k 層のレジスタ)。Ecw[n] = 暗号文の語（連鎖したとき次のレジスタになる）。語 0・1 の仮説を全部ためし、語 2 以降を予想と照合。
       返り値: 通った仮説の一覧 [(h0, h1, [語 2 以降で決まった (b7,b6)…])]"""
    nh = len(HYP); i0 = np.repeat(np.arange(nh), nh); i1 = np.tile(np.arange(nh), nh)
    f0 = HF[i0]; K0 = np.uint64(X[0]) ^ (f0 * np.uint64(R0)); v0 = vfrom(K0, HB7[i0], HB6[i0], lo); r1 = np.where(f0 == 1, np.uint64(Ecw[0]), np.uint64(R0))
    f1 = HF[i1]; K1 = np.uint64(X[1]) ^ (f1 * r1); v1 = vfrom(K1, HB7[i1], HB6[i1], lo); r2 = np.where(f1 == 1, np.uint64(Ecw[1]), r1)
    usq = (v1 - C * v0) / Q; ok = usq >= 0
    i0, i1, v0, v1, usq, reg = i0[ok], i1[ok], v0[ok], v1[ok], usq[ok], r2[ok]
    hist = [np.stack([HB7[i0], HB6[i0]], 1), np.stack([HB7[i1], HB6[i1]], 1)]
    vprev, vcur = v0, v1
    for n in range(2, 2 + nchk):
        usq = usq * (A - Q * vprev) ** 2; vp = C * vcur + Q * usq
        m24p, b6p, b7p = m24_of_val(vp); fp = np.where(b7p == 0x40, 1, (b6p <= 0xE6).astype(np.int64)).astype(np.uint64)
        Kn = np.uint64(X[n]) ^ (fp * reg); m24i = m24_of_word(Kn); d = np.abs(m24p - m24i); keep = (d <= 2) | (d >= (1 << 24) - 2)
        if not keep.any(): return []
        usq, vprev, reg = usq[keep], vcur[keep], np.where(fp[keep] == 1, np.uint64(Ecw[n]), reg[keep]); hist = [h[keep] for h in hist]
        vcur = vfrom(Kn[keep], b7p[keep], b6p[keep], lo); hist.append(np.stack([b7p[keep], b6p[keep]], 1))
    return [tuple((int(h[j][0]), int(h[j][1])) for h in hist) for j in range(len(hist[0]))]
def main():
    t0 = time.time()
    print('E: 全 %d バイト。2 ブロック目 %d バイト（%d 語と端数 %d バイト）' % (len(EFULL), len(E2), NW2, len(TAIL)))
    print('ブロック 1 の終わりの連鎖レジスタ: k 層 = E の語 %d = %08x（2 バイト目 %02x）、v 層 = 内側の暗号文の語 %d = %08x（2 バイト目 %02x）' % (lk, RK0, (RK0 >> 16) & 0xFF, lv, RV0, (RV0 >> 16) & 0xFF))
    if '-selftest' in sys.argv:
        # ブロック 1 の真の値で、仮説の照合が正解 1 つだけを残すか
        Q, A, C = u2d(V.K_Q), u2d(V.K_A), u2d(V.K_C)
        for n0 in (40, 100, 137, 200):
            reg = V.K_IV
            for n in range(n0):
                if F[n]: reg = Ew1[n]
            X = [Ew1[n] ^ Cp1[n] for n in range(n0, n0 + 6)]
            res = consist(X, Ew1[n0:n0 + 6], reg, Q, A, C, V.K_LO, nchk=3)
            print('  自己検査 語 %d から: 通った仮説 %d 個 %s' % (n0, len(res), res[:3]))
        return
    kp = k_params2((RK0 >> 16) & 0xFF); vp = v_params2((RV0 >> 16) & 0xFF)
    print('鍵 k: k7b0 ≈ %.15f（バイト 5 = %02x → %02x）、k7b8 ≈ %.15f。ブロック 2 の定数の近似 Q\' = %.15f、A\' = %.15f、C\' = %.15f' % (kp['k7b0'], kp['old'], (RK0 >> 16) & 0xFF, kp['k7b8'], kp['Q'], kp['A'], kp['C']))
    print('鍵 v: k7b0 ≈ %.15f（バイト 5 = %02x → %02x）。ブロック 2 の定数の近似 p\' = %.15f' % (vp['k7b0'], vp['old'], (RV0 >> 16) & 0xFF, vp['p']))
    cands = pickle.load(open(os.path.join(D, 'vcand.pkl'), 'rb')); found = []; ntest = 0
    p48c = d2u(vp['p']) >> 16
    for lo, Sv, warm, kw, pts in cands:
        for idx in sorted(set(abs(int(i)) for i in pts)):
            yinit = u2d((idx << 16) | lo)
            for dp in (-1, 0, 1, -2, 2):
                p2 = u2d(((p48c + dp) << 16) | lo); ks, fl = vstream2(yinit, p2, lo, 8); cin, _ = inner_ct(Mw[:8], ks, fl, RV0)
                X = [Ew2[n] ^ cin[n] for n in range(8)]; ntest += 1
                res = consist(X, Ew2[:8], RK0, kp['Q'], kp['A'], kp['C'], V.K_LO, nchk=4)
                for r in res: found.append((lo, Sv, idx, dp, r)); print('  ★ 通った: lo_v = %04x、y_init の idx = %d（%.17g）、p\' の上位 48 ビット = 近似 %+d、k の語の (b7,b6) = %s' % (lo, idx, yinit, dp, r), flush=True)
    print('試した v の候補 %d 通り。通ったもの %d、%.0f 秒' % (ntest, len(found), time.time() - t0))
    pickle.dump(dict(found=found, kp=kp, vp=vp, p48c=p48c), open(os.path.join(D, 'eblock2_found.pkl'), 'wb'))
if __name__ == '__main__': main()
