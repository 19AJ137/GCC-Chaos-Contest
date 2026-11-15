"""decrypt_all.py — 復元した 2 つの鍵の状態（鍵ストリームと連鎖フラグ）で A・B・E を復号し、B を作り直して照合する（2026-11-12）
  A = E_k(E_v(M_A))、E = E_k(E_v(M_E))、B = E_v(E_k(M_A))。
  連鎖の初期値 IV_k・IV_v（鍵から決まる 32 ビット）は状態からは出ないので、先頭 2 語の平文（"Is t"・"rue?"）から逆算する。
  python3 -I decrypt_all.py"""
import sys, os, pickle
D = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, D)
CT = os.path.join(D, '..', '..', '..', 'ciphertext')   # 懸賞の暗号文（writeup/ciphertext の A.bin・B.bin・E.bin）
def W(d, i): return int.from_bytes(d[4 * i:4 * i + 4], 'big')
def dec_layer(Cw, ks, fl, iv):
    """暗号化の逆: P_n = C_n ⊕ K_n ⊕ f_n·(直前の連鎖語の暗号文、なければ IV)"""
    out = []; prev = None
    for n, c in enumerate(Cw):
        p = c ^ ks[n]
        if fl[n]:
            p ^= (iv if prev is None else Cw[prev]); prev = n
        out.append(p)
    return out
def enc_layer(Pw, ks, fl, iv):
    out = []; prev = None
    for n, p in enumerate(Pw):
        c = p ^ ks[n]
        if fl[n]:
            c ^= (iv if prev is None else out[prev]); prev = n
        out.append(c)
    return out
def tobytes(ws): return b''.join(w.to_bytes(4, 'big') for w in ws)
def show(t):
    out = []; i = 0
    while i < len(t):
        c = t[i]
        if c == 0x0d: out.append('\\r'); i += 1
        elif c == 0x0a: out.append('\\n\n'); i += 1
        elif c == 0x09: out.append('\\t'); i += 1
        elif 0x20 <= c < 0x7f: out.append(chr(c)); i += 1
        elif (0x81 <= c <= 0x9f or 0xe0 <= c <= 0xfc) and i + 1 < len(t):
            try: out.append(t[i:i + 2].decode('cp932')); i += 2
            except Exception: out.append('<%02x>' % c); i += 1
        elif 0xa1 <= c <= 0xdf: out.append(bytes([c]).decode('cp932')); i += 1
        else: out.append('<%02x>' % c); i += 1
    return ''.join(out)
def main():
    A = open(os.path.join(CT, 'A.bin'), 'rb').read(); B = open(os.path.join(CT, 'B.bin'), 'rb').read(); E = open(os.path.join(CT, 'E.bin'), 'rb').read()[:1024]
    Aw = [W(A, n) for n in range(256)]; Bw = [W(B, n) for n in range(256)]; Ew = [W(E, n) for n in range(256)]
    kf = pickle.load(open(os.path.join(D, 'k_family.pkl'), 'rb'))['fam']; vf = pickle.load(open(os.path.join(D, 'v_family.pkl'), 'rb'))
    M0 = int.from_bytes(b'Is t', 'big'); M1 = int.from_bytes(b'rue?', 'big')
    best = None
    for ik, kv in enumerate(kf):
        for iv_, vv in enumerate(vf):
            K, F = kv['ks'], kv['flags']; KP, FP = vv['ks'], vv['flags']
            # IV の逆算（外側の最初の連鎖語と内側の最初の連鎖語）
            n1 = F.index(1); n1p = FP.index(1)
            # まず IV を 0 として外側を外し、内側の最初の連鎖語の式から IV_v、外側の最初の連鎖語の式から IV_k を解く（先頭 2 語の平文を使う）
            res = None
            for ivk_guess in (0,):
                C1 = dec_layer(Aw, K, F, 0)          # C'_n（語 n1 だけ IV_k の分ずれている）
                # 語 0・1 の平文から: 2 つの未知数 IV_k・IV_v を解く
                # M_n = C'_n ⊕ K'_n ⊕ f'_n·(C'_{p'(n)} or IV_v)。C'_{n1} = C1[n1] ⊕ IV_k
                def plain(ivk, ivv):
                    C = list(C1); C[n1] ^= ivk
                    return dec_layer(C, KP, FP, ivv)
                P00 = plain(0, 0)
                # 語 0・1 に IV がどう入るかを 1 ビットずつでなく語として見る（XOR なので線形）
                e_k = [a ^ b for a, b in zip(plain(0xFFFFFFFF, 0), P00)][:2]; e_v = [a ^ b for a, b in zip(plain(0, 0xFFFFFFFF), P00)][:2]
                # e_k[i]・e_v[i] は 0 か 0xFFFFFFFF（その語に IV が入るかどうか）
                need0 = P00[0] ^ M0; need1 = P00[1] ^ M1
                sol = None
                for ivk_in0, ivv_in0, ivk_in1, ivv_in1 in [(bool(e_k[0]), bool(e_v[0]), bool(e_k[1]), bool(e_v[1]))]:
                    # need0 = ivk·[k が語 0 に入る] ⊕ ivv·[v が語 0 に入る]、need1 も同様
                    import itertools
                    a, b, c, d = ivk_in0, ivv_in0, ivk_in1, ivv_in1
                    det_ = (a and d) != (b and c)
                    if a and not b and not c and d: sol = (need0, need1)
                    elif b and not a and c and not d: sol = (need1, need0)
                    elif b and not a and c and d: sol = (need1 ^ need0, need0)
                    elif a and not b and c and d: sol = (need0, need1 ^ need0)
                    elif a and b and c and not d: sol = (need1, need0 ^ need1)
                    elif a and b and d and not c: sol = (need0 ^ need1, need1)
                    elif b and not a and not c and not d: sol = (0, need0) if need1 == 0 else None
                    else: sol = None
                res = sol
            if res is None: continue
            ivk, ivv = res
            MA = plain(ivk, ivv)
            # B を作り直す: B = E_v(E_k(M_A))
            G = enc_layer(MA, K, F, ivk); B2 = enc_layer(G, KP, FP, ivv)
            nb = sum(1 for n in range(256) if B2[n] == Bw[n])
            if best is None or nb > best[0]: best = (nb, ik, iv_, ivk, ivv, MA, [n for n in range(256) if B2[n] != Bw[n]])
            print('k の並び %d × v の並び %d: IV_k = %08x、IV_v = %08x、作り直した B と実物の一致 %d / 256 語' % (ik, iv_, ivk, ivv, nb))
    nb, ik, iv_, ivk, ivv, MA, bad = best
    K, F = kf[ik]['ks'], kf[ik]['flags']; KP, FP = vf[iv_]['ks'], vf[iv_]['flags']
    print('\n採用: k の並び %d、v の並び %d。B の一致 %d / 256 語（合わない語 %s）' % (ik, iv_, nb, bad))
    C1 = dec_layer(Ew, K, F, ivk); ME = dec_layer(C1, KP, FP, ivv)
    ma = tobytes(MA); me = tobytes(ME)
    print('\n===== A・B の平文（1,024 バイト）=====\n' + show(ma)); print('\n===== E の平文（先頭 1,024 バイト）=====\n' + show(me))
    open(os.path.join(D, 'plain_A.bin'), 'wb').write(ma); open(os.path.join(D, 'plain_E_1024.bin'), 'wb').write(me)
    pickle.dump(dict(K=K, F=F, KP=KP, FP=FP, ivk=ivk, ivv=ivv, ik=ik, iv=iv_), open(os.path.join(D, 'keystreams.pkl'), 'wb'))
    print('\n外側の鍵 k: 連鎖なしの語 %d 個。内側の鍵 v: 連鎖なしの語 %d 個。両方とも連鎖なしの語: %s' % (F.count(0), FP.count(0), [n for n in range(256) if F[n] == 0 and FP[n] == 0]))
    print('A = B の語: %s' % [n for n in range(256) if Aw[n] == Bw[n]])
if __name__ == '__main__': main()
