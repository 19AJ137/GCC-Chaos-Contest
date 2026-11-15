"""check_key.py — 鍵（パスワード）の候補を、復元した内部状態（TARGETS.json）と懸賞の暗号文 3 本で確かめる。

  使い方（このフォルダで）:
    python3 -I -B check_key.py --k "鍵kの文字列" --v "鍵vの文字列"     文字列は cp932（Shift_JIS）のバイト列にして渡す
    python3 -I -B check_key.py --khex 6162... --vhex ...                16 進のバイト列
    python3 -I -B check_key.py --k "..."                                片方だけでもよい（その鍵だけ照合する）
    python3 -I -B check_key.py --k "..." --v "..." --emu16              A・E の復号を 16bit DLL のエミュレーターでも行う（unicorn が要る）
    python3 -I -B check_key.py --selftest [--emu16]                     既知の鍵 2 本で目標と暗号文を作り直し、手順が当たりを出すか確かめる

  判定:
    1. 鍵スケジュールの出力（gcc6/gcc.py。32bit DLL と全件一致を確認済みの再実装）を TARGETS.json と比べる:
       写像の番号、S（= IV）、T の下位 16 ビット（b0b1）が範囲内、k780・k788・k7b0・k7b8 が 1e-9 以内
       （鍵 v は k780 と k7b0 だけ。k780 は 1 − k780 でもよい）
    2. 最終判定: その鍵で作った鍵ストリームの語 1〜8 が復元したものと一致するか。語 0 も表示する
    3. 2 本そろっていれば: A・E を「k で復号 → v で復号」、B を「v で復号 → k で復号」して、../plaintext/A.txt・E.txt と比べ、
       平文を暗号化し直して A・B・E に戻るかも見る
    --emu16: A・E は 16bit 版で作られたので、16bit DLL（gcc6/dll/gcc6.dll。gcc6/extract.py で取り出す）を実モードで動かして
       精度 64 ビット（制御語 0x37F）と 53 ビット（0x27F）の両方で復号する（gcc6/emu16f.py。pip install unicorn）

  例: python3 -I -B check_key.py --khex 8ea590ce814594ad8cf5814549495349 --vhex d9e88b7e8be072ab26c423dbd356b3dde8c18df96a
"""
import sys, os, json, struct, argparse, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'gcc6'))
import gcc

TOL = 1e-9
u2d = lambda h: struct.unpack('<d', struct.pack('<Q', int(h, 16)))[0]


def load_targets(fn=os.path.join(HERE, 'TARGETS.json')):
    T = json.load(open(fn, encoding='utf-8'))
    out = {}
    for r in ('k', 'v'):
        t = T[r]
        d = dict(mode=t['mode'], S=int(t['S'], 16), lo=tuple(int(x, 16) for x in t['b0b1_range_pwk_h']),
                 ks=[int(x, 16) for x in t['keystream_words_0_8']], f={})
        for f in ('k780', 'k788', 'k7b0', 'k7b8'):
            if f + '_bits' in t: d['f'][f] = u2d(t[f + '_bits'])
        out[r] = d
    return out


def ks_words(key, n=9):
    st = gcc.key_schedule(key)
    return [int.from_bytes(gcc.ks_word(gcc.step(st))[0], 'big') for _ in range(n)]


def check_one(key, t, role):
    """1 本の鍵を、役割 role（'k' か 'v'）の目標と比べる。(最終判定, 項目ごとの結果) を返す"""
    try:
        st = gcc.key_schedule(key)
    except ValueError as e:
        return False, {'鍵スケジュール': '失敗 (%s)' % e}
    r = {}
    r['写像の番号'] = (st.mode == t['mode'], st.mode)
    r['S'] = (st.S == t['S'], '%08x' % st.S)
    lo = st.T & 0xffff
    r['b0b1'] = (t['lo'][0] <= lo <= t['lo'][1], '%04x' % lo)
    got = dict(k780=st.k780, k788=st.k788, k7b0=st.k7b0_saved, k7b8=st.k7b8)
    for f, x in t['f'].items():
        g = got[f]
        ok = abs(g - x) < TOL or (role == 'v' and f == 'k780' and abs(g - (1.0 - x)) < TOL)
        r[f] = (ok, repr(g))
    w = ks_words(key)
    r['鍵ストリーム語0'] = (w[0] == t['ks'][0], '%08x' % w[0])
    r['鍵ストリーム語1〜8'] = (w[1:9] == t['ks'][1:9], ' '.join('%08x' % x for x in w[1:9]))
    final = r['写像の番号'][0] and r['鍵ストリーム語1〜8'][0]
    return final, r


def show(role, key, final, r):
    print('== 鍵 %s の候補 %r（%d バイト）: %s' % (role, key, len(key), '当たり' if final else '外れ'))
    for k, v in r.items():
        if isinstance(v, tuple): print('   %-14s %s  %s' % (k, '一致' if v[0] else '不一致', v[1]))
        else: print('   %-14s %s' % (k, v))


CTDIR = os.path.join(HERE, '..', 'ciphertext'); PTDIR = os.path.join(HERE, '..', 'plaintext')


def read_contest(cdir=CTDIR, pdir=PTDIR):
    rd = lambda d, f: open(os.path.join(d, f), 'rb').read()
    return dict(A=rd(cdir, 'A.bin'), B=rd(cdir, 'B.bin'), E=rd(cdir, 'E.bin'),
                MA=rd(pdir, 'A.txt'), ME=rd(pdir, 'E.txt'))


def check_pair(kk, vk, C, emu16=False):
    """2 本そろったときの復号と、暗号化し直しの照合。全部合えば True"""
    fc = gcc.file_cryption
    MA = fc(vk, fc(kk, C['A'], False), False)
    ME = fc(vk, fc(kk, C['E'], False), False)
    MB = fc(kk, fc(vk, C['B'], False), False)
    res = {
        'A を k→v で復号 = 平文A': MA == C['MA'],
        'E を k→v で復号 = 平文E': ME == C['ME'],
        'B を v→k で復号 = 平文A': MB == C['MA'],
        '平文A を v→k で暗号化 = A': fc(kk, fc(vk, C['MA'], True), True) == C['A'],
        '平文E を v→k で暗号化 = E': fc(kk, fc(vk, C['ME'], True), True) == C['E'],
        '平文A を k→v で暗号化 = B': fc(vk, fc(kk, C['MA'], True), True) == C['B'],
    }
    if emu16:
        from emu16f import Emu16F, DLL16
        e = Emu16F(DLL16)
        for cw in (0x37F, 0x27F):
            for name, c, m in (('A', C['A'], C['MA']), ('E', C['E'], C['ME'])):
                x, _ = e.file_cryption(kk, c, mode=2, cw=cw)
                y, _ = e.file_cryption(vk, x, mode=2, cw=cw)
                res['16bit DLL（制御語 %#x）で %s を k→v で復号' % (cw, name)] = (y == m)
    for k, v in res.items(): print('   %-44s %s' % (k, '一致' if v else '不一致'))
    return all(res.values())


def parse_key(s, h):
    if h is not None: return bytes.fromhex(h)
    if s is not None: return s.encode('cp932')
    return None


def selftest(emu16=False):
    """既知の鍵 2 本（k 役 = 写像 2、v 役 = 写像 1）で目標と暗号文 3 本を作り直し、同じ手順で当たりになるかを見る"""
    kk = next(c for c in (b'SelfTestOuterKey%02d' % i for i in range(100)) if gcc.key_schedule(c).mode == 2)
    vk = next(c for c in (b'Self test inner key, rather long, no.%02d' % i for i in range(100)) if gcc.key_schedule(c).mode == 1)
    T = {}
    for role, key in (('k', kk), ('v', vk)):
        st = gcc.key_schedule(key); lo = st.T & 0xffff
        f = dict(k780=st.k780, k7b0=st.k7b0_saved) if role == 'v' else dict(k780=st.k780, k788=st.k788, k7b0=st.k7b0_saved, k7b8=st.k7b8)
        T[role] = dict(mode=st.mode, S=st.S, lo=(max(0, lo - 20), lo + 20), ks=ks_words(key), f=f)
    real = read_contest()
    fc = gcc.file_cryption
    C = dict(MA=real['MA'], ME=real['ME'])
    C['A'] = fc(kk, fc(vk, C['MA'], True), True); C['E'] = fc(kk, fc(vk, C['ME'], True), True); C['B'] = fc(vk, fc(kk, C['MA'], True), True)
    ok = True
    for role, key in (('k', kk), ('v', vk)):
        fin, r = check_one(key, T[role], role); show(role, key, fin, r); ok &= fin
    fin, r = check_one(kk[::-1], T['k'], 'k'); show('k（わざと違う鍵）', kk[::-1], fin, r); ok &= not fin
    ok &= check_pair(kk, vk, C, emu16=emu16)
    print('自己検査:', '合格' if ok else '不合格')
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--k'); ap.add_argument('--v'); ap.add_argument('--khex'); ap.add_argument('--vhex')
    ap.add_argument('--emu16', action='store_true'); ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--ciphertext', default=CTDIR); ap.add_argument('--plaintext', default=PTDIR)
    a = ap.parse_args()
    if a.selftest: sys.exit(0 if selftest(a.emu16) else 1)
    kk, vk = parse_key(a.k, a.khex), parse_key(a.v, a.vhex)
    if kk is None and vk is None: ap.error('--k / --khex か --v / --vhex のどちらかが要る')
    T = load_targets(); hit = {}
    for role, key in (('k', kk), ('v', vk)):
        if key is None: continue
        fin, r = check_one(key, T[role], role); show(role, key, fin, r); hit[role] = fin
        other = 'v' if role == 'k' else 'k'
        fin2, _ = check_one(key, T[other], other)
        if fin2: print('   ※ この鍵は鍵 %s の目標に当たる（役割が逆）' % other)
    if kk is not None and vk is not None:
        print('== 2 本での復号の照合')
        ok = check_pair(kk, vk, read_contest(a.ciphertext, a.plaintext), emu16=a.emu16)
        print('総合:', '当たり（平文 3 本が一致）' if ok else '外れ')
        sys.exit(0 if ok else 1)
    sys.exit(0 if all(hit.values()) else 1)


if __name__ == '__main__':
    main()
