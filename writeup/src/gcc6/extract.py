"""extract.py — software/ の配布物（LHA の自己解凍 exe）から、GCC6 の DLL と懸賞の暗号文を取り出す。

  使い方:  pip install lhafile
           python3 -I extract.py [software フォルダ] [出力フォルダ]
             既定は ../../../software と ./dll（このスクリプトのあるフォルダの下）

  取り出すもの（SZDD 圧縮の *.??_ は展開する）:
    Gcc6_32.dll  32bit 版 GCC6（1997-01-23 ビルド）     net9729.exe の infgd/Gcc6_32.dl_    emu.py が使う
    gcc6.dll     16bit 版 GCC6（1996-09-02）             Net9718.exe の INFGD/GCC6.DL_       emu16.py・emu16f.py が使う
    contest_A.txt・contest_B.txt・contest_E.txt          3 本の配布物の CONTEST.TXT（CONTEST.TX_）。writeup/ciphertext/ の A.bin・B.bin・E.bin と同じもの

  注意: 配布物の再配布の条件は「自己解凍前の状態で」なので、取り出したファイルはリポジトリに入れない（dll/ は .gitignore 済み）。
"""
import sys, os, io, re, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from szdd import expand

# (配布物, 書庫の中の名前（大文字・小文字を区別しない）, 出力名, 展開後の SHA-256)
WANT = [
    ('net9729.exe',  'infgd/gcc6_32.dl_',  'Gcc6_32.dll',   '501159f0c2865ca423d66c9d1adedd82c8f964c0e850214ddd449a91af223763'),
    ('Net9718.exe',  'net9718/infgd/gcc6.dl_', 'gcc6.dll',  'c74624e38eba95fbc9216c3bf262bb843e8cba6661bdd2bd3c228326bd9487fb'),
    ('Net9718.exe',  'net9718/infgd/contest.txt', 'contest_A.txt', 'c0d74bd3bfd758a099bc10028866fee58e7e2911d1d2339f63c9ffa43be36d7b'),
    ('net9729.exe',  'infgd/contest.tx_',  'contest_B.txt', '7a707ef0c7399ad4d0e8c69633a72759c20001e49b0b8d0dda43c80842cc5fbf'),
    ('iisi9605.exe', 'iisi9605/infgd/contest.txt', 'contest_E.txt', '4891751a2e7484fd195e788376faeb516b8d096bcb4a8276aa23a7f176d20042'),
]


def open_sfx(path):
    import lhafile
    raw = open(path, 'rb').read()
    m = re.search(rb'..-lh[0-7d]-', raw, re.S)          # 自己解凍部のあとに続く LHA の最初のヘッダ
    if not m: raise SystemExit('LHA の書庫が見つからない: %s' % path)
    return lhafile.Lhafile(io.BytesIO(raw[m.start():]))


def main():
    sw = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, '..', '..', '..', 'software')
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'dll')
    os.makedirs(out, exist_ok=True)
    cache = {}; ok = True
    for exe, name, dst, sha in WANT:
        if exe not in cache: cache[exe] = open_sfx(os.path.join(sw, exe))
        lf = cache[exe]
        ent = [i.filename for i in lf.infolist() if i.filename.replace('\\', '/').lower() == name]
        if not ent: print('見つからない: %s の %s' % (exe, name)); ok = False; continue
        d = lf.read(ent[0])
        if d[:4] == b'SZDD': d = expand(d)[0]
        h = hashlib.sha256(d).hexdigest(); good = (h == sha); ok &= good
        open(os.path.join(out, dst), 'wb').write(d)
        print('%-14s %-28s -> %-14s %6d B  SHA-256 %s…  %s' % (exe, ent[0], dst, len(d), h[:16], '一致' if good else '不一致'))
    print('全部一致' if ok else '一致しないものがある')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
