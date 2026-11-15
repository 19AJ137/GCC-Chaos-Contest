# コード

GCC6 の再実装と、懸賞の暗号文の解読・パスワード復元に使ったコードです。説明は [`../docs/`](../docs/) にあります。

```
src/
├── verify.py        復元した内部状態だけで、懸賞の暗号文 3 本を復号して確かめる（約 1 秒）
├── check_key.py     鍵（パスワード）の候補を、内部状態と暗号文 3 本で確かめる
├── TARGETS.json     2 本の鍵 k・v の「鍵スケジュールの出力」（check_key.py が読む目標値）
├── gcc6/            GCC6 の参照実装と、実物の DLL を動かすエミュレーター
├── keyless/         鍵なしの解読（内部状態の復元と平文の復号）       … docs/writeup.md
└── keyrec/          内部状態からの鍵（パスワード）の復元             … docs/password.md
```

暗号文は [`../ciphertext/`](../ciphertext/)、平文は [`../plaintext/`](../plaintext/) を読みます（スクリプトの場所からの相対パスなので、どこから実行してもよい）。

## 必要なもの

| もの | 使うところ |
|---|---|
| Python 3.8 以上（3.10 と 3.13 で確認） | 全部 |
| `numpy` | `keyless/`（`bform/`・`init/`） |
| `unicorn`・`pefile` | `gcc6/` の実物の DLL のエミュレーター（`emu.py`・`emu16f.py`）、`check_key.py --emu16` |
| `lhafile` | `gcc6/extract.py`（配布物から DLL を取り出す） |
| `wordfreq` | `keyrec/` の辞書を使う探索（`vdict.py` など） |
| gcc（OpenMP） | `keyless/`・`keyrec/` の C のプログラム（`make`） |

```
pip install numpy unicorn pefile lhafile wordfreq
```

- 例では `python3 -I -B`（隔離モード）で実行している。`pip install --user` で入れたパッケージは隔離モードでは読めないので、その場合は `-I` を外す。
- C のプログラムは、浮動小数点の計算順を DLL と同じにするため、FMA の縮約と fast-math を禁止してコンパイルする（`-ffp-contract=off -fno-fast-math`。`Makefile` に入れてある）。

## すぐ確かめる

`writeup/` フォルダで実行する。

```sh
# 鍵を使わずに、復元した内部状態だけで 3 本を復号し、B を作り直す
python3 -I -B src/verify.py

# 鍵 k と、鍵 v と等価な 21 バイトの鍵で 3 本を復号する（「総合: 当たり」になる）
python3 -B src/check_key.py --khex 8ea590ce814594ad8cf5814549495349 --vhex d9e88b7e8be072ab26c423dbd356b3dde8c18df96a
python3 -B src/check_key.py --k "磁石・発光・IISI" --v "ﾙ闍~金rｫ&ﾄ#ﾛﾓVｳﾝ霖笹j"      # 文字列は cp932 にして渡す

# 鍵 v の候補を 1 本だけ確かめる（当たりなら、鍵ストリームの語 1〜8 が一致する）
python3 -B src/check_key.py --v "候補の文字列"

# check_key.py の自己検査（既知の鍵 2 本で目標と暗号文を作り直す）
python3 -B src/check_key.py --selftest
```

## gcc6/ — 参照実装

| ファイル | 内容 |
|---|---|
| `gcc.py` | **参照実装**。`key_schedule(key)`・`str_cryption(key, data, encrypt)`・`file_cryption(key, data, encrypt)`（`encrypt=False` で復号）。`python3 gcc.py 鍵` で鍵スケジュールの出力と例を表示 |
| `gccfast.h`・`gcc_vec.c` | C の移植（`gcc -O2 -ffp-contract=off -fno-fast-math gcc_vec.c -lm`。標準入力の 16 進の鍵ごとに、写像の番号・空回し・S・T と、0 を復号した 16 語を出す） |
| `extract.py`・`szdd.py` | 配布物（`../../software/` の自己解凍 exe）から DLL と `CONTEST.TXT` を取り出す |
| `emu.py`・`test_vs_emu.py`・`test_file.py` | 32bit 版の DLL（`Gcc6_32.dll`）を Unicorn で動かし、`gcc.py` と比べる（`Str_cryption` 600 件、`File_cryption` 60 件） |
| `emu16.py`・`emu16f.py`・`ne.py` | 16bit 版の DLL（`gcc6.dll`）を Unicorn の実モードで動かす。精度 64 ビット（制御語 0x37F）と 53 ビット（0x27F） |

実物の DLL で確かめる手順:

```sh
python3 -B src/gcc6/extract.py          # software/ から src/gcc6/dll/ に Gcc6_32.dll・gcc6.dll・contest_*.txt を取り出す
python3 -B src/gcc6/test_vs_emu.py      # tests 600 mismatches 0
python3 -B src/gcc6/test_file.py        # File tests mismatches: 0
python3 -B src/gcc6/emu16f.py           # 16bit 版の File_cryption と gcc.py の比較（cases 88 mismatch 0）
```

配布物の再配布の条件は「自己解凍前の状態で」なので、取り出したファイル（`gcc6/dll/`）はリポジトリに入れない（`.gitignore` 済み）。

## keyless/ — 鍵なしの解読（[docs/writeup.md](../docs/writeup.md)）

ライトアップの節の順に並んでいる。時間は 2 コアの計算機での目安。

| 段 | 場所 | 実行 |
|---|---|---|
| §3〜§4 鍵なしの観測と、内側の連鎖の検査 | `puzzles/` | `python3 -I -B puzzles/puzzles.py`、`nest.py`〜`nest5.py`、`nest_null.py`（帰無）、`mxdump.py`（2 本の平文の XOR を `mx_*.bin` に書く）。どれも 1 秒前後。連鎖フラグの読みは `real_flags_*.txt` |
| §5 平文とフラグの同時探索（言語モデル） | `ttp/` | `python3 ttp/lm2.py`（言語モデルを作る）→ `python3 ttp/fulldec2.py -crib ttp/cribE1.txt -hi 132 -o fd2_head.pkl`（約 5 分）。**コーパスが要る**（下の注意）。結果は `bform/head132_lm.json` に入れてある |
| §6 B の形の判定 | `bform/` | `python3 -I -B bform/bform.py bform/head132_fixed.json`（数秒） |
| §7 鍵ストリームの連立と写像との照合 | `bform/` | `cd bform`、`python3 -I -B solveks.py head132_fixed.json Ek,Ev`、`python3 -I -B solveks.py head132_fixed.json Dv,Dk`、`python3 -I -B maptest.py ks_Ek_Ev.pkl,ks_Dv_Dk.pkl head132_fixed.json`（数十秒〜数分） |
| §8.1 鍵 v（写像 1）の状態 | `bform/` | `python3 -I -B solve_v.py ks_Ek_Ev.pkl head132_fixed.json` → `python3 -I -B v_family.py`（数秒） |
| §8.2 鍵 k（写像 2）の状態 | `bform/` | `make -C ..` でビルドし、`./ksolve2 k_words.txt 33`（近似。約 2 分）→ 出力の Q・A・\|u\| とバイト 6/7 を渡して `./kexact2 k_words.txt 33 0.840524930010221 1.808451713696830 1.265805319365 e0 3f f6 3f 06 40 f6 3f > kexact_out.txt`（厳密。十数分。同じ出力を `kexact_out.txt` として同梱）→ `python3 -I -B kstate.py kexact_out.txt` |
| §9 先頭 1,024 バイトの復号 | `bform/` | `python3 -I -B decrypt_all.py`、`python3 -I -B facts.py` |
| §10 初期値と E の 2 ブロック目 | `init/` | `cd ../init`、`vcand.py` → `eblock2.py` → `eblock2x.py` → `edec2.py` → `kinit.py` → `mkstate.py` の順（各数秒）。最後の `state.json` が `verify.py` に入れた値 |

- `keystreams.tsv` は、語ごとの鍵ストリーム K_n・K′_n と連鎖フラグ f_n・f′_n（1 ブロック目 256 語と E の 2 ブロック目）。
- `bform/head132_lm.json` は §5 の言語モデルが読んだ語 0〜131 の平文とフラグ（バイト 335・345 に読み違いがある）、`head132_fixed.json` はそれを直したもの（§8.3）。
- §5 の言語モデルの材料（一般の英文・日本語の文書と、IISI の公開ページ・配布物の説明書）は、著作物なので入れていない。やり直すときは、`ttp/lm.py`・`ttp/lm2.py` の説明に従って用意する（一般の英文は環境変数 `LM_CORPUS`、IISI のページは `IISI_PAGES`）。§6 以降は `head132_fixed.json` から始められる。
- `init/state0.py` は途中段階の値で、鍵 v の IV が仮の値（`b4 11 b9 57`）になっている。正しい値は `81 34 9d 3a`（`vcand.py` が決める）。先頭 1,024 バイトの復号結果は同じ。
- スクリプトは途中の結果を `*.pkl` としてスクリプトと同じフォルダに書く（`.gitignore` 済み）。

## keyrec/ — 鍵（パスワード）の復元（[docs/password.md](../docs/password.md)）

ファイルごとの説明は [docs/password.md](../docs/password.md) の §11。よく使うもの:

```sh
cd src/keyrec
make                                   # C のプログラムをビルド
./spjoin2 selftest                     # 突き合わせの自己検査
python3 -B xcheck_gcc.py 300           # 再実装 gks.py と参照実装 gcc.py の突き合わせ
python3 -B kall.py control 4 sjis 2    # 対照（乱数の全角 16 バイトの鍵を、内部状態だけから復元する）
python3 -B kall.py real 4              # 鍵 k の 16 バイトの探索（2 コアで約 11 分。磁石・発光・IISI の並べ替えが出る）
python3 -B check_short_v.py            # 21〜28 バイトの等価な鍵 38 本を実物の 3 本で確かめる
python3 -B short/vshort_bound.py       # 鍵 v は 20 バイト以下ではありえないことの確かめ
```

鍵 v の仮説を当てるときは、`vpart.py` の `build_hyp`・`run_vres`（「既知のバイト列＋未知の区間」の形。§8.4）か、`check_key.py --v`（文字列そのもの）を使う。
