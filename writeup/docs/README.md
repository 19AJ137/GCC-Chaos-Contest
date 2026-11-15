# 解読の記録

GCC カオス暗号の懸賞問題を、暗号文だけから解いた記録です。**平文と鍵のネタバレを含みます。**

## 結果

| | |
|---|---|
| 暗号の構造 | 懸賞の暗号文 3 本は、**2 本の鍵 k・v による二重暗号化**だった。A = E_k(E_v(M_A))、E = E_k(E_v(M_E))、B = E_v(E_k(M_A))（B は A と同じ平文で、鍵の順序が逆）。どの層も `File_cryption` の暗号化 |
| 平文 | 3 本とも全バイト（[`../plaintext/`](../plaintext/)）。A・B は「Is true? Do you think you can decrypte GCC? No, you cann't! Because this is GCC6 only!」で始まる英語と日本語の文章、E は「Hi, I am Gao@IISI.」で始まる英語の挨拶と GCC の説明 |
| 内部状態 | 2 本の鍵の「鍵スケジュールの出力」（DLL が鍵から作る値の一式）を、鍵なしで全部復元した（[`../src/TARGETS.json`](../src/TARGETS.json)、[`../src/verify.py`](../src/verify.py)） |
| 鍵 k（外側） | **`磁石・発光・IISI`**（Shift_JIS で 16 バイト。[`../password/k.txt`](../password/k.txt)） |
| 鍵 v（内側） | **元の文字列は未特定。** 同じ鍵ストリームを出す 21 バイトの鍵 `ﾙ闍~金rｫ&ﾄ#ﾛﾓVｳﾝ霖笹j` を作った（[`../password/v_equivalent.txt`](../password/v_equivalent.txt)）。元の鍵は 21 バイト以上（33〜40 バイトが有力） |

鍵 k と、鍵 v と等価な鍵があれば、懸賞の 3 本はインフォガード（[`../../software/`](../../software/)）と同じ手順、つまり「復号化」を 2 回行うだけで読めます（実物の DLL で確認）。A・E は「鍵 k で復号化 → 鍵 v（等価な鍵）で復号化」、B は「鍵 v（等価な鍵）で復号化 → 鍵 k で復号化」です。

## 読む順

1. [algorithm.md](algorithm.md) — GCC6 の仕様（DLL から復元したアルゴリズム、実装の注意、版）
2. [writeup.md](writeup.md) — 解読のライトアップ（鍵なしで内部状態を復元して、平文を読むまで）
3. [password.md](password.md) — 内部状態からのパスワード復元（鍵 k の特定、鍵 v の等価な鍵、鍵 v についての制約と調べた範囲）

コードの使い方は [`../src/README.md`](../src/README.md) にあります。

## 確かめ方

`writeup/` フォルダで:

```sh
python3 -I -B src/verify.py                    # 内部状態だけで 3 本を復号し、B を作り直す（鍵は使わない）
python3 -B src/check_key.py --khex 8ea590ce814594ad8cf5814549495349 --vhex d9e88b7e8be072ab26c423dbd356b3dde8c18df96a
                                               # 鍵 k と等価な鍵 v で 3 本を復号する（「総合: 当たり」）
```

## 未解決の問題: 鍵 v の元の文字列

鍵 v は、鍵ストリームからは 1 つに決まりません（条件が約 104〜106 ビットしかなく、同じ鍵ストリームを出す鍵が無数にある）。元の文字列を決めるには、鍵ストリームの外の手がかりが要ります。分かっていることと、調べた範囲・アイディアは [password.md](password.md) の §8〜§10 にまとめました。

- 鍵 k と同じ形（全角の語を「・」で区切り、`IISI` で終わる）ではありえない。
- 4・8・12…バイト目の和が小さい（漢字やひらがなが続く鍵ではなく、「・」・カタカナ・半角の文字が多い鍵）。
- 英小文字だけの鍵は 38 文字以下には存在しない。

候補を思いついたら、`python3 -B src/check_key.py --v "候補"` で一瞬で確かめられます。
