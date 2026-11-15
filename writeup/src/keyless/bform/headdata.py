"""headdata.py — 語 0〜131 の平文と 2 組の連鎖フラグ（ttp/ の言語モデルの探索の結果）を読む
  head132_lm.json    言語モデルが読んだもの（バイト 335・345 に読み違いがある）
  head132_fixed.json その 2 バイトを直したもの
  .pkl（ttp/fulldec2.py などの出力を同じ形にしたもの）も読める"""
import json, pickle


def load_head(path):
    """dict(flags={語: (f, f')}, MA=bytes, ME=bytes, known=[語, ...]) を返す"""
    if path.endswith('.json'):
        j = json.load(open(path, encoding='utf-8'))
        return dict(flags={int(n): tuple(v) for n, v in j['flags'].items()}, MA=bytes.fromhex(j['MA']), ME=bytes.fromhex(j['ME']), known=j['known'])
    return pickle.load(open(path, 'rb'))
