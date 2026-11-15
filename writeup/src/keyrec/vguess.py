# 鍵 v（または k）の候補を目標値に当てる。python3 vguess.py [--k] cand1 cand2 ... / 標準入力から 1 行 1 候補
import sys
from gks import key_schedule, h2d
TV = dict(S=0x3a9d3481, m=1, k780=h2d('3fe99b19be354000'), k7b0=h2d('3fcf8c8bdb300652'), lo=0xab0f, hi=0xabc5)
TK = dict(S=0xa678044d, m=2, k780=h2d('3fe500287933f358'), k7b0=h2d('3fc3ca0d1ca9ebf4'), lo=0x7ec1, hi=0x7eec)
def score(key, T):
    r = key_schedule(key)
    k780, k7b0 = r['k'][0], r['k'][2]
    d1 = min(abs(k780 - T['k780']), abs(k780 - (1 - T['k780']))) if T is TV else abs(k780 - T['k780'])
    d3 = abs(k7b0 - T['k7b0'])
    hit = (r['S'] == T['S']) + (r['m'] == T['m']) + (d1 < 1e-9) + (d3 < 1e-9) + (T['lo'] <= r['b0b1'] <= T['hi'])
    return hit, r, d1, d3
if __name__ == '__main__':
    args = sys.argv[1:]
    T = TV
    if args and args[0] == '--k': T = TK; args = args[1:]
    cands = args if args else [l.rstrip('\r\n') for l in sys.stdin if l.strip()]
    for s in cands:
        key = s.encode('cp932')
        hit, r, d1, d3 = score(key, T)
        print('%d/5 len %2d m %d S %08x b0b1 %04x k780 %.13f k7b0 %.13f | %s' % (hit, len(key), r['m'], r['S'], r['b0b1'], r['k'][0], r['k'][2], s))
