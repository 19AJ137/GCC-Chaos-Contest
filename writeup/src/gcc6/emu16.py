"""16bit 版 gcc6.dl_（NE、Win16）を Unicorn の実モードで直接動かす。

- セグメントを段落境界に置き、セレクタの再配置には段落番号を書く（実モードなのでセレクタ＝段落）。
- 取り込み関数（KERNEL・WIN87EM）は専用の段の番地に割り当て、呼ばれたら Python で処理する。
- 呼び出す前に fninit; fldcw で FPU 制御語を設定する（DLL 自身は精度を設定しない）。
STR_CRYPTION(mode, in:far, len, out:far, key:far) は Pascal 呼び出し（引数 16 バイト、callee が片付ける）。
"""
import os as _os
# 16bit 版の DLL（gcc6.dll、1996-09-02）。extract.py で software/Net9718.exe か iisi9605.exe から取り出す。環境変数 GCC6_DLL16 で別の場所を指定できる
DLL16 = _os.environ.get('GCC6_DLL16', _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'dll', 'gcc6.dll'))
HERE = _os.path.dirname(_os.path.abspath(__file__))
import struct, sys
from unicorn import *
from unicorn.x86_const import *
sys.path.insert(0, HERE)
from ne import NE

PARA = {'code': 0x1000, 'data': 0x2000, 'stack': 0x3000, 'buf': 0x4000, 'stub': 0x5000, 'drv': 0x6000}

class Emu16:
    def __init__(self, path, trace=False):
        self.ne = n = NE(path)
        self.trace = trace
        uc = self.uc = Uc(UC_ARCH_X86, UC_MODE_16)
        uc.mem_map(0, 0x100000)
        self.segpara = {1: PARA['code'], 2: PARA['data']}
        for i, s in enumerate(n.segs, 1):
            uc.mem_write(self.segpara[i] << 4, s['data'])
        self.stubs = {}       # (mod, ord) -> stub offset
        self.stubaddr = {}    # linear -> (mod, ord)
        self._apply_relocs()
        self.calls = {}
        uc.hook_add(UC_HOOK_CODE, self._stubhook, begin=PARA['stub'] << 4, end=(PARA['stub'] << 4) + 0xffff)
    # ----------------------------------------------------------- relocations
    def _stub(self, mod, ordn):
        k = (mod, ordn)
        if k not in self.stubs:
            off = 0x10 * (len(self.stubs) + 1)
            self.stubs[k] = off
            self.stubaddr[(PARA['stub'] << 4) + off] = k
            self.uc.mem_write((PARA['stub'] << 4) + off, b'\x90' * 0x10)
        return PARA['stub'], self.stubs[k]
    CONST = {('KERNEL', 113): 12, ('KERNEL', 114): 0x1000, ('KERNEL', 178): 0x0400 | 0x0008}  # __AHSHIFT, __AHINCR, __WINFLAGS(80x87|386)
    def _apply_relocs(self):
        n, uc = self.ne, self.uc
        self.osfix = []
        for i, s in enumerate(n.segs, 1):
            base = self.segpara[i] << 4
            for (src, typ, off, a1, a2) in s['relocs']:
                kind = typ & 3; additive = bool(typ & 4)
                if kind == 3:
                    self.osfix.append((i, off, a1)); continue
                if kind == 0:
                    if a1 == 0xff:
                        sn, eo, _ = n.entries[a2]; seg, val = self.segpara[sn], eo
                    else:
                        seg, val = self.segpara[a1], a2
                elif kind == 1:
                    mod = n.mods[a1 - 1]
                    if (mod, a2) in self.CONST:
                        seg, val = 0, self.CONST[(mod, a2)]
                    else:
                        seg, val = self._stub(mod, a2)
                else:
                    raise NotImplementedError('import by name')
                # 再配置の鎖をたどる（加算型でなければ、書き込み先の語が次の位置）
                p = off
                while True:
                    a = base + p
                    cur = struct.unpack('<H', uc.mem_read(a + (2 if src == 3 else 0), 2) if src != 3 else uc.mem_read(a, 2))[0]
                    if src == 2:      # selector
                        newv = struct.pack('<H', (seg + (cur if additive else 0)) & 0xffff) if kind != 1 or (mod, a2) not in self.CONST else struct.pack('<H', val)
                        if kind == 1 and (mod, a2) in self.CONST: newv = struct.pack('<H', val)
                        uc.mem_write(a, newv)
                    elif src == 3:    # far pointer
                        uc.mem_write(a, struct.pack('<HH', (val + (cur if additive else 0)) & 0xffff, seg))
                    elif src == 5:    # offset
                        uc.mem_write(a, struct.pack('<H', (val + (cur if additive else 0)) & 0xffff))
                    elif src == 0:    # low byte
                        uc.mem_write(a, bytes([val & 0xff]))
                    else:
                        raise NotImplementedError('src type %d' % src)
                    if additive or cur == 0xffff or src == 0:
                        break
                    p = cur
    # ----------------------------------------------------------- imports
    ARGB = {('KERNEL', 1): 2, ('KERNEL', 3): 0, ('KERNEL', 4): 6, ('KERNEL', 15): 6, ('KERNEL', 16): 8,
            ('KERNEL', 17): 2, ('KERNEL', 18): 2, ('KERNEL', 19): 2, ('KERNEL', 20): 2, ('KERNEL', 23): 2,
            ('KERNEL', 24): 2, ('KERNEL', 48): 2, ('KERNEL', 131): 0, ('KERNEL', 137): 6, ('KERNEL', 91): 0,
            ('KERNEL', 102): 0, ('WIN87EM', 1): 0}
    def _stubhook(self, uc, addr, size, ud):
        k = self.stubaddr.get(addr)
        if k is None: return
        self.calls[k] = self.calls.get(k, 0) + 1
        sp = uc.reg_read(UC_X86_REG_SP); ss = uc.reg_read(UC_X86_REG_SS)
        ip, cs = struct.unpack('<HH', uc.mem_read((ss << 4) + sp, 4))
        ret_ax = 0
        if k == ('WIN87EM', 1):
            ret_ax = 0  # __fpmath: 何もしない（FPU は呼び出し前に設定済み）
        elif k == ('KERNEL', 3):
            ret_ax = 0x0A03
        else:
            raise RuntimeError('unexpected import call %s' % (k,))
        uc.reg_write(UC_X86_REG_AX, ret_ax)
        uc.reg_write(UC_X86_REG_SP, sp + 4 + self.ARGB.get(k, 0))
        uc.reg_write(UC_X86_REG_CS, cs); uc.reg_write(UC_X86_REG_IP, ip)
    # ----------------------------------------------------------- call
    def str_cryption(self, key: bytes, data: bytes, mode=1, cw=0x037F, count=50_000_000):
        uc = self.uc
        bufl = PARA['buf'] << 4
        IN, OUT, KEY = 0x0000, 0x8000, 0xF000
        uc.mem_write(bufl + IN, data)
        uc.mem_write(bufl + OUT, b'\0' * (len(data) + 8))
        uc.mem_write(bufl + KEY, key + b'\0')
        drv = PARA['drv'] << 4
        sn, eo, _ = self.ne.entries[4]
        def push(v): return b'\x68' + struct.pack('<H', v & 0xffff)
        code = b'\xdb\xe3'                                  # fninit
        code += b'\x2e\xd9\x2e' + struct.pack('<H', 0x100)  # fldcw cs:[0x100]
        code += push(mode) + push(PARA['buf']) + push(IN) + push(len(data)) + push(PARA['buf']) + push(OUT) + push(PARA['buf']) + push(KEY)
        code += b'\x9a' + struct.pack('<HH', eo, self.segpara[sn])   # call far
        code += b'\x2e\xd9\x3e' + struct.pack('<H', 0x102)  # fnstcw cs:[0x102]
        end = len(code)
        code += b'\xf4'
        uc.mem_write(drv, code)
        uc.mem_write(drv + 0x100, struct.pack('<HH', cw, 0))
        uc.reg_write(UC_X86_REG_SS, PARA['stack']); uc.reg_write(UC_X86_REG_SP, 0xFFF0)
        uc.reg_write(UC_X86_REG_DS, PARA['data']); uc.reg_write(UC_X86_REG_ES, PARA['data'])
        uc.reg_write(UC_X86_REG_CS, PARA['drv'])
        uc.emu_start(drv, drv + end, count=count)
        ax = uc.reg_read(UC_X86_REG_AX)
        return bytes(uc.mem_read(bufl + OUT, len(data))), ax

if __name__ == '__main__':
    sys.path.insert(0, HERE)
    import gcc, random
    e = Emu16(DLL16)
    print('osfixups', len(e.osfix), 'stubs', e.stubs)
    data = bytes(range(256)) * 4
    for key in (b'testkey1', b'123', b'gao'):
        out, ax = e.str_cryption(key, data)
        ref = gcc.str_cryption(key, data)
        print(key, 'ax', ax, 'match gcc.py' if out == ref else 'DIFF at %d' % next(i for i in range(len(out)) if out[i] != ref[i]), 'calls', e.calls)
