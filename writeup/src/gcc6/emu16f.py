"""emu16 に FILE_CRYPTION 用のファイル入出力（DOS3Call = INT 21h 相当）と GlobalAlloc 系を足したもの。"""
import os as _os
# 16bit 版の DLL（gcc6.dll、1996-09-02）。extract.py で software/Net9718.exe か iisi9605.exe から取り出す。環境変数 GCC6_DLL16 で別の場所を指定できる
DLL16 = _os.environ.get('GCC6_DLL16', _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), 'dll', 'gcc6.dll'))
HERE = _os.path.dirname(_os.path.abspath(__file__))
import struct, sys
from unicorn import *
from unicorn.x86_const import *
sys.path.insert(0, HERE)
from emu16 import Emu16, PARA

class Emu16F(Emu16):
    HEAP0 = 0x7000
    def __init__(self, path):
        super().__init__(path)
        self.files = {}; self.fh = {}; self.nexth = 5; self.gnext = self.HEAP0; self.gsize = {}
        self.dos = {}
        self.uc.hook_add(UC_HOOK_INTR, self._intr)
    def _cstr(self, seg, off):
        a = (seg << 4) + off; s = b''
        while True:
            c = bytes(self.uc.mem_read(a + len(s), 1))
            if c == b'\0': return s.decode('latin1')
            s += c
    def _stubhook(self, uc, addr, size, ud):
        k = self.stubaddr.get(addr)
        if k is None: return
        self.calls[k] = self.calls.get(k, 0) + 1
        sp = uc.reg_read(UC_X86_REG_SP); ss = uc.reg_read(UC_X86_REG_SS)
        stk = lambda o, n=2: int.from_bytes(uc.mem_read((ss << 4) + sp + o, n), 'little')
        ip, cs = stk(0), stk(2)
        R = lambda r: uc.reg_read(r)
        argb = self.ARGB.get(k, 0)
        fl = uc.reg_read(UC_X86_REG_EFLAGS) & ~1
        if k == ('KERNEL', 102):                       # DOS3Call
            fl = self._dos(uc, fl)
        elif False:
            ah = R(UC_X86_REG_AH); al = R(UC_X86_REG_AL); bx = R(UC_X86_REG_BX); cx = R(UC_X86_REG_CX)
            dx = R(UC_X86_REG_DX); ds = R(UC_X86_REG_DS)
            self.dos[ah] = self.dos.get(ah, 0) + 1
            ax = 0; err = False
            if ah in (0x3c, 0x3d):
                name = self._cstr(ds, dx)
                if ah == 0x3c: self.files[name] = bytearray()
                if name not in self.files: ax, err = 2, True
                else:
                    h = self.nexth; self.nexth += 1; self.fh[h] = [name, 0]; ax = h
            elif ah == 0x3e:
                self.fh.pop(bx, None)
            elif ah == 0x3f:
                name, pos = self.fh[bx]; d = self.files[name][pos:pos + cx]
                uc.mem_write((ds << 4) + dx, bytes(d)); self.fh[bx][1] += len(d); ax = len(d)
            elif ah == 0x40:
                name, pos = self.fh[bx]; d = bytes(uc.mem_read((ds << 4) + dx, cx))
                f = self.files[name]; f[pos:pos + cx] = d; self.fh[bx][1] += cx; ax = cx
            elif ah == 0x42:
                name, pos = self.fh[bx]; off = (cx << 16) | dx
                if off >= 0x80000000: off -= 1 << 32
                base = {0: 0, 1: pos, 2: len(self.files[name])}[al]; npos = base + off
                self.fh[bx][1] = npos; ax = npos & 0xffff; uc.reg_write(UC_X86_REG_DX, (npos >> 16) & 0xffff)
            elif ah == 0x44:
                ax = 0; uc.reg_write(UC_X86_REG_DX, 0)     # ファイル（装置ではない）
            elif ah == 0x30:
                ax = 0x0005
            else:
                raise RuntimeError('DOS function %02x' % ah)
            uc.reg_write(UC_X86_REG_AX, ax & 0xffff)
            if err: fl |= 1
        elif k == ('KERNEL', 15):                      # GlobalAlloc(flags, size)
            size = stk(4, 4); seg = self.gnext; self.gnext += (size + 15) >> 4; self.gsize[seg] = size
            uc.mem_write(seg << 4, b'\0' * size); uc.reg_write(UC_X86_REG_AX, seg)
        elif k == ('KERNEL', 16):                      # GlobalReAlloc(h, size, flags)
            h = stk(8); size = stk(4, 4); seg = self.gnext; self.gnext += (size + 15) >> 4
            old = bytes(uc.mem_read(h << 4, min(self.gsize.get(h, 0), size))); uc.mem_write(seg << 4, old + b'\0' * (size - len(old)))
            self.gsize[seg] = size; uc.reg_write(UC_X86_REG_AX, seg)
        elif k == ('KERNEL', 18):                      # GlobalLock(h) -> DX:AX
            uc.reg_write(UC_X86_REG_AX, 0); uc.reg_write(UC_X86_REG_DX, stk(4))
        elif k == ('KERNEL', 20):
            s = self.gsize.get(stk(4), 0); uc.reg_write(UC_X86_REG_AX, s & 0xffff); uc.reg_write(UC_X86_REG_DX, s >> 16)
        elif k in (('KERNEL', 17), ('KERNEL', 19)):
            uc.reg_write(UC_X86_REG_AX, 0)
        elif k in (('KERNEL', 23), ('KERNEL', 24)):
            uc.reg_write(UC_X86_REG_AX, stk(4))
        elif k == ('WIN87EM', 1):
            uc.reg_write(UC_X86_REG_AX, 0)
        elif k == ('KERNEL', 3):
            uc.reg_write(UC_X86_REG_AX, 0x0A03)
        else:
            raise RuntimeError('unexpected import call %s' % (k,))
        uc.reg_write(UC_X86_REG_EFLAGS, fl)
        uc.reg_write(UC_X86_REG_SP, sp + 4 + argb)
        uc.reg_write(UC_X86_REG_CS, cs); uc.reg_write(UC_X86_REG_IP, ip)
    def _dos(self, uc, fl):
        R = lambda r: uc.reg_read(r)
        ah = R(UC_X86_REG_AH); al = R(UC_X86_REG_AL); bx = R(UC_X86_REG_BX); cx = R(UC_X86_REG_CX)
        dx = R(UC_X86_REG_DX); ds = R(UC_X86_REG_DS)
        self.dos[ah] = self.dos.get(ah, 0) + 1
        ax = 0; err = False
        if ah in (0x3c, 0x3d):
            name = self._cstr(ds, dx)
            if ah == 0x3c: self.files[name] = bytearray()
            if name not in self.files: ax, err = 2, True
            else:
                h = min(x for x in range(5, 64) if x not in self.fh); self.fh[h] = [name, 0]; ax = h
        elif ah == 0x3e:
            self.fh.pop(bx, None)
        elif ah == 0x3f:
            name, pos = self.fh[bx]; d = self.files[name][pos:pos + cx]
            uc.mem_write((ds << 4) + dx, bytes(d)); self.fh[bx][1] += len(d); ax = len(d)
        elif ah == 0x40:
            name, pos = self.fh[bx]; d = bytes(uc.mem_read((ds << 4) + dx, cx))
            f = self.files[name]; f[pos:pos + cx] = d; self.fh[bx][1] += cx; ax = cx
        elif ah == 0x42:
            name, pos = self.fh[bx]; off = (cx << 16) | dx
            if off >= 0x80000000: off -= 1 << 32
            base = {0: 0, 1: pos, 2: len(self.files[name])}[al]; npos = base + off
            self.fh[bx][1] = npos; ax = npos & 0xffff; uc.reg_write(UC_X86_REG_DX, (npos >> 16) & 0xffff)
        elif ah == 0x44:
            ax = 0; uc.reg_write(UC_X86_REG_DX, 0)
        elif ah == 0x30:
            ax = 0x0005
        elif ah == 0x43:                                # 属性の取得・設定
            name = self._cstr(ds, dx)
            if name in self.files: uc.reg_write(UC_X86_REG_CX, 0x20)
            else: ax, err = 2, True
        elif ah == 0x41:                                # 削除
            name = self._cstr(ds, dx)
            if self.files.pop(name, None) is None: ax, err = 2, True
        else:
            raise RuntimeError('DOS function %02x' % ah)
        if not err and ah == 0x43:
            ax = R(UC_X86_REG_AX)
        uc.reg_write(UC_X86_REG_AX, ax & 0xffff)
        return fl | 1 if err else fl & ~1
    def _intr(self, uc, intno, ud):
        if intno != 0x21:
            raise RuntimeError('int %#x at %04x:%04x' % (intno, uc.reg_read(UC_X86_REG_CS), uc.reg_read(UC_X86_REG_IP)))
        fl = uc.reg_read(UC_X86_REG_EFLAGS)
        uc.reg_write(UC_X86_REG_EFLAGS, self._dos(uc, fl))
    def file_cryption(self, key, data, mode=1, cw=0x037F, count=200_000_000):
        uc = self.uc
        self.files = {'C:\\IN.BIN': bytearray(data)}; self.fh = {}; self.gnext = self.HEAP0
        bufl = PARA['buf'] << 4
        IN, OUT, KEY = 0x0000, 0x0100, 0x0200
        uc.mem_write(bufl + IN, b'C:\\IN.BIN\0'); uc.mem_write(bufl + OUT, b'C:\\OUT.BIN\0'); uc.mem_write(bufl + KEY, key + b'\0')
        drv = PARA['drv'] << 4
        sn, eo, _ = self.ne.entries[2]
        push = lambda v: b'\x68' + struct.pack('<H', v & 0xffff)
        code = b'\xdb\xe3' + b'\x2e\xd9\x2e' + struct.pack('<H', 0x100)
        code += push(mode) + push(PARA['buf']) + push(IN) + push(PARA['buf']) + push(OUT) + push(PARA['buf']) + push(KEY)
        code += b'\x9a' + struct.pack('<HH', eo, self.segpara[sn])
        end = len(code); code += b'\xf4'
        uc.mem_write(drv, code); uc.mem_write(drv + 0x100, struct.pack('<HH', cw, 0))
        uc.reg_write(UC_X86_REG_SS, PARA['stack']); uc.reg_write(UC_X86_REG_SP, 0xFFF0)
        uc.reg_write(UC_X86_REG_DS, PARA['data']); uc.reg_write(UC_X86_REG_ES, PARA['data'])
        uc.reg_write(UC_X86_REG_CS, PARA['drv'])
        uc.emu_start(drv, drv + end, count=count)
        return bytes(self.files.get('C:\\OUT.BIN', b'')), uc.reg_read(UC_X86_REG_AX)

if __name__ == '__main__':
    sys.path.insert(0, HERE)
    import gcc, random
    ARGB = Emu16.ARGB
    e = Emu16F(DLL16)
    rnd = random.Random(3)
    bad = 0; tot = 0
    for n in (1, 3, 4, 7, 1023, 1024, 1025, 1338, 2048, 2051, 3072):
        for key in (b'testkey1', b'123', b'gao', b'n2j8c0twldod'):
            data = bytes(rnd.randrange(256) for _ in range(n))
            for mode, enc in ((1, True), (2, False)):
                out, ax = e.file_cryption(key, data, mode=mode)
                ref = gcc.file_cryption(key, data, enc); tot += 1
                if out != ref: bad += 1; print('DIFF', n, key, mode, len(out), ax)
    print('cases', tot, 'mismatch', bad, 'dos', e.dos, 'calls', e.calls)
