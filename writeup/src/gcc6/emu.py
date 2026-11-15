"""Run the real Gcc6_32.dll code under Unicorn (x86-32) to get ground-truth outputs."""
import struct
import pefile
from unicorn import Uc, UC_ARCH_X86, UC_MODE_32, UC_HOOK_MEM_UNMAPPED, UC_HOOK_CODE
from unicorn.x86_const import *

import os
# 32bit 版の DLL（Gcc6_32.dll、1997-01-23 ビルド）。extract.py で software/net9729.exe から取り出す。環境変数 GCC6_DLL32 で別の場所を指定できる
DLL = os.environ.get('GCC6_DLL32', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dll', 'Gcc6_32.dll'))
BASE = 0x10000000
STACK = 0x00200000
HEAP = 0x00400000
STOP = 0x00100000


class Emu:
    def __init__(self, path=DLL, fpcw=0x027F):
        pe = pefile.PE(path)
        self.img = pe.get_memory_mapped_image()
        self.fpcw = fpcw
        self.exports = {e.name.decode(): BASE + e.address for e in pe.DIRECTORY_ENTRY_EXPORT.symbols}

    def _fresh(self):
        uc = Uc(UC_ARCH_X86, UC_MODE_32)
        size = (len(self.img) + 0xFFF) & ~0xFFF
        uc.mem_map(BASE, size)
        uc.mem_write(BASE, bytes(self.img))
        uc.mem_map(STACK, 0x100000)
        uc.mem_map(HEAP, 0x400000)
        uc.mem_map(STOP, 0x1000)
        uc.mem_write(STOP, b'\xf4')  # hlt

        def unmapped(uc, access, addr, size, value, ud):
            eip = uc.reg_read(UC_X86_REG_EIP)
            raise RuntimeError('unmapped access %x at eip %x' % (addr, eip))
        uc.hook_add(UC_HOOK_MEM_UNMAPPED, unmapped)
        uc.reg_write(UC_X86_REG_FPCW, self.fpcw)
        return uc

    def call(self, uc, func, args):
        sp = STACK + 0x100000 - 0x100
        for a in reversed(args):
            sp -= 4
            uc.mem_write(sp, struct.pack('<I', a & 0xFFFFFFFF))
        sp -= 4
        uc.mem_write(sp, struct.pack('<I', STOP))
        uc.reg_write(UC_X86_REG_ESP, sp)
        uc.emu_start(func, STOP, count=50_000_000)
        return uc.reg_read(UC_X86_REG_EAX)

    def str_cryption(self, key: bytes, data: bytes, encrypt=True, want_state=False):
        uc = self._fresh()
        kp, ip, op = HEAP, HEAP + 0x1000, HEAP + 0x200000
        uc.mem_write(kp, key + b'\0')
        uc.mem_write(ip, data)
        uc.mem_write(op, b'\0' * (len(data) + 8))
        self.call(uc, self.exports['Str_cryption'], [1 if encrypt else 2, ip, len(data), op, kp])
        out = bytes(uc.mem_read(op, len(data)))
        if want_state:
            return out, uc
        return out

    def key_state(self, key: bytes):
        uc = self._fresh()
        kp = HEAP
        uc.mem_write(kp, key + b'\0')
        self.call(uc, 0x100020e0, [kp])  # cdecl key schedule; returns S in st0
        g = lambda a: struct.unpack('<d', bytes(uc.mem_read(a, 8)))[0]
        d = lambda a: struct.unpack('<I', bytes(uc.mem_read(a, 4)))[0]
        return {
            'mode': d(0x1000c790), 'x': g(0x1000c750), 'r': g(0x1000c760), 'y': g(0x1000c768),
            'u': g(0x1000c710), 'v': g(0x1000c720), 'A': g(0x1000c770), 'C': g(0x1000c778),
            'b0': uc.mem_read(0x1000c740, 1)[0], 'b1': uc.mem_read(0x1000c741, 1)[0],
            'k780': g(0x1000c780), 'k788': g(0x1000c788), 'k7b0': g(0x1000c7b0), 'k7b8': g(0x1000c7b8),
        }


if __name__ == '__main__':
    e = Emu()
    k = b'password'
    print(e.key_state(k))
    ct = e.str_cryption(k, b'Hello, chaos world! 0123456789')
    print(ct.hex())
    print(e.str_cryption(k, ct, encrypt=False))
