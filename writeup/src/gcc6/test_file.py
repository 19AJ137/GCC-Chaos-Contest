import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import struct, random
from unicorn import UC_HOOK_CODE
from unicorn.x86_const import *
from emu import Emu, HEAP
import gcc

class FileEmu(Emu):
    def file_cryption(self, key, data, encrypt=True):
        uc = self._fresh()
        files = {}
        state = {'fds': {}, 'next': 3}
        names = {HEAP+0x100: 'in', HEAP+0x200: 'out'}
        uc.mem_write(HEAP, key+b'\0'); uc.mem_write(HEAP+0x100, b'in\0'); uc.mem_write(HEAP+0x200, b'out\0')
        files['in'] = bytearray(data); 
        def arg(i): return struct.unpack('<I', bytes(uc.mem_read(uc.reg_read(UC_X86_REG_ESP)+4+4*i, 4)))[0]
        def ret(v):
            sp = uc.reg_read(UC_X86_REG_ESP); ra = struct.unpack('<I', bytes(uc.mem_read(sp,4)))[0]
            uc.reg_write(UC_X86_REG_EAX, v & 0xffffffff); uc.reg_write(UC_X86_REG_ESP, sp+4); uc.reg_write(UC_X86_REG_EIP, ra)
        def hook(uc, addr, size, ud):
            if addr == 0x10009130:   # _open(name, flags, mode)
                nm = names[arg(0)]; fd = state['next']; state['next'] += 1
                if nm == 'out': files['out'] = bytearray()
                state['fds'][fd] = [nm, 0]; ret(fd)
            elif addr == 0x10008e80: # _read(fd, buf, n)
                fd, buf, n = arg(0), arg(1), arg(2); nm, pos = state['fds'][fd]
                chunk = bytes(files[nm][pos:pos+n]); uc.mem_write(buf, chunk); state['fds'][fd][1] += len(chunk); ret(len(chunk))
            elif addr == 0x10008500: # _write(fd, buf, n)
                fd, buf, n = arg(0), arg(1), arg(2); files[state['fds'][fd][0]] += bytes(uc.mem_read(buf, n)); ret(n)
            elif addr == 0x10008dd0: # _eof(fd)
                fd = arg(0); nm, pos = state['fds'][fd]; ret(1 if pos >= len(files[nm]) else 0)
            elif addr == 0x100087c0: # _close
                ret(0)
            elif addr == 0x10002890: # errno ptr
                uc.mem_write(HEAP+0x300, struct.pack('<I', 0)); ret(HEAP+0x300)
        uc.hook_add(UC_HOOK_CODE, hook, begin=0x10002890, end=0x10009200)
        r = self.call(uc, self.exports['File_cryption'], [1 if encrypt else 2, HEAP+0x100, HEAP+0x200, HEAP])
        return bytes(files['out']), r & 0xffff

e = FileEmu(); rnd = random.Random(5); bad = 0
for t in range(60):
    key = bytes(rnd.choice(b'abcdefghijklmnopqrstuvwxyz0123456789') for _ in range(rnd.randrange(1,17)))
    data = bytes(rnd.randrange(256) for _ in range(rnd.choice([1024, 2048, 3000, 5000, 100, 1027])))
    ref, rc = e.file_cryption(key, data)
    mine = gcc.file_cryption(key, data)
    if ref != mine:
        bad += 1; print('MISMATCH', key, len(data), rc, len(ref), [i for i in range(min(len(ref),len(mine))) if ref[i]!=mine[i]][:3])
    else:
        dec, _ = e.file_cryption(key, ref, encrypt=False)
        assert dec == data and gcc.file_cryption(key, ref, encrypt=False) == data
print('file tests mismatches:', bad, 'rc', rc)
