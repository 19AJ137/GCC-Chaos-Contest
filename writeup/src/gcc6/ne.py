import os as _os
HERE = _os.path.dirname(_os.path.abspath(__file__))
# NE（Win16）DLL の簡易パーサ：セグメント、再配置、エントリ、取り込み名
import struct
class NE:
    def __init__(self, path):
        d = self.d = open(path, 'rb').read()
        o = self.o = struct.unpack_from('<I', d, 0x3c)[0]
        assert d[o:o+2] == b'NE'
        h = struct.unpack_from('<2sBBHHIHHHHHHHHHHHHHHHHHHH', d, o)
        (self.sig, _, _, self.entoff, self.entlen, _, self.flags, self.autodata, self.heap, self.stack,
         self.csip_ip, self.csip_cs, self.sssp_sp, self.sssp_ss, self.nseg, self.nmod, self.nnres,
         self.segoff, self.resoff, self.resnoff, self.modoff, self.impoff, *_r) = h[:22] + h[22:]
        self.nonres_off = struct.unpack_from('<I', d, o + 0x2c)[0]
        self.align = struct.unpack_from('<H', d, o + 0x32)[0]
        self.segs = []
        for i in range(self.nseg):
            so, sl, fl, ms = struct.unpack_from('<HHHH', d, o + self.segoff + 8 * i)
            fo = so << self.align
            ln = sl or 0x10000
            ms = ms or 0x10000
            data = d[fo:fo + ln] if so else b''
            relocs = []
            if fl & 0x100:
                r = fo + ln
                n = struct.unpack_from('<H', d, r)[0]
                for k in range(n):
                    relocs.append(struct.unpack_from('<BBHHH', d, r + 2 + 8 * k))
            self.segs.append(dict(fileoff=fo, len=ln, flags=fl, minalloc=ms, data=data, relocs=relocs))
        # module names
        self.mods = []
        for i in range(self.nmod):
            p = struct.unpack_from('<H', d, o + self.modoff + 2 * i)[0]
            self.mods.append(self.impname(p))
        self.entries = self._entries()
        self.names = self._names(o + self.resnoff) + self._names(self.nonres_off)
    def impname(self, p):
        a = self.o + self.impoff + p
        n = self.d[a]
        return self.d[a+1:a+1+n].decode('latin1')
    def _names(self, a):
        out = []
        first = True
        while True:
            n = self.d[a]
            if n == 0: break
            s = self.d[a+1:a+1+n].decode('latin1'); ordn = struct.unpack_from('<H', self.d, a+1+n)[0]
            out.append((s, ordn)); a += 3 + n
        return out
    def _entries(self):
        a = self.o + self.entoff; end = a + self.entlen; ordn = 1; ent = {}
        while a < end:
            cnt, seg = self.d[a], self.d[a+1]; a += 2
            if cnt == 0: break
            for k in range(cnt):
                if seg == 0: pass
                elif seg == 0xff:
                    fl, _, sn, off = struct.unpack_from('<BHBH', self.d, a); a += 6; ent[ordn] = (sn, off, fl)
                else:
                    fl, off = struct.unpack_from('<BH', self.d, a); a += 3; ent[ordn] = (seg, off, fl)
                ordn += 1
            if seg == 0: ordn += 0  # 空きの番号
        return ent
if __name__ == '__main__':
    import sys
    n = NE(sys.argv[1])
    print('flags %#x autodata %d heap %#x stack %#x  CS:IP %d:%#x  segs %d mods %s' % (n.flags, n.autodata, n.heap, n.stack, n.csip_cs, n.csip_ip, n.nseg, n.mods))
    for i, s in enumerate(n.segs, 1):
        print('seg %d: off %#x len %#x flags %#x minalloc %#x relocs %d' % (i, s['fileoff'], s['len'], s['flags'], s['minalloc'], len(s['relocs'])))
    print('names', n.names)
    print('entries', {k: (v[0], hex(v[1]), v[2]) for k, v in n.entries.items()})
    from collections import Counter
    for i, s in enumerate(n.segs, 1):
        c = Counter()
        for (src, typ, off, a1, a2) in s['relocs']:
            kind = typ & 3
            if kind == 1: c[('imp', n.mods[a1-1], a2, src)] += 1
            elif kind == 2: c[('impname', n.mods[a1-1], n.impname(a2), src)] += 1
            elif kind == 3: c[('OSFIXUP', a1, src)] += 1
            else: c[('internal', a1, src)] += 1
        for k, v in sorted(c.items(), key=str): print('  seg', i, k, v)
