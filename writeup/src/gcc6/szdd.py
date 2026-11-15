import sys, struct, os
def expand(data):
    assert data[:8] == b"SZDD\x88\xf0\x27\x33", "not SZDD"
    mode = data[8]; missing = chr(data[9]); size = struct.unpack("<I", data[10:14])[0]
    win = bytearray(b" " * 4096); pos = 4096 - 16
    out = bytearray(); i = 14
    while i < len(data) and len(out) < size:
        ctrl = data[i]; i += 1
        for bit in range(8):
            if i >= len(data) or len(out) >= size: break
            if ctrl & (1 << bit):
                c = data[i]; i += 1
                out.append(c); win[pos] = c; pos = (pos + 1) & 0xFFF
            else:
                if i + 1 >= len(data): break
                b0, b1 = data[i], data[i+1]; i += 2
                mp = b0 | ((b1 & 0xF0) << 4); ml = (b1 & 0x0F) + 3
                for k in range(ml):
                    c = win[(mp + k) & 0xFFF]
                    out.append(c); win[pos] = c; pos = (pos + 1) & 0xFFF
    return bytes(out[:size]), missing, size
if __name__ == "__main__":
    src, dst = sys.argv[1], sys.argv[2]
    os.makedirs(dst, exist_ok=True)
    for fn in sorted(os.listdir(src)):
        p = os.path.join(src, fn)
        if not os.path.isfile(p): continue
        d = open(p, "rb").read()
        if d[:4] != b"SZDD": continue
        out, miss, size = expand(d)
        name = (fn[:-1] + miss) if (fn.endswith("_") and miss.isalnum()) else fn + ".out"
        open(os.path.join(dst, name), "wb").write(out)
        print(f"{fn:16s} -> {name:16s} {size:8d} {'OK' if len(out)==size else 'SHORT '+str(len(out))}")
