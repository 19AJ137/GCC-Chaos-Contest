# 鍵なしの解読（bform/ の写像 2 の状態復元）に使う C プログラム。
# 浮動小数点の計算順を DLL と同じにするため、FMA の縮約と fast-math は禁止（-ffp-contract=off -fno-fast-math）。
CC      ?= gcc
CFLAGS  ?= -O2 -march=native -ffp-contract=off -fno-fast-math
LDLIBS  = -lm
PROGS   = bform/ksolve2 bform/kexact2

all: $(PROGS)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(PROGS)

.PHONY: all clean
