# 鍵（パスワード）の復元に使う C プログラム。
# 浮動小数点の計算順を DLL と同じにするため、FMA の縮約と fast-math は禁止（-ffp-contract=off -fno-fast-math）。
CC      ?= gcc
CFLAGS  ?= -O2 -march=native -ffp-contract=off -fno-fast-math -fopenmp
LDLIBS  = -lm
PROGS   = spenum spenum_count spjoin2 tail3 vres kmitm k6band mktable spjoin short/hiscan

all: $(PROGS)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

k6band: k6band.c
	$(CC) $(CFLAGS) -O3 -o $@ $< $(LDLIBS)

clean:
	rm -f $(PROGS)

.PHONY: all clean
