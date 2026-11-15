/* reads lines "hexkey" from stdin, prints mode warm S T and first 16 dwords of decrypt(zeros) = keystream-ish */
#include <stdio.h>
#include <stdlib.h>
#include "gccfast.h"
int main(void) {
    char line[512];
    while (fgets(line, sizeof line, stdin)) {
        uint8_t key[256]; int n = 0;
        for (char *p = line; p[0] && p[1] && p[0] != '\n'; p += 2) { unsigned v; sscanf(p, "%2x", &v); key[n++] = v; }
        gcc_state s;
        if (gcc_key_schedule(key, n, &s)) { printf("ERR\n"); continue; }
        uint32_t c[16] = {0}, p[16];
        gcc_decrypt_words(&s, c, p, 16);
        printf("%d %d %08x %08x", s.mode, s.warm, s.S, s.T);
        for (int i = 0; i < 16; i++) printf(" %08x", p[i]);
        printf("\n");
    }
    return 0;
}
