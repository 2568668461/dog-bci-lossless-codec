/* drice.c — Delta + 自适应 Golomb-Rice 无损压缩器（CCSDS-121 / Neuralink 挑战赛风格参考模型）
 *
 * 面向 FPGA 移植设计：
 *   - 逐通道独立编码（无跨通道依赖，通道可完全并行）
 *   - 一阶/二阶差分预测（1 bit 选择）+ Zigzag + Golomb-Rice
 *   - 每 BLK 个样本一个块，块头 1 字节 + 起始样本原始值，块间字节对齐
 *   - 无查找表、无乘除法（k 用移位实现）、无跨块状态 → 流式天然友好
 *
 * 块格式（每通道每块）：
 *   header 1 字节: 0xFF = 本块 16bit 原样转义
 *                  否则 bit7 = 预测器(0:一阶差分, 1:二阶差分), bit0..6 = k (0..15)
 *   delta1: 2 字节起始样本; delta2: 4 字节两个起始样本
 *   随后样本的 GR 码（MSB-first 比特流），块尾 flush 到字节边界
 *
 * 文件格式: 12 字节头(nch, ns, blk, reserved) + nch*8 长度表 + 各通道码流
 * 编译: gcc -O3 -o drice drice.c
 * 用法: drice c in.bin out.bin nch [blk]
 *       drice d in.bin out.bin        (参数在文件头里)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define ESC 0xFF
#define DEF_BLK 1024

/* ---------- 比特流（直接写底层缓冲） ---------- */
typedef struct { uint8_t *buf; size_t len; uint32_t acc; int nb; } BW;
static inline void bw_put(BW *w, uint32_t bit) {
    w->acc = (w->acc << 1) | (bit & 1);
    if (++w->nb == 8) { w->buf[w->len++] = (uint8_t)w->acc; w->acc = 0; w->nb = 0; }
}
static inline void bw_bits(BW *w, uint32_t v, int n) { for (int i = n - 1; i >= 0; i--) bw_put(w, (v >> i) & 1); }
static inline void bw_flush(BW *w) { while (w->nb) bw_put(w, 0); }
/* 块首字节对齐写入（调用时保证 nb==0） */
static inline void bw_byte(BW *w, uint8_t b) { w->buf[w->len++] = b; }

/* ---------- 比特流读 ---------- */
typedef struct { const uint8_t *buf; size_t pos; uint32_t acc; int nb; } BR;
static inline uint32_t br_get(BR *r) {
    if (r->nb == 0) { r->acc = r->buf[r->pos++]; r->nb = 8; }
    r->nb--; return (r->acc >> r->nb) & 1;
}
static inline uint32_t br_bits(BR *r, int n) { uint32_t v = 0; for (int i = 0; i < n; i++) v = (v << 1) | br_get(r); return v; }
static inline uint32_t br_unary(BR *r) { uint32_t q = 0; while (!br_get(r)) q++; return q; }

static inline uint16_t zigzag16(int16_t e) { return (uint16_t)(((uint16_t)e << 1) ^ (uint16_t)(e >> 15)); }
static inline int16_t unzigzag16(uint16_t z) { return (int16_t)((z >> 1) ^ (0u - (z & 1))); }

/* ---------- 编码一个通道 ---------- */
static size_t encode_ch(const int16_t *x, size_t n, int blk, uint8_t *out) {
    BW w = { out, 0, 0, 0 };
    size_t nblk = (n + blk - 1) / blk;
    uint32_t *zz = malloc((size_t)blk * sizeof(uint32_t));
    for (size_t b = 0; b < nblk; b++) {
        size_t s = b * (size_t)blk, m = (n - s < (size_t)blk) ? n - s : (size_t)blk;

        /* 预测器选择: 比较一阶/二阶残差平均幅度 */
        double m1 = 0, m2 = 0; int pred = 0; int nraw, nres;
        for (size_t i = 1; i < m; i++) m1 += abs((int)x[s + i] - (int)x[s + i - 1]);
        m1 /= (m > 1) ? (double)(m - 1) : 1.0;
        if (m > 2) {
            for (size_t i = 2; i < m; i++) { int p = 2 * (int)x[s + i - 1] - (int)x[s + i - 2]; m2 += abs((int)x[s + i] - p); }
            m2 /= (double)(m - 2);
            if (m2 < m1 * 0.95) pred = 1;
        }
        if (pred == 0) {
            nraw = 1; nres = (int)m - 1;
            for (size_t i = 1; i < m; i++) zz[i - 1] = zigzag16((int16_t)((int)x[s + i] - (int)x[s + i - 1]));
        } else {
            nraw = 2; nres = (int)m - 2;
            for (size_t i = 2; i < m; i++) { int p = 2 * (int)x[s + i - 1] - (int)x[s + i - 2]; zz[i - 2] = zigzag16((int16_t)((int)x[s + i] - p)); }
        }

        /* k 完整搜索（0..15），并与整块转义比较 */
        int bestk = 0; double bestbits = 1e18;
        for (int k = 0; k <= 15; k++) {
            double bits = 8.0 + 16.0 * nraw;
            for (int i = 0; i < nres; i++) bits += 1.0 + (double)(zz[i] >> k) + k;
            if (bits < bestbits) { bestbits = bits; bestk = k; }
        }
        if (8.0 + 16.0 * (double)m < bestbits) { /* 转义块 */
            bw_byte(&w, ESC);
            for (size_t i = 0; i < m; i++) { uint16_t rv = (uint16_t)x[s + i]; bw_byte(&w, rv & 0xFF); bw_byte(&w, rv >> 8); }
            continue;
        }
        bw_byte(&w, (uint8_t)((pred << 7) | bestk));
        for (int i = 0; i < nraw; i++) { uint16_t rv = (uint16_t)x[s + i]; bw_byte(&w, rv & 0xFF); bw_byte(&w, rv >> 8); }
        for (int i = 0; i < nres; i++) {
            uint32_t q = zz[i] >> bestk;
            while (q--) bw_put(&w, 0);
            bw_put(&w, 1);
            if (bestk) bw_bits(&w, zz[i] & ((1u << bestk) - 1), bestk);
        }
        bw_flush(&w);
    }
    free(zz);
    return w.len;
}

/* ---------- 解码一个通道 ---------- */
static void decode_ch(const uint8_t *in, size_t n, int16_t *x, size_t ns, int blk) {
    size_t nblk = (ns + blk - 1) / blk, p = 0;
    for (size_t b = 0; b < nblk; b++) {
        size_t s = b * (size_t)blk, m = (ns - s < (size_t)blk) ? ns - s : (size_t)blk;
        uint8_t h = in[p++];
        if (h == ESC) { for (size_t i = 0; i < m; i++) { x[s + i] = (int16_t)(in[p] | (in[p + 1] << 8)); p += 2; } continue; }
        int pred = (h >> 7) & 1, k = h & 0x7F;
        if (pred == 0) {
            x[s] = (int16_t)(in[p] | (in[p + 1] << 8)); p += 2;
            BR r = { in, p, 0, 0 };
            for (size_t i = 1; i < m; i++) {
                uint32_t q = br_unary(&r);
                uint32_t lo = k ? br_bits(&r, k) : 0;
                uint16_t z = (uint32_t)((q << k) | lo);
                x[s + i] = (int16_t)((int)x[s + i - 1] + (int)unzigzag16(z));
            }
            p = r.pos;
        } else {
            x[s] = (int16_t)(in[p] | (in[p + 1] << 8)); x[s + 1] = (int16_t)(in[p + 2] | (in[p + 3] << 8)); p += 4;
            BR r = { in, p, 0, 0 };
            for (size_t i = 2; i < m; i++) {
                uint32_t q = br_unary(&r);
                uint32_t lo = k ? br_bits(&r, k) : 0;
                uint16_t z = (uint32_t)((q << k) | lo);
                int pp = 2 * (int)x[s + i - 1] - (int)x[s + i - 2];
                x[s + i] = (int16_t)(pp + (int)unzigzag16(z));
            }
            p = r.pos;
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s c in.bin out.bin nch [blk] | %s d in.bin out.bin\n", argv[0], argv[0]); return 1; }
    char mode = argv[1][0];
    FILE *fi = fopen(argv[2], "rb");
    if (!fi) { perror(argv[2]); return 1; }

    if (mode == 'c') {
        if (argc < 5) { fprintf(stderr, "need nch\n"); return 1; }
        int nch = atoi(argv[4]);
        int blk = (argc > 5) ? atoi(argv[5]) : DEF_BLK;
        fseek(fi, 0, SEEK_END); long fsz = ftell(fi); fseek(fi, 0, SEEK_SET);
        size_t ns = (size_t)fsz / (nch * 2);
        int16_t *x = malloc(fsz); fread(x, 1, fsz, fi); fclose(fi);

        FILE *fo = fopen(argv[3], "wb");
        uint32_t hdr[3] = { (uint32_t)nch, (uint32_t)ns, (uint32_t)blk };
        fwrite(hdr, 4, 3, fo); fwrite("\0\0\0\0", 1, 4, fo);
        size_t *lens = malloc((size_t)nch * sizeof(size_t));
        uint8_t **obufs = malloc((size_t)nch * sizeof(uint8_t *));
        size_t chcap = (size_t)fsz / nch + 65536;
        for (int c = 0; c < nch; c++) obufs[c] = malloc(chcap);
        for (int c = 0; c < nch; c++) {
            lens[c] = encode_ch(x + (size_t)c * ns, ns, blk, obufs[c]);
            fwrite(&lens[c], sizeof(size_t), 1, fo);
        }
        for (int c = 0; c < nch; c++) fwrite(obufs[c], 1, lens[c], fo);
        long osz = ftell(fo); fclose(fo);
        fprintf(stderr, "in=%ld out=%ld CR=%.4f (nch=%d ns=%zu blk=%d)\n", fsz, osz, (double)fsz / osz, nch, ns, blk);
    } else {
        uint32_t hdr[3]; uint32_t pad;
        if (fread(hdr, 4, 3, fi) != 3 || fread(&pad, 4, 1, fi) != 1) { fprintf(stderr, "bad header\n"); return 1; }
        int nch = hdr[0]; size_t ns = hdr[1]; int blk = hdr[2];
        fseek(fi, 0, SEEK_END); long fsz = ftell(fi); fseek(fi, 0, SEEK_SET);
        uint8_t *in = malloc(fsz); fread(in, 1, fsz, fi); fclose(fi);
        size_t *lens = malloc((size_t)nch * sizeof(size_t));
        memcpy(lens, in + 16, (size_t)nch * sizeof(size_t));
        long off = 16 + nch * (long)sizeof(size_t);
        int16_t *x = malloc((size_t)nch * ns * 2);
        for (int c = 0; c < nch; c++) { decode_ch(in + off, lens[c], x + (size_t)c * ns, ns, blk); off += (long)lens[c]; }
        FILE *fo = fopen(argv[3], "wb");
        fwrite(x, 2, (size_t)nch * ns, fo); fclose(fo);
    }
    return 0;
}
