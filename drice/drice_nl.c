/* drice_nl.c — drice 的近无损扩展：闭环差分预测 + 逐样本死区量化
 *
 * 相对 drice.c 的改动：
 *   - 每通道每块有一个 NEAR 参数（0 = 无损；>0 = 死区宽度 δ）
 *   - 预测基准改为"上一重建值"（闭环预测），保证残差量化误差不累积
 *   - 死区量化：|r| ≤ δ → q = 0；否则 q = sign(r) * (|r| - δ)
 *   - 重建：x'[t] = prev + (q == 0 ? 0 : sign(q) * (|q| + δ))   = prev + sign(r)*|r| when |r|>δ
 *   - 帧头扩展：12B + reserved(4B) → 16B
 *       头字段：nch(4), ns(4), blk(4), near(4)           ← near 是 fallback 值
 *       后续：nch * (8B 长度表 + 4B per-ch NEAR)           ← per-ch NEAR 缺省 = near
 *
 * 块格式：与 drice.c 一致，但预测基准 prev_recon 替代了原始 x[t-1]。
 *
 * 编译: gcc -O3 -o drice_nl drice_nl.c
 * 用法: drice_nl c in.bin out.bin nch [blk] [near]
 *       drice_nl d in.bin out.bin
 *       drice_nl c ... --per-near=file   （每通道 NEAR 文件：nch 个 uint16 LSB）
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define ESC 0xFF
#define DEF_BLK 1024

/* ---------- 比特流 ---------- */
typedef struct { uint8_t *buf; size_t len; uint32_t acc; int nb; } BW;
static inline void bw_put(BW *w, uint32_t bit) {
    w->acc = (w->acc << 1) | (bit & 1);
    if (++w->nb == 8) { w->buf[w->len++] = (uint8_t)w->acc; w->acc = 0; w->nb = 0; }
}
static inline void bw_bits(BW *w, uint32_t v, int n) { for (int i = n - 1; i >= 0; i--) bw_put(w, (v >> i) & 1); }
static inline void bw_flush(BW *w) { while (w->nb) bw_put(w, 0); }
static inline void bw_byte(BW *w, uint8_t b) { w->buf[w->len++] = b; }

typedef struct { const uint8_t *buf; size_t pos; uint32_t acc; int nb; } BR;
static inline uint32_t br_get(BR *r) {
    if (r->nb == 0) { r->acc = r->buf[r->pos++]; r->nb = 8; }
    r->nb--; return (r->acc >> r->nb) & 1;
}
static inline uint32_t br_bits(BR *r, int n) { uint32_t v = 0; for (int i = 0; i < n; i++) v = (v << 1) | br_get(r); return v; }
static inline uint32_t br_unary(BR *r) { uint32_t q = 0; while (!br_get(r)) q++; return q; }

static inline uint16_t zigzag16(int16_t e) { return (uint16_t)(((uint16_t)e << 1) ^ (uint16_t)(e >> 15)); }
static inline int16_t unzigzag16(uint16_t z) { return (int16_t)((z >> 1) ^ (0u - (z & 1))); }

/* 死区量化 + 重建映射。
 * 返回值 q（量化残差），通过 *recon 写出重建增量。
 *   |r| ≤ δ : q = 0,  recon = 0        (误差=|r|≤δ)
 *   |r| > δ : q = sign(r)*(|r|-δ), recon = r (重建完全等于原始)
 */
static inline int32_t deadzone(int32_t r, int32_t delta, int32_t *recon) {
    if ((r <= delta) && (r >= -delta)) { *recon = 0; return 0; }
    int32_t s = (r > 0) ? 1 : -1;
    int32_t mag = (r > 0) ? r : -r;
    int32_t q = s * (mag - delta);
    *recon = r;
    return q;
}

/* ---------- 编码一个通道 ---------- */
static size_t encode_ch(const int16_t *x, size_t n, int blk, int32_t near, uint8_t *out) {
    BW w = { out, 0, 0, 0 };
    size_t nblk = (n + blk - 1) / blk;
    uint32_t *zz = malloc((size_t)blk * sizeof(uint32_t));
    for (size_t b = 0; b < nblk; b++) {
        size_t s = b * (size_t)blk, m = (n - s < (size_t)blk) ? n - s : (size_t)blk;
        /* 用重建值前样本模拟：闭环不影响首样本，且内部循环里我们将 prev 起点置为 raw prev，再每步更新到 recon */
        /* 直接计算"假定使用重建值的"残差序列（最简做法：本地 prev 直接 = raw 前一样本起，每步改为 recon 即可） */
        /* 为简单起见，每块重建基准用原始首样本起，块内闭环 */
        int16_t prev = x[s];
        /* 计算残差（含死区） */
        int32_t *res = malloc((size_t)blk * sizeof(int32_t));
        int32_t *rec = malloc((size_t)blk * sizeof(int32_t));
        for (size_t i = 1; i < m; i++) {
            int32_t r = (int32_t)x[s + i] - (int32_t)prev;
            int32_t q = deadzone(r, near, &rec[i]);
            res[i] = q;
            zz[i - 1] = zigzag16((int16_t)q);
            prev = (int16_t)((int32_t)prev + rec[i]);
        }
        /* k 搜索 */
        int bestk = 0; double bestbits = 1e18;
        for (int k = 0; k <= 15; k++) {
            double bits = 8.0 + 16.0;
            for (size_t i = 1; i < m; i++) bits += 1.0 + (double)(zz[i - 1] >> k) + k;
            if (bits < bestbits) { bestbits = bits; bestk = k; }
        }
        if (8.0 + 16.0 * (double)m < bestbits) {
            bw_byte(&w, ESC);
            for (size_t i = 0; i < m; i++) { uint16_t rv = (uint16_t)x[s + i]; bw_byte(&w, rv & 0xFF); bw_byte(&w, rv >> 8); }
            free(res); free(rec); continue;
        }
        bw_byte(&w, (uint8_t)((0u << 7) | bestk));   /* 仅一阶预测，bit7=0 保留 */
        { uint16_t rv = (uint16_t)x[s]; bw_byte(&w, rv & 0xFF); bw_byte(&w, rv >> 8); }
        for (size_t i = 1; i < m; i++) {
            uint32_t q = zz[i - 1] >> bestk;
            while (q--) bw_put(&w, 0);
            bw_put(&w, 1);
            if (bestk) bw_bits(&w, zz[i - 1] & ((1u << bestk) - 1), bestk);
        }
        bw_flush(&w);
        free(res); free(rec);
    }
    free(zz);
    return w.len;
}

/* ---------- 解码一个通道 ---------- */
static void decode_ch(const uint8_t *in, size_t n, int16_t *x, size_t ns, int blk, int32_t near) {
    size_t nblk = (ns + blk - 1) / blk, p = 0;
    for (size_t b = 0; b < nblk; b++) {
        size_t s = b * (size_t)blk, m = (ns - s < (size_t)blk) ? ns - s : (size_t)blk;
        uint8_t h = in[p++];
        if (h == ESC) { for (size_t i = 0; i < m; i++) { x[s + i] = (int16_t)(in[p] | (in[p + 1] << 8)); p += 2; } continue; }
        int k = h & 0x7F;
        x[s] = (int16_t)(in[p] | (in[p + 1] << 8)); p += 2;
        BR r = { in, p, 0, 0 };
        int16_t prev = x[s];
        for (size_t i = 1; i < m; i++) {
            uint32_t qun = br_unary(&r);
            uint32_t lo = k ? br_bits(&r, k) : 0;
            int16_t z = unzigzag16((uint16_t)((qun << k) | lo));
            /* 死区闭式重建：编码端 r>δ 时 q=sign(r)*(|r|-δ) → 重建 r' = sign(q)*(|q|+δ) = r（无损）
               |r|≤δ 时 q=0 → 重建 r' = 0（误差 = |r| ≤ δ） */
            int32_t val;
            if (z == 0) {
                val = (int32_t)prev;
            } else {
                int32_t recon = (z > 0) ? ((int32_t)z + near) : ((int32_t)z - near);
                val = (int32_t)prev + recon;
            }
            if (val > 32767) val = 32767; if (val < -32768) val = -32768;
            x[s + i] = (int16_t)val;
            prev = x[s + i];
        }
        p = r.pos;
    }
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s c in.bin out.bin nch [blk] [near] [--per-near=file] | %s d in.bin out.bin\n", argv[0], argv[0]); return 1; }
    char mode = argv[1][0];
    FILE *fi = fopen(argv[2], "rb");
    if (!fi) { perror(argv[2]); return 1; }

    /* 解析可选参数 --per-near= */
    const char *per_near_file = NULL;
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "--per-near=", 11) == 0) per_near_file = argv[i] + 11;
    }

    if (mode == 'c') {
        if (argc < 5) { fprintf(stderr, "need nch\n"); return 1; }
        int nch = atoi(argv[4]);
        int blk = DEF_BLK; int32_t near = 0;
        int blk_set = 0, near_set = 0;
        for (int i = 5; i < argc; i++) {
            if (argv[i][0] != '-') {
                if (!blk_set) { blk = atoi(argv[i]); blk_set = 1; }
                else if (!near_set) { near = atoi(argv[i]); near_set = 1; }
            }
        }
        /* 读 per-channel NEAR 表 */
        int32_t *near_c = calloc((size_t)nch, sizeof(int32_t));
        for (int c = 0; c < nch; c++) near_c[c] = near;
        if (per_near_file) {
            FILE *fp = fopen(per_near_file, "rb"); if (!fp) { perror(per_near_file); return 1; }
            for (int c = 0; c < nch; c++) { uint16_t v; if (fread(&v, 2, 1, fp) != 1) { fprintf(stderr, "near file short\n"); return 1; } near_c[c] = v; }
            fclose(fp);
        }

        fseek(fi, 0, SEEK_END); long fsz = ftell(fi); fseek(fi, 0, SEEK_SET);
        size_t ns = (size_t)fsz / (nch * 2);
        int16_t *x = malloc(fsz); fread(x, 1, fsz, fi); fclose(fi);

        FILE *fo = fopen(argv[3], "wb");
        uint32_t hdr[4] = { (uint32_t)nch, (uint32_t)ns, (uint32_t)blk, (uint32_t)near };
        fwrite(hdr, 4, 4, fo);
        size_t *lens = malloc((size_t)nch * sizeof(size_t));
        uint8_t **obufs = malloc((size_t)nch * sizeof(uint8_t *));
        size_t chcap = (size_t)fsz / nch + 65536;
        for (int c = 0; c < nch; c++) obufs[c] = malloc(chcap);
        for (int c = 0; c < nch; c++) {
            lens[c] = encode_ch(x + (size_t)c * ns, ns, blk, near_c[c], obufs[c]);
            fwrite(&lens[c], sizeof(size_t), 1, fo);
        }
        for (int c = 0; c < nch; c++) { uint16_t nv = (uint16_t)near_c[c]; fwrite(&nv, 2, 1, fo); }
        for (int c = 0; c < nch; c++) fwrite(obufs[c], 1, lens[c], fo);
        long osz = ftell(fo); fclose(fo);
        fprintf(stderr, "in=%ld out=%ld CR=%.4f (nch=%d ns=%zu blk=%d near=%d%s)\n",
            fsz, osz, (double)fsz / osz, nch, ns, blk, near, per_near_file ? " per-ch" : "");
        free(near_c);
    } else {
        uint32_t hdr[4];
        if (fread(hdr, 4, 4, fi) != 4) { fprintf(stderr, "bad header\n"); return 1; }
        int nch = hdr[0]; size_t ns = hdr[1]; int blk = hdr[2]; int32_t near = (int32_t)hdr[3];
        fseek(fi, 0, SEEK_END); long fsz = ftell(fi); fseek(fi, 0, SEEK_SET);
        uint8_t *in = malloc(fsz); fread(in, 1, fsz, fi); fclose(fi);
        size_t *lens = malloc((size_t)nch * sizeof(size_t));
        uint16_t *near_c = malloc((size_t)nch * sizeof(uint16_t));
        memcpy(lens, in + 16, (size_t)nch * sizeof(size_t));
        memcpy(near_c, in + 16 + (size_t)nch * sizeof(size_t), (size_t)nch * sizeof(uint16_t));
        long off = 16 + (long)nch * (sizeof(size_t) + sizeof(uint16_t));
        int16_t *x = malloc((size_t)nch * ns * 2);
        for (int c = 0; c < nch; c++) { decode_ch(in + off, lens[c], x + (size_t)c * ns, ns, blk, (int32_t)near_c[c]); off += (long)lens[c]; }
        FILE *fo = fopen(argv[3], "wb");
        fwrite(x, 2, (size_t)nch * ns, fo); fclose(fo);
        free(near_c);
    }
    return 0;
}