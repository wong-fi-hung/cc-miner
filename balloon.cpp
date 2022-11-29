#include <stdbool.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <openssl/aes.h>
#include <openssl/evp.h>
#include <openssl/sha.h>
#include "balloon.h"
#include <malloc.h>
#include "cuda_helper.h"
#include "miner.h"

//#ifdef __cplusplus
//extern "C" {
//#endif
static uint32_t *d_KNonce[MAX_GPUS];
__constant__ uint32_t pTarget[8];
__constant__ uint2 c_PaddedMessage80[10];

#define TPB 512
#define NONCES_PER_THREAD 32

static void balloon_init(struct balloon_options *opts, int64_t s_cost, int32_t t_cost) {
	opts->s_cost = s_cost;
	opts->t_cost = t_cost;
}

void balloon_hash(unsigned char *input, unsigned char *output) {
	//balloon(input, output);
}

void balloon(int thr_id, uint32_t threads, const uint8_t* input, uint8_t* output, uint32_t *h_nounce) {
		
    struct balloon_options opts;
	struct hash_state s;
	
	balloon_init(&opts, (int64_t)128, (int32_t)4);
	hash_state_init(&s, &opts, input);             
	hash_state_fill(&s, input, 80);
	hash_state_mix(&s);
	uint8_t *b = block_index(&s, 4095);
	memcpy((char *)output, (const char *)b, 32);
	hash_state_free(&s);
}

__host__
void balloon_128_cuda(int thr_id, uint32_t threads, const uint8_t* input, uint8_t* output, uint32_t *h_nounce) {

	CUDA_SAFE_CALL(cudaMemsetAsync(d_KNonce[thr_id], 0xff, 2 * sizeof(uint32_t), gpustream[thr_id]));
	const uint32_t threadsperblock = 512;
	struct balloon_options opts;
	struct hash_state s;
	dim3 grid((threads + TPB*NONCES_PER_THREAD - 1) / TPB / NONCES_PER_THREAD);
	dim3 block(TPB);


	balloon_init(&opts, (int64_t)128, (int32_t)4); //ok
	hash_state_init(&s, &opts, input); //ok 

	hash_state_fill_cuda(&s, input, 80);
	hash_state_mix_cuda(&s);
	uint8_t *b = block_index(&s, 4095);
	memcpy((char *)output, (const char *)b, 32);
	hash_state_free(&s);
	CUDA_SAFE_CALL(cudaDeviceSynchronize());
	CUDA_SAFE_CALL(cudaMemcpy(h_nounce, d_KNonce[thr_id], 2 * sizeof(uint32_t), cudaMemcpyDeviceToHost));
}


static inline void bitstream_init(struct bitstream *b) {
	SHA256_Init(&b->c);
#if   OPENSSL_VERSION_NUMBER >= 0x10100000L
	b->ctx = EVP_CIPHER_CTX_new();
	EVP_CIPHER_CTX_init(b->ctx);
#else
	//EVP_CIPHER_CTX_init(&b->ctx);
#endif
	b->zeros = (uint8_t*)malloc(512);
	memset(b->zeros, 0, 512);
}

static inline void bitstream_free(struct bitstream *b) {
	uint8_t out[16];
	int outl;
#if   OPENSSL_VERSION_NUMBER >= 0x10100000L
	EVP_EncryptFinal(b->ctx, out, &outl);
	EVP_CIPHER_CTX_free(b->ctx);
#else
	//EVP_EncryptFinal(&b->ctx, out, &outl);
	//EVP_CIPHER_CTX_cleanup(&b->ctx);
#endif
	free(b->zeros);
}

static inline void bitstream_seed_add(struct bitstream *b, const void *seed, size_t seedlen) {
	SHA256_Update(&b->c, seed, seedlen);
}

static inline void bitstream_seed_finalize(struct bitstream *b) {
	uint8_t key_bytes[32];
	SHA256_Final(key_bytes, &b->c);
	uint8_t iv[16];
	memset(iv, 0, 16);
#if   OPENSSL_VERSION_NUMBER >= 0x10100000L
	EVP_EncryptInit(b->ctx, EVP_aes_128_ctr(), key_bytes, iv);
#else
	//	EVP_EncryptInit(&b->ctx, EVP_aes_128_ctr(), key_bytes, iv);
#endif
}

void bitstream_fill_buffer(struct bitstream *b, void *out, size_t outlen) {
	int encl;
#if   OPENSSL_VERSION_NUMBER >= 0x10100000L
	EVP_EncryptUpdate(b->ctx, (unsigned char*)out, &encl, b->zeros, 8);
#else
	//	EVP_EncryptUpdate(&b->ctx, out, &encl, b->zeros, 8);
#endif
}

static void expand(uint64_t *counter, uint8_t *buf) {
	const uint8_t *blocks[1] = { buf };
	uint8_t *cur = buf + 32;
	uint8_t hashmash[40];
	int i;
	for (i = 1; i < 4096; i++) {
		SHA256_CTX ctx;
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], counter, 8);
		memcpy(&hashmash[8], blocks[0], 32);
		SHA256_Update(&ctx, hashmash, 40);
		SHA256_Final(cur, &ctx);
		*counter += 1;
		blocks[0] += 32;
		cur += 32;
	}
}

uint8_t* block_index(const struct hash_state *s, size_t i) {
	return s->buffer + (32 * i);
}

void hash_state_init(struct hash_state *s, const struct balloon_options *opts, const uint8_t salt[32]) {
	s->counter = 0;
	s->buffer = (uint8_t*)malloc(131072);
	s->opts = opts;
	bitstream_init(&s->bstream);
	bitstream_seed_add(&s->bstream, salt, 32);
	bitstream_seed_add(&s->bstream, &opts->s_cost, 8);
	bitstream_seed_add(&s->bstream, &opts->t_cost, 4);
	bitstream_seed_finalize(&s->bstream);
}

void hash_state_free(struct hash_state *s) {
	bitstream_free(&s->bstream);
	free(s->buffer);
}

void hash_state_fill(struct hash_state *s, const uint8_t *in, size_t inlen) {
	uint8_t hashmash[132];
	SHA256_CTX c;
	SHA256_Init(&c);
	memcpy(&hashmash[0], &s->counter, 8);
	memcpy(&hashmash[8], in, 32);
	memcpy(&hashmash[40], in, 80);
	memcpy(&hashmash[120], &s->opts->s_cost, 8);
	memcpy(&hashmash[128], &s->opts->t_cost, 4);
	SHA256_Update(&c, hashmash, 132);
	SHA256_Final(s->buffer, &c);
	s->counter++;
	expand(&s->counter, s->buffer);
}
void hash_state_fill_cuda(struct hash_state *s, const uint8_t *in, size_t inlen) {
	uint8_t hashmash[132];
	SHA256_CTX c;
	SHA256_Init(&c);
	memcpy(&hashmash[0], &s->counter, 8);
	memcpy(&hashmash[8], in, 32);
	memcpy(&hashmash[40], in, 80);
	memcpy(&hashmash[120], &s->opts->s_cost, 8);
	memcpy(&hashmash[128], &s->opts->t_cost, 4);
	SHA256_Update(&c, hashmash, 132);
	SHA256_Final(s->buffer, &c);
	s->counter++;
	expand(&s->counter, s->buffer);
}
void hash_state_mix(struct hash_state *s) {
	SHA256_CTX ctx;
	uint8_t buf[8];
	uint8_t hashmash[168];
	int i;

	// round = 0
	uint64_t neighbor;
	for (i = 0; i < 4096; i++) {
		uint8_t *cur_block = s->buffer + (32 * i);
		const uint8_t *blocks[5];
		const uint8_t *prev_block = i ? cur_block - 32 : block_index(s, 4095);
		blocks[0] = prev_block;
		blocks[1] = cur_block;
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[2] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[3] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[4] = block_index(s, neighbor % 4096);
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], &s->counter, 8);
		for (int j = 0; j<5; j++)
			memcpy(&hashmash[8 + (j * 32)], blocks[j], 32);
		SHA256_Update(&ctx, hashmash, 168);
		SHA256_Final(cur_block, &ctx);
		s->counter += 1;
	}
	// round = 1
	for (i = 0; i < 4096; i++) {
		uint8_t *cur_block = s->buffer + (32 * i);
		const uint8_t *blocks[5];
		const uint8_t *prev_block = i ? cur_block - 32 : block_index(s, 4095);
		blocks[0] = prev_block;
		blocks[1] = cur_block;
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[2] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[3] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[4] = block_index(s, neighbor % 4096);
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], &s->counter, 8);
		for (int j = 0; j<5; j++)
			memcpy(&hashmash[8 + (j * 32)], blocks[j], 32);
		SHA256_Update(&ctx, hashmash, 168);
		SHA256_Final(cur_block, &ctx);
		s->counter += 1;
	}
	// round = 2
	for (i = 0; i < 4096; i++) {
		uint8_t *cur_block = s->buffer + (32 * i);
		const uint8_t *blocks[5];
		const uint8_t *prev_block = i ? cur_block - 32 : block_index(s, 4095);
		blocks[0] = prev_block;
		blocks[1] = cur_block;
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[2] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[3] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[4] = block_index(s, neighbor % 4096);
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], &s->counter, 8);
		for (int j = 0; j<5; j++)
			memcpy(&hashmash[8 + (j * 32)], blocks[j], 32);
		SHA256_Update(&ctx, hashmash, 168);
		SHA256_Final(cur_block, &ctx);
		s->counter += 1;
	}
	// round = 3
	for (i = 0; i < 4096; i++) {
		uint8_t *cur_block = s->buffer + (32 * i);
		const uint8_t *blocks[5];
		const uint8_t *prev_block = i ? cur_block - 32 : block_index(s, 4095);
		blocks[0] = prev_block;
		blocks[1] = cur_block;
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[2] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[3] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[4] = block_index(s, neighbor % 4096);
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], &s->counter, 8);
		for (int j = 0; j<5; j++)
			memcpy(&hashmash[8 + (j * 32)], blocks[j], 32);
		SHA256_Update(&ctx, hashmash, 168);
		SHA256_Final(cur_block, &ctx);
		s->counter += 1;
	}
}
void hash_state_mix_cuda(struct hash_state *s) {
	SHA256_CTX ctx;
	uint8_t buf[8];
	uint8_t hashmash[168];
	int i;

	// round = 0
	uint64_t neighbor;
	for (i = 0; i < 4096; i++) {
		uint8_t *cur_block = s->buffer + (32 * i);
		const uint8_t *blocks[5];
		const uint8_t *prev_block = i ? cur_block - 32 : block_index(s, 4095);
		blocks[0] = prev_block;
		blocks[1] = cur_block;
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[2] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[3] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[4] = block_index(s, neighbor % 4096);
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], &s->counter, 8);
		for (int j = 0; j<5; j++)
			memcpy(&hashmash[8 + (j * 32)], blocks[j], 32);
		SHA256_Update(&ctx, hashmash, 168);
		SHA256_Final(cur_block, &ctx);
		s->counter += 1;
	}
	// round = 1
	for (i = 0; i < 4096; i++) {
		uint8_t *cur_block = s->buffer + (32 * i);
		const uint8_t *blocks[5];
		const uint8_t *prev_block = i ? cur_block - 32 : block_index(s, 4095);
		blocks[0] = prev_block;
		blocks[1] = cur_block;
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[2] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[3] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[4] = block_index(s, neighbor % 4096);
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], &s->counter, 8);
		for (int j = 0; j<5; j++)
			memcpy(&hashmash[8 + (j * 32)], blocks[j], 32);
		SHA256_Update(&ctx, hashmash, 168);
		SHA256_Final(cur_block, &ctx);
		s->counter += 1;
	}
	// round = 2
	for (i = 0; i < 4096; i++) {
		uint8_t *cur_block = s->buffer + (32 * i);
		const uint8_t *blocks[5];
		const uint8_t *prev_block = i ? cur_block - 32 : block_index(s, 4095);
		blocks[0] = prev_block;
		blocks[1] = cur_block;
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[2] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[3] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[4] = block_index(s, neighbor % 4096);
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], &s->counter, 8);
		for (int j = 0; j<5; j++)
			memcpy(&hashmash[8 + (j * 32)], blocks[j], 32);
		SHA256_Update(&ctx, hashmash, 168);
		SHA256_Final(cur_block, &ctx);
		s->counter += 1;
	}
	// round = 3
	for (i = 0; i < 4096; i++) {
		uint8_t *cur_block = s->buffer + (32 * i);
		const uint8_t *blocks[5];
		const uint8_t *prev_block = i ? cur_block - 32 : block_index(s, 4095);
		blocks[0] = prev_block;
		blocks[1] = cur_block;
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[2] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[3] = block_index(s, neighbor % 4096);
		bitstream_fill_buffer(&s->bstream, buf, 8);
		neighbor = 0;
		neighbor |= buf[7]; neighbor <<= 8; neighbor |= buf[6]; neighbor <<= 8;
		neighbor |= buf[5]; neighbor <<= 8; neighbor |= buf[4]; neighbor <<= 8;
		neighbor |= buf[3]; neighbor <<= 8; neighbor |= buf[2]; neighbor <<= 8;
		neighbor |= buf[1]; neighbor <<= 8; neighbor |= buf[0];
		blocks[4] = block_index(s, neighbor % 4096);
		SHA256_Init(&ctx);
		memcpy(&hashmash[0], &s->counter, 8);
		for (int j = 0; j<5; j++)
			memcpy(&hashmash[8 + (j * 32)], blocks[j], 32);
		SHA256_Update(&ctx, hashmash, 168);
		SHA256_Final(cur_block, &ctx);
		s->counter += 1;
	}
}
__host__
void balloon256_cpu_init(int thr_id, uint32_t threads)
{
	CUDA_SAFE_CALL(cudaMalloc(&d_KNonce[thr_id], 2 * sizeof(uint32_t)));
}

__host__
void balloon_setBlock_80(int thr_id, void *pdata, const void *pTargetIn)
{
	unsigned char PaddedMessage[80];
	memcpy(PaddedMessage, pdata, 80);
	CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(pTarget, pTargetIn, 8 * sizeof(uint32_t), 0, cudaMemcpyHostToDevice, gpustream[thr_id]));
	CUDA_SAFE_CALL(cudaMemcpyToSymbolAsync(c_PaddedMessage80, PaddedMessage, 10 * sizeof(uint64_t), 0, cudaMemcpyHostToDevice, gpustream[thr_id]));
	if (opt_debug)
		CUDA_SAFE_CALL(cudaDeviceSynchronize());
}
//#ifdef __cplusplus
//}
//#endif