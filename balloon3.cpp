/*
* balloon algorithm
*
*/
#include "miner.h"
#include <string.h>
#include <stdint.h>

#include <openssl/sha.h>

#include "balloon.h"
#include "cuda_helper.h"
extern void balloon256_cpu_init(int thr_id, uint32_t threads);
extern void balloon_setBlock_80(int thr_id, void *pdata, const void *ptarget);
//int scanhash_balloon(int thr_id, struct work *work, uint32_t max_nonce, uint32_t *hashes_done)

int scanhash_balloon(int thr_id, uint32_t *pdata,
	uint32_t *ptarget, uint32_t max_nonce,
	uint32_t *hashes_done)
{
	static THREAD uint32_t *h_nounce = nullptr;  //monk1

	uint32_t _ALIGN(128) hash32[8];
	uint32_t _ALIGN(128) endiandata[20];
	//uint32_t *pdata = work->data;
	//uint32_t *ptarget = work->target;

	const uint32_t Htarg = ptarget[7];
	const uint32_t first_nonce = pdata[19];

	uint32_t intensity = (device_sm[device_map[thr_id]] > 500) ? 1 << 28 : 1 << 27;;
	uint32_t throughputmax = device_intensity(device_map[thr_id], __func__, intensity); // 256*4096
	uint32_t throughput = min(throughputmax, max_nonce - first_nonce) & 0xfffffc00;

	uint32_t n = first_nonce;
	static THREAD volatile bool init = false;
	if (!init)
	{
		if (throughputmax == intensity)
			applog(LOG_INFO, "GPU #%d: using default intensity %.3f", device_map[thr_id], throughput2intensity(throughputmax));
		CUDA_SAFE_CALL(cudaSetDevice(device_map[thr_id]));
		CUDA_SAFE_CALL(cudaDeviceReset());
		CUDA_SAFE_CALL(cudaSetDeviceFlags(cudaschedule));
		CUDA_SAFE_CALL(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
		CUDA_SAFE_CALL(cudaStreamCreate(&gpustream[thr_id]));
		CUDA_SAFE_CALL(cudaMallocHost(&h_nounce, 2 * sizeof(uint32_t)));
		
		balloon256_cpu_init(thr_id, (int)throughputmax); 

		mining_has_stopped[thr_id] = false;
		init = true;
	}
	for (int i = 0; i < 19; i++) {
		be32enc(&endiandata[i], pdata[i]);
	};
	balloon_setBlock_80(thr_id, (void*)endiandata, ptarget);
	do {
		be32enc(&endiandata[19], n);
		//ballon_128_cuda(thr_id, (int)throughput, (unsigned char *)endiandata, (unsigned char *)hash32, h_nounce);
		balloon(thr_id, (int)throughput, (unsigned char *)endiandata, (unsigned char *)hash32, h_nounce);
		if (stop_mining) { mining_has_stopped[thr_id] = true; cudaStreamDestroy(gpustream[thr_id]); pthread_exit(nullptr); }

		if (hash32[7] < Htarg && fulltest(hash32, ptarget)) 
		{
			//work_set_target_ratio(work, hash32);
			*hashes_done = n - first_nonce + throughput;

			pdata[19] = n;
			return true;
		}
		n++;
		pdata[19] += throughput; CUDA_SAFE_CALL(cudaGetLastError());
	} while (!work_restart[thr_id].restart && ((uint64_t)max_nonce > ((uint64_t)(pdata[19]) + (uint64_t)throughput)));

	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;

	return 0;
}