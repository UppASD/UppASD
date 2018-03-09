#ifndef __GPU_ERROR_CHK_HPP__
#define __GPU_ERROR_CHK_HPP__

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);}

#include <cuda.h>
#include <stdio.h>

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}
#endif

	//hostMatrix<unsigned int,2> loca;
	//loca.initiate(ex.neighbourCount.dimension_size(0),ex.neighbourCount.dimension_size(1));
	//ex.neighbourCount.write(loca);
	//cout<<"Here: "<<loca(0,0)<<" "<<loca(1,0)<<endl;
	//exit(0);
