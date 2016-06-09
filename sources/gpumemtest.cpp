#define ENABLE_PROF
#include <chrono>
#include <vector>
#include "clcpp.h"


int main()
{
	try {
		context_with_program ctx;
		ctx.load(R"xxx(
		__kernel void mem_copy(int N,__global float4 *tgt,__global float4 const *src)
		{
			int id = get_global_id(0);
			if(id < N)
				tgt[id]=src[id];
		}
		)xxx");
		std::cout << "Device:" << ctx.name() << std::endl;
		kernel mem_copy(ctx,"mem_copy");
		size_t total=1024ULL*1024*1024/4;
		int N=total/sizeof(cl_float4);
		cl_float4 fl={1.0,2.0,3.0,4.0};
		std::vector<cl_float4> sample(N,fl);
		memory_object<cl_float4> src(ctx,sample.data(),N);
		memory_object<cl_float4> tgt(ctx,N);
		mem_copy.local(128);
		int sampling=100;
		for(int i=0;i<sampling;i++) {
			mem_copy(N,src,tgt);
		}
		double time = mem_copy.prof().specific();
		mem_copy.prof().reset();
		double throughput = 2 * total / time / (1024Ull*1024*1024);
		std::cout << "Copy by Kernel GB/s:" <<throughput << std::endl;
		{
			profiler prf("");
			for(int i=0;i<sampling;i++) {
				prf.start();
				cl_event ev;
				int err=0;
				if(err != clEnqueueCopyBuffer(ctx.queue(),src.data(),tgt.data(),0,0,total,0,NULL,&ev))
					throw cl_error("Failed to copy memory",err);
				clWaitForEvents(1,&ev);
				prf.stop();
			}
			double time = prf.specific();
			prf.reset();
			double throughput = 2 * total / time / (1024Ull*1024*1024);
			std::cout << "Copy using device GB/s:" <<throughput << std::endl;
		}
	}
	catch(std::exception const &e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}
	return 0;
}
