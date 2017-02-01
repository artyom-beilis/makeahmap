#pragma once
#include <CL/cl.h>
#include <stdexcept>
#include <ostream>
#include <iostream>
#include <vector>
#include <sstream>
#include <chrono>

//#define PRINT_DEVICE
//#define ENABLE_PROF

#ifdef ENABLE_PROF
class profiler {
public:
	profiler(std::string const &name=std::string()) : name_(name), total_time_(), total_calls_() {}
	void set_name(std::string const &name) { name_ = name; }
	void start()
	{
		start_ = std::chrono::high_resolution_clock::now();
	}
	void stop()
	{
		auto end = std::chrono::high_resolution_clock::now();
	    	total_time_ += std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end-start_).count();
		total_calls_ ++;
	}
	void add_time(long long ns)
	{
		total_time_ += ns * 1e-9;
		total_calls_++;
	}
	double specific() const
	{
		return total_time_ / total_calls_;
	}
	double total() const
	{
		return total_time_;
	}
	void reset()
	{
		total_calls_=0;
		total_time_ = 0;
	}
	~profiler()
	{
		if(total_calls_ > 0)
			fprintf(stderr,"Kernel %20s total:%8.3f ms specific:%8.3f us calls:%d\n",name_.c_str(),total_time_* 1e3,total_time_ * 1e6 / total_calls_,total_calls_);
	}
private:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_;
	std::string name_;
	double total_time_;
	int total_calls_;
};
#else
class profiler {
public:
	profiler() {}
	profiler(std::string const &) {}
	void set_name(std::string const &) {}
	void start(){}
	void stop(){}
};
#endif

class cl_error : public std::runtime_error {
public:
	cl_error(std::string const &msg) : std::runtime_error("OpenCL: " + msg) {}
	cl_error(std::string const &msg,int code) : std::runtime_error("OpenCL: " + msg + ": " + std::to_string(code)) {}
};


class context_with_program {
public:
	context_with_program() :
		context_(0),
		queue_(0),
		program_(0),
		device_id_()
	{
	}
	~context_with_program()
	{
		clFinish(queue_);
		clReleaseProgram(program_);
		clReleaseCommandQueue(queue_);
		clReleaseContext(context_);
	}
	bool is_cpu()
	{
		cl_device_type dt;
		clGetDeviceInfo(device_id_,CL_DEVICE_TYPE,sizeof(dt),&dt,0);
		return dt == CL_DEVICE_TYPE_CPU;
	}
	std::string name()
	{
		std::string type;
		cl_device_type dt;
		clGetDeviceInfo(device_id_,CL_DEVICE_TYPE,sizeof(dt),&dt,0);
		if(dt == CL_DEVICE_TYPE_GPU) 
			type = "GPU";
		else if(dt == CL_DEVICE_TYPE_CPU) 
			type = "CPU";
		else
			type = "Unknown type";

		char n[256]={};
		clGetDeviceInfo(device_id_,CL_DEVICE_NAME,256,n,0);
		unsigned int mcu=0;
		clGetDeviceInfo(device_id_,CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(mcu),&mcu,0);
		std::ostringstream ss;
		ss << type << ": " <<n << "; computing units=" << mcu;
		return ss.str();
	}
	void load(char const *prg)
	{
		int err = 0;
		cl_platform_id pid=0;
		if((err = clGetPlatformIDs(1,&pid,NULL))!=CL_SUCCESS)
			throw cl_error("Failed to get platform ID",err);
		cl_device_id avail[4];
		unsigned avail_no = 0;
		if((err = clGetDeviceIDs(pid, CL_DEVICE_TYPE_ALL, 4, avail, &avail_no))!=CL_SUCCESS)
			throw cl_error("Failed to get device ID",err);
		if(avail_no == 0)
			throw cl_error("Failed No devices");
		device_id_ = avail[0];
		for(size_t i=0;i<avail_no;i++) {
#ifdef PRINT_DEVICE
			char name[256]={};
			clGetDeviceInfo(avail[i],CL_DEVICE_NAME,256,name,0);
			unsigned int mcu=0;
			clGetDeviceInfo(avail[i],CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(mcu),&mcu,0);
			std::cout << "-- Device supporting OCL found: " << name << " compute units=" << mcu;
#endif			
			cl_device_type dt;
			clGetDeviceInfo(avail[i],CL_DEVICE_TYPE,sizeof(dt),&dt,0);
			if(dt == CL_DEVICE_TYPE_GPU) {
#ifdef PRINT_DEVICE
				std::cout << "; device is GPU, using it";
#endif				
				device_id_ = avail[i];
			}
		}

		if(!(context_ = clCreateContext(0,1,&device_id_,NULL,NULL,&err)))
			throw cl_error("Failed to create context",err);
		int flags = 0;
#ifdef ENABLE_PROF
		flags |= CL_QUEUE_PROFILING_ENABLE;
#endif		
		if(!(queue_ = clCreateCommandQueue(context_,device_id_,flags,&err)))
			throw cl_error("Failed to create queue",err);	
		if(!(program_ = clCreateProgramWithSource(context_,1,&prg,NULL,&err)))
			throw cl_error("Failed to create program",err);	
		err = clBuildProgram(program_, 0, NULL, NULL, NULL, NULL);
		if(err != CL_SUCCESS) {
			size_t len=0;
			char buffer[16384];
			clGetProgramBuildInfo(program_, device_id_, CL_PROGRAM_BUILD_LOG, sizeof(buffer)-1, buffer, &len);
			buffer[len]=0;
			throw cl_error("Failed to build program:" + std::string(buffer),err);
		}
	}
	cl_device_id device_id() { return device_id_; }
	cl_context context() { return context_; }
	cl_command_queue queue() { return queue_; }
	cl_program program() { return program_; }
	void finish()
	{
		int err = clFinish(queue_);
		if(err != CL_SUCCESS)
			throw cl_error("Failed to finish clQueue",err);
	}

private:
	cl_context context_;
	cl_command_queue queue_;
	cl_program program_;
	cl_device_id device_id_;
};

template<typename Item>
class memory_object {
	memory_object(memory_object const &) = delete;
	void operator=(memory_object const &) = delete;
public:
	memory_object(context_with_program &ctx,size_t len,int flags=CL_MEM_READ_WRITE) :
		ctx_(ctx)
	{
		int err = 0;
		data_ = clCreateBuffer(ctx.context(),flags,sizeof(Item)*len,NULL,&err);
		if(!data_)
			throw cl_error("Failed to allocate device memory",err);
		size_ = len;
	}
	memory_object(context_with_program &ctx,Item const *ptr,size_t len,int flags=CL_MEM_READ_WRITE) :
		ctx_(ctx),
		data_()
	{
		int err = 0;
		data_ = clCreateBuffer(ctx.context(),flags | CL_MEM_COPY_HOST_PTR,sizeof(Item)*len,const_cast<void *>(static_cast<void const *>(ptr)),&err);
		if(!data_)
			throw cl_error("Failed to allocate device memory",err);
		size_ = len;
	}
	cl_mem &data() { return data_; }
	void copy_betweem(memory_object &out,int len)
	{
		int err = clEnqueueCopyBuffer(ctx_.queue(),data_,out.data_,0,0,sizeof(Item) * len,0,NULL,NULL);
		if(err != CL_SUCCESS)
			throw cl_error("Failed to copy memory",err);
	
	}
	void copy_from_host(Item const *data,int len)
	{
		int err = clEnqueueWriteBuffer(ctx_.queue(),data_,CL_TRUE,0,sizeof(Item) * len, data,0,NULL,NULL);
		if(err != CL_SUCCESS)
			throw cl_error("Failed to copy to device",err);
	}
	void copy_to_host(Item *data,int len)
	{
		int err = clEnqueueReadBuffer(ctx_.queue(),data_,CL_TRUE,0,sizeof(Item) * len, data,0,NULL,NULL);
		if(err != CL_SUCCESS)
			throw cl_error("Failed to copy to device",err);
	}
	void free()
	{
		if(data_)
			clReleaseMemObject(data_);
		data_=nullptr;
	}
	void print(std::ostream &out)
	{
		std::vector<Item> tmp(size_);
		copy_to_host(tmp.data(),size_);
		out << '[';
		for(auto ptr = tmp.begin();ptr!=tmp.end();++ptr) {
			if(ptr != tmp.begin())
				out << ',';
			out << *ptr;
		}
		out << "]" << std::endl;
	}
	~memory_object() 
	{
		free();
	}
private:
	friend class kernel;
	context_with_program &ctx_;
	cl_mem data_;
	size_t size_;
};

template<typename T>
std::ostream &operator<<(std::ostream &out,memory_object<T>  &o)
{
	o.print(out);
	return out;
}

struct local_proxy { 
	size_t size;
};

inline local_proxy local_memory(size_t len)
{
	local_proxy lp={len};
	return lp;
}

class kernel {
	kernel(kernel const &) = delete;
	void operator=(kernel const &) = delete;
public:

	kernel(context_with_program &ctx) :
		ctx_(ctx), kernel_(),local_(),param_(),lwg_(0)
	{
	}
	kernel(context_with_program &ctx,std::string const &name) : 
		ctx_(ctx),
		kernel_(),
		local_(),
		param_(),
		name_(name),
		prof_(name)
	{
		int err=0;
		kernel_ = clCreateKernel(ctx.program(), name.c_str(), &err);
		if (!kernel_ || err != CL_SUCCESS)
		{
			if(kernel_)
				clReleaseKernel(kernel_);
			throw cl_error("Failed to create kernel " + name,err);
		}
		name_ = name;
	}
	void assign(std::string const &name)
	{
		prof_.set_name(name);
		name_ = name;
		if(kernel_)
			clReleaseKernel(kernel_);
		int err =0;
		kernel_ = clCreateKernel(ctx_.program(), name.c_str(), &err);
		if (!kernel_ || err != CL_SUCCESS)
		{
			if(kernel_) {
				clReleaseKernel(kernel_);
				kernel_=0;
			}
			throw cl_error("Failed to create kernel " + name,err);
		}
	}
	void local(size_t l) { local_ = l; }
	void reset() { param_=0; }
	void bind(local_proxy const &lp)
	{
		int err = clSetKernelArg(kernel_,param_ ++,lp.size,0);
		if(err!=CL_SUCCESS)
			throw cl_error("Failed to set local kernel argument of size " + std::to_string(lp.size) ,err);
		
	}
	template<typename Item>
	void bind(memory_object<Item> const &mem)
	{
		int err = clSetKernelArg(kernel_,param_ ++,sizeof(cl_mem),&mem.data_);
		if(err != CL_SUCCESS)
			throw cl_error("Failed to set kernel argument #" + std::to_string(param_-1) + "items=" + std::to_string(mem.size_) + " elememt-size=" +std::to_string(sizeof(Item)) + " to " + name_,err);
			
	}
	size_t get_local_workgroup()
	{
		if(lwg_!=0)
			return lwg_;
		size_t local=0;
		int err = clGetKernelWorkGroupInfo(kernel_, ctx_.device_id(), CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
		if(err!=CL_SUCCESS)
			throw cl_error("Failed to get kernel WG",err);
		lwg_ = local;
		return local;
	}
	void bind_all() {}
	template<typename T1,typename ...Targs>
	void bind_all(T1 const &p,Targs const &...args)
	{
		bind(p);
		bind_all(args...);
	}
	template<typename ...T>
	void operator()(int items,T const &...args)
	{
		reset();
		bind_all(items,args...);
		enqueue(items);
	}
	void enqueue(size_t d1,size_t d2) {
		size_t sizes[]={d1,d2};
		int err = clEnqueueNDRangeKernel(ctx_.queue(),kernel_, 2, NULL, sizes, 0 , 0, NULL, 0);
		if(err!=CL_SUCCESS)
			throw cl_error("Failed to enqueue kernel",err);
	}
	void enqueue(size_t global_items)
	{
		if(local_ != 0) 
			global_items = (global_items + local_ - 1)/local_*local_;
#ifdef ENABLE_PROF 	
		cl_event ev;
		cl_event *ev_ptr = &ev;
#else
		cl_event *ev_ptr = 0;
#endif
		int err = clEnqueueNDRangeKernel(ctx_.queue(),kernel_, 1, NULL, &global_items, (local_ ? &local_ : 0 ), 0, NULL, ev_ptr);
		if(err!=CL_SUCCESS)
			throw cl_error("Failed to enqueue kernel",err);
#ifdef ENABLE_PROF 	
		clWaitForEvents(1,&ev);
		cl_ulong start=0,stop=0;
		clGetEventProfilingInfo(ev,CL_PROFILING_COMMAND_START,sizeof(start),&start,0);
		clGetEventProfilingInfo(ev,CL_PROFILING_COMMAND_END,sizeof(stop),&stop,0);
		prof_.add_time(stop-start);
#endif		
	}
	template<typename Item>
	void bind(Item &it)
	{
		int err = clSetKernelArg(kernel_,param_ ++,sizeof(it),&it);
		if(err != CL_SUCCESS)
			throw cl_error("Failed to set kernel argument #" + std::to_string(param_-1),err);

	}
	~kernel()
	{
		if(kernel_)
			clReleaseKernel(kernel_);
	}
	profiler &prof()
	{
		return prof_;
	}
private:
	context_with_program &ctx_;
	cl_kernel kernel_;
	size_t local_;
	int param_;
	std::string name_;
	size_t lwg_;
	profiler prof_;
};


