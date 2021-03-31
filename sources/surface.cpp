#include "surface.h"
#include "solver.h"
#include "surface_solver.h"
#include <iostream>
#include <sstream>
#include <chrono>
#include <thread>

#if defined _WIN32 || defined WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#include <libgen.h>
#include <unistd.h>
#endif

extern "C" {
	typedef surface_solver_base *(*get_gpu_solver_type)(char *,size_t,int);
}

std::unique_ptr<surface_solver_base> load_gpu_solver(int pid)
{
	get_gpu_solver_type get_gpu_solver = nullptr;

#if defined _WIN32 || defined WIN32
	UINT oldMode = SetErrorMode(0);
	SetErrorMode(oldMode | SEM_FAILCRITICALERRORS | SEM_NOOPENFILEERRORBOX);
	HMODULE h = LoadLibrary("libmakeahmap_gpu.dll");
	SetErrorMode(oldMode);
	if(!h)
		throw std::runtime_error("Failed to load libmakeahmap_gpu.dll");
	FARPROC sym = GetProcAddress(h,"get_gpu_solver");
	if(!sym)
		throw std::runtime_error("Failed to resolve get_gpu_solver");
#else
	char path[4096]={};
	int size = readlink("/proc/self/exe",path,sizeof(path)-256);
	if(size < 0)
		std::runtime_error("Failed to read /proc/self/exe");
	char *p=dirname(path);
	strcat(p,"/libmakeahmap_gpu.so");
		
	void *h = dlopen(p,RTLD_LAZY | RTLD_GLOBAL);
	if(!h) {
		throw std::runtime_error("Failed to load " + std::string(p));
	}
	void *sym = dlsym(h,"get_gpu_solver");
	if(!sym) {
		throw std::runtime_error("Failed to resolve get_gpu_solver");
	}
#endif	
	get_gpu_solver = reinterpret_cast<get_gpu_solver_type>(sym);

	char buf[4096]="Unknown";
	std::unique_ptr<surface_solver_base> r(get_gpu_solver(buf,sizeof(buf),pid));
	if(!r.get()) {
		throw std::runtime_error(buf);
	}
	return std::move(r);
}


std::unique_ptr<surface_solver_base> get_solver(surface_solver_options const &opt)
{
	std::unique_ptr<surface_solver_base> ptr;
	
	if(!opt.force_cpu) {
		try {
			ptr = std::move(load_gpu_solver(opt.platform_id));
			if(!(ptr->is_cpu())) 
				return std::move(ptr);
			if(opt.allow_cpu) {
				std::cout <<"   WARNING: Only CPU OpenCL support avalible (suboptimal), using it as per user request " << std::endl;
				return std::move(ptr);
			}
			std::cout <<"   NOTE: No GPU OpenCL support avalible, falling back to CPU" << std::endl;
		}
		catch(std::exception const &e) {
			std::cout <<"   WARMING: Failed to create OpenCL solver, falling back to CPU: " << e.what() << std::endl;
		}
	}
	
	ptr.reset(new surface_solver<cpu::eq_solver>());
	return std::move(ptr);

}

