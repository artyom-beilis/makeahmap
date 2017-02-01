#include "solver_ocl.h"
#include "surface_solver.h"

extern "C" surface_solver_base *get_gpu_solver(char *msg,size_t len)
{
	try {
		surface_solver_base *ptr = new surface_solver<gpu::eq_solver>();
		return ptr;
	}
	catch(std::exception const &e) {
		snprintf(msg,len,"%s",e.what());
		return 0;
	}
	catch(...) {
		snprintf(msg,len,"Unknown error");
		return 0;
	}
}
