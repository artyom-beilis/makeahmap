#include "surface.h"
#include "solver.h"

#ifdef ENABLE_OCL
#include "solver_ocl.h"
#endif

#include <iostream>
#include <sstream>
#include <chrono>


template<typename Solver>
class surface_solver : public surface_solver_base {
public:
	surface_solver() {}
	virtual ~surface_solver() {}
	virtual std::string name()
	{
		return slv.name();
	}
	stats run(std::vector<std::vector<char> > const &bmask,std::vector<std::vector<float> > &bvalues,float thresh)
	{
		int N=bmask.size();
		std::vector<std::pair<int,int> > index;
		std::vector<std::vector<int> > rindex(N,std::vector<int>(N,0));
		index.reserve(N*N);
		int variables=0;
		for(int i=0;i<N;i++) {
			for(int j=0;j<N;j++) {
				if(!bmask[i][j]) {
					index.push_back(std::make_pair(i,j));
					rindex[i][j]=variables;
					variables++;
				}
				else
					rindex[i][j]=-1;
			}
		}

		slv.init_matrix(variables);
		std::vector<float> y(variables,0);
		std::vector<float> x0(variables,0);

		for(int i=0;i<variables;i++) {
			int r=index[i].first;
			int c=index[i].second;
			x0[i]=bvalues[r][c];
			int drs[4]={0,0,-1,+1};
			int dcs[4]={1,-1,0,0};
			for(int p=0;p<4;p++) {
				int r2=r+drs[p];
				int c2=c+dcs[p];
				if(bmask[r2][c2])
					y[i]+= bvalues[r2][c2] * 0.25;
				else 
					slv.add(i,rindex[r2][c2]);
			}
		}

		auto start = std::chrono::high_resolution_clock::now();
		std::pair<int,double> st = slv.solve(y.data(),x0.data(),thresh,variables);
		auto end = std::chrono::high_resolution_clock::now();

		double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end-start).count();

		for(int i=0;i<variables;i++) {
			bvalues[index[i].first][index[i].second]=x0[i];
		}

		stats res;
		res.iterations = st.first;
		res.bandwidth = double(slv.get_bytes_per_it()) * res.iterations * variables /  st.second / (1024*1024*1024);
		res.time = time;
		return res;
	}
	Solver slv;
};

std::unique_ptr<surface_solver_base> get_solver(surface_solver_options const &opt)
{
	std::unique_ptr<surface_solver_base> ptr;
#ifdef ENABLE_OCL
	if(!opt.force_cpu) {
		try {
			std::unique_ptr<surface_solver<gpu::eq_solver> > rptr (new surface_solver<gpu::eq_solver>());
			if(rptr->slv.is_cpu()){
				if(opt.allow_cpu) {
					std::cout <<"   WARNING: Only CPU OpenCL support avalible (suboptimal), using it as per user request " << std::endl;
					ptr = std::move(rptr);
				}
				else {
					std::cout <<"   NOTE: No GPU OpenCL support avalible, falling back to CPU" << std::endl;
				}
			}
		}
		catch(std::exception const &e) {
			std::cout <<"   WARMING: Failed to create OpenCL solver, falling back to CPU: " << e.what() << std::endl;
		}
		if(ptr)
			return std::move(ptr);
	}
#endif
	ptr.reset(new surface_solver<cpu::eq_solver>());
	return std::move(ptr);

}

