#pragma once
#include <string>
#include <vector>
#include <chrono>
#include <thread>
#include "surface.h"

template<typename Solver>
class surface_solver : public surface_solver_base {
public:
	template<typename ...T>
	surface_solver(T... args) : slv(args...) {}
	virtual ~surface_solver() {}
	virtual std::string name()
	{
		return slv.name();
	}
	virtual bool is_cpu()
	{
		return slv.is_cpu();
	}
	stats run(std::vector<std::vector<char> > const &bmask,std::vector<std::vector<float> > &bvalues,float thresh)
	{
		int N=bmask.size();
		std::vector<std::pair<int,int> > index;
		std::vector<std::vector<int> > rindex(N,std::vector<int>(N,0));
		index.reserve(N*N);
		int variables = 0;
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
		
		int align_factor = 4 * std::thread::hardware_concurrency();
		int variables_aligned = (variables + align_factor - 1) / align_factor * align_factor;

		slv.init_matrix(variables_aligned);
		std::vector<float> y(variables_aligned,0);
		std::vector<float> x0(variables_aligned,0);

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
		std::pair<int,double> st = slv.solve(y.data(),x0.data(),thresh,variables_aligned);
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


