#include "surface.h"
#ifdef USE_CPU
#ifdef USE_OMP
#include <omp.h>
#endif
#include "solver.h"
#else
#include "solver_ocl.h"
#endif
#include <iostream>
#include <sstream>
#include <chrono>

struct surface_solver::data {
#ifndef USE_CPU	
	eq_solver slv;
#endif	
};

surface_solver::surface_solver() : impl(new data()) {}
surface_solver::~surface_solver() {}

#ifdef USE_CPU
std::string surface_solver::name()
{
	return solver_name();
}
#else
std::string surface_solver::name()
{
	return "OpenCL: " + impl->slv.context().name();
}
#endif

std::pair<int,double> surface_solver::run(std::vector<std::vector<char> > const &bmask,std::vector<std::vector<float> > &bvalues,float thresh)
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
#ifdef USE_CPU
	sparce_matrix X;
	X.reserve(variables);
	sparce_matrix &Matrix = X;
#else
	eq_solver &Matrix = impl->slv;
	Matrix.init_matrix(variables);
#endif	
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
				Matrix.add(i,rindex[r2][c2]);
		}
	}

	auto start = std::chrono::high_resolution_clock::now();
#ifdef USE_CPU	
	int iterations = solve(variables,X,y.data(),x0.data(),thresh,variables);
#else	
	int iterations = Matrix.solve(y.data(),x0.data(),thresh,variables);
#endif
	auto end = std::chrono::high_resolution_clock::now();
	double time = std::chrono::duration_cast<std::chrono::duration<float, std::ratio<1> > >(end-start).count();

	for(int i=0;i<variables;i++) {
		bvalues[index[i].first][index[i].second]=x0[i];
	}
	return std::pair<int,double>(iterations,time);
}


