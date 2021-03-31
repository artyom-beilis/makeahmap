#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <chrono>
#include <memory>

#ifdef USE_CPU

#include "solver.h"

#elif defined USE_VIE

#include "vie.h"

#else
#include "solver_ocl.h"
#endif

void solve_surface(std::vector<std::vector<char> > const &bmask,std::vector<std::vector<float> > &bvalues,float thresh,int platform_id)
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
#elif defined USE_VIE
	sp Matrix;
	Matrix.init(variables);
#else
	static std::unique_ptr<eq_solver> slv_ptr;
	if(!slv_ptr)
		slv_ptr.reset(new eq_solver(platform_id));	
	eq_solver &Matrix = *slv_ptr;
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
	solve(variables,X,y.data(),x0.data(),thresh,variables);
#elif defined USE_VIE
	solve(variables,Matrix,y,x0,thresh,N);	
#else	
	Matrix.solve(y.data(),x0.data(),thresh,variables);
#endif
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::duration<float, std::ratio<1> > >(end-start).count() << std::endl;
		
	for(int i=0;i<variables;i++) {
		bvalues[index[i].first][index[i].second]=x0[i];
	}
}


int main(int argc,char **argv)
{
	auto start = std::chrono::high_resolution_clock::now();
	float thresh = argc>= 4 ? float(atof(argv[3])) : 1e-4f;
	int NS=argc >= 3 ? atoi(argv[2]) : 128;
	int NM=argc >= 2 ? atoi(argv[1]) : 128;
	std::vector<std::vector<float> > prev_values;
	int N;
	for(N=NS;N<=NM;N*=2) {
		std::vector<std::vector<float> > bvalues(N,std::vector<float>(N,0));
		if(!prev_values.empty() && prev_values.size()*2==size_t(N)) {
			for(int r=0;r<N-1;r++) {
				for(int c=0;c<N-1;c++) {
					if((r&1)==0 && (c&1)==0)
						bvalues[r][c]=prev_values[r/2][c/2];
					else if((r&1) == 0) 
						bvalues[r][c]=(prev_values[r/2][c/2] + prev_values[r/2][c/2+1])/2;
					else if((c&1) == 0) 
						bvalues[r][c]=(prev_values[r/2][c/2] + prev_values[r/2+1][c/2])/2;
					else
						bvalues[r][c]=(prev_values[r/2][c/2] + prev_values[r/2][c/2+1] +
						               prev_values[r/2][c/2] + prev_values[r/2+1][c/2])/4;
				}
			}
		}
		std::vector<std::vector<char> >   bmask(N,std::vector<char>(N,0));

		for(int i=0;i<N;i++) {
			bmask[i][0]=1;
			bmask[i][N-1]=1;
			bmask[0][i]=1;
			bmask[N-1][i]=1;
		}
		for(int r=0;r<N;r++) {
			for(int c=0;c<N;c++) {
				float y=float(r)/(N-1)-0.5;
				float x=float(c)/(N-1)-0.5;
				if(x*x+y*y < 0.25*0.25) {
					bmask[r][c]=1;
					bvalues[r][c]=1.0*fabs(x) * 2;
				}
			}
		}
		std::cout << "Solving for " << N << std::endl;
		solve_surface(bmask,bvalues,thresh);
		prev_values.swap(bvalues);
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Total Time: "<< std::chrono::duration_cast<std::chrono::duration<float, std::ratio<1> > >(end-start).count() << std::endl;

	N=prev_values.size();

	std::vector<char> res(N*N,0);
	int p=0;
	std::ofstream csv("out.csv");
	for(int r=0;r<N;r++) {
		for(int c=0;c<N;c++) {
			res[p++]=std::max(0,std::min(255,int(prev_values[r][c]*255)));
			csv << prev_values[r][c];
			if(c+1 < N)
				csv<<',';
		}
		csv << '\n';
	}
	std::ofstream pgm("out.pgm",std::fstream::binary);
	pgm<<"P5 " << N << " " << N << " 255\n";
	pgm.write(res.data(),N*N);
}
