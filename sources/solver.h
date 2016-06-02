#pragma once
#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef USE_OMP


class sparce_matrix {
	struct row {
		int   c[4];
		float v[4];
	};
public:
	void reserve(size_t n) 
	{
		rows_.resize(n);
		for(size_t i=0;i<n;i++) {
			for(int j=0;j<4;j++) {
				rows_[i].c[j]=i;
				rows_[i].v[j]=0.0f;
			}
		}
	}
	void add(int r,int c)
	{
		row &ind = rows_[r];
		int i;
		for(i=0;i<4;i++) {
			if(ind.c[i]==r) {
				ind.c[i]=c;
				ind.v[i]=-0.25;
				break;
			}
		}
		assert(i<4);
	}
	void mpl(int N,float const *vin,float *vout) const
	{
		int size = rows_.size();
		#pragma omp for
		for(int i=0;i<size;i++) {
			row const &r = rows_[i];
			float v[4]= { 
				vin[r.c[0]],  
				vin[r.c[1]],  
				vin[r.c[2]],  
				vin[r.c[3]]
			};
			float sum = r.v[0]*v[0]+r.v[1]*v[1]+r.v[2]*v[2]+r.v[3]*v[3];
			vout[i]=sum + vin[i];
		}
	}
private:
	std::vector<row> rows_;
};


void solve(int N,sparce_matrix const &A,float const *b,float *x,float thresh=1e-8,int limit=-1)
{
	std::vector<float> rv(N,0);
	std::vector<float> pv(N,0);
	std::vector<float> Apv(N,0);

	float *r=rv.data();
	float *p=pv.data();
	float *Ap = Apv.data();

	if(limit==-1)
		limit=N;
	
	float th2=thresh*thresh;

	A.mpl(N,x,r);
	float rsold=0;
	#pragma omp parallel for reduction(+:rsold)
	for(int i=0;i<N;i++) {
		p[i] = r[i]=b[i]-r[i];
		rsold+=r[i]*r[i];
	}
	int it;
	float rsnew = 0;	
	
	for(it=0;it<limit;it++) {
		A.mpl(N,p,Ap);
		float pAp=0;
		#pragma omp parallel for reduction(+:pAp)
		for(int i=0;i<N;i++)
			pAp+=p[i]*Ap[i];
		float alpha = rsold / pAp;
		rsnew=0;
		#pragma omp for
		for(int i=0;i<N;i++)  {
			x[i]+=alpha*p[i];
			r[i]-=alpha*Ap[i];
		}
		#pragma omp parallel for reduction(+:rsnew)
		for(int i=0;i<N;i++) {
			rsnew +=r[i]*r[i];
		}
		if(rsnew < th2)
			break;
		float factor = rsnew / rsold;
		#pragma omp for
		for(int i=0;i<N;i++)
			p[i]=r[i]+factor*p[i];
		rsold = rsnew;


	}
	std::cout << "Counted in " << it << " with error=" << sqrt(rsnew) << std::endl;
}

#else

class sparce_matrix {
	struct row {
		int   c[4];
		float v[4];
	};
public:
	void reserve(size_t n) 
	{
		rows_.resize(n);
		for(size_t i=0;i<n;i++) {
			for(int j=0;j<4;j++) {
				rows_[i].c[j]=i;
				rows_[i].v[j]=0.0f;
			}
		}
	}
	void add(int r,int c)
	{
		row &ind = rows_[r];
		int i;
		for(i=0;i<4;i++) {
			if(ind.c[i]==r) {
				ind.c[i]=c;
				ind.v[i]=-0.25;
				break;
			}
		}
		assert(i<4);
	}
	void mpl(int N,float const *vin,float *vout) const
	{
		int size = rows_.size();
		for(int i=0;i<size;i++) {
			row const &r = rows_[i];
			float v[4]= { 
				vin[r.c[0]],  
				vin[r.c[1]],  
				vin[r.c[2]],  
				vin[r.c[3]]
			};
			float sum = r.v[0]*v[0]+r.v[1]*v[1]+r.v[2]*v[2]+r.v[3]*v[3];
			vout[i]=sum + vin[i];
		}
	}
private:
	std::vector<row> rows_;
};



void solve(int N,sparce_matrix const &A,float const *b,float *x,float thresh=1e-8,int limit=-1)
{
	std::vector<float> rv(N,0);
	std::vector<float> pv(N,0);
	std::vector<float> Apv(N,0);

	float *r=rv.data();
	float *p=pv.data();
	float *Ap = Apv.data();

	if(limit==-1)
		limit=N;
	
	float th2=thresh*thresh;

	A.mpl(N,x,r);
	float rsold=0;
	for(int i=0;i<N;i++) {
		p[i] = r[i]=b[i]-r[i];
		rsold+=r[i]*r[i];
	}
	int it;
	float rsnew = 0;	


	for(it=0;it<limit;it++) {
		A.mpl(N,p,Ap);
		float pAp=0;
		for(int i=0;i<N;i++)
			pAp+=p[i]*Ap[i];
		float alpha = rsold / pAp;
		for(int i=0;i<N;i++) 
			x[i]+=alpha*p[i];
		rsnew=0;
		for(int i=0;i<N;i++) {
			r[i]-=alpha*Ap[i];
			rsnew +=r[i]*r[i];
		}
		if(rsnew < th2)
			break;
		float factor = rsnew / rsold;
		for(int i=0;i<N;i++)
			p[i]=r[i]+factor*p[i];
		rsold = rsnew;


	}
	std::cout << "Counted in " << it << " with error=" << sqrt(rsnew) << std::endl;
}

#endif

