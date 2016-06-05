#pragma once
#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#ifdef CPU_BARRIERS
#include <thread>
#include <pthread.h>
#endif


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
	void mpl(int N,float const *vin,float *vout) const;
#ifdef CPU_BARRIERS
	void mpl(int from,int to,float const *vin,float *vout) const;
#endif	

private:
	std::vector<row> rows_;
};


#ifdef USE_SMP

#include <thread>
#include <future>

std::string solver_name()
{
	int threads = std::thread::hardware_concurrency();
	std::ostringstream ss;
	ss << "CPU using up to " << threads << " threads" << std::endl;
	return ss.str();
}

template<typename Func>
void run(Func const &f,int range)
{
	static int threads = std::min<int>(16,std::thread::hardware_concurrency());
	std::future<void> res[16];
	int chunk=(range+threads-1)/threads;
	int start = 0;
	for(int i=0;i<threads;i++) {
		int next = start + chunk;
		if(next >= range)
			next = range;
		res[i]=std::move(std::async(f,start,next));
		start = next;
	}
	for(int i=0;i<threads;i++)
		res[i].get();
}

template<typename Func>
float reduce(Func const &f,int range)
{
	static int threads = std::min<int>(16,std::thread::hardware_concurrency());
	std::future<float> res[16];
	int chunk=(range+threads-1)/threads;
	int start = 0;
	for(int i=0;i<threads;i++) {
		int next = start + chunk;
		if(next >= range)
			next = range;
		res[i]= std::move(std::async(f,start,next));
		start = next;
	}
	float result = 0;
	for(int i=0;i<threads;i++)
		result += res[i].get();
	return result;
}

void sparce_matrix::mpl(int N,float const *vin,float *vout) const
{
	int size = rows_.size();
	run([&](int start,int end) {
		for(int i=start;i<end;i++) {
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
	},size);
}


/// CMP SMP 
int solve(int N,sparce_matrix const &A,float const *b,float *x,float thresh=1e-8,int limit=-1)
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
		float pAp=reduce([&](int start,int end) {
			float res=0;
			for(int i=start;i<end;i++)
				res+=p[i]*Ap[i];
			return res;
		},N);

		float alpha = rsold / pAp;

		run([&](int begin,int end){
			for(int i=begin;i<end;i++) 
				x[i]+=alpha*p[i];
		},N);
		rsnew=reduce([&](int start,int end) {
			float res=0;
			for(int i=start;i<end;i++) {
				r[i]-=alpha*Ap[i];
				res +=r[i]*r[i];
			}
			return res;
		},N);
		if(rsnew < th2)
			break;
		float factor = rsnew / rsold;
		run([&](int start,int end) {
			for(int i=start;i<end;i++)
				p[i]=r[i]+factor*p[i];
		},N);
		rsold = rsnew;

	}
	return it;
}


#elif defined USE_OMP
std::string solver_name()
{
	return "CPU using OpenMP";
}

void sparce_matrix::mpl(int N,float const *vin,float *vout) const
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
		vout[i]=vin[i] + (v[0]*r.v[0]+v[1]*r.v[1]+v[2]*r.v[2]+v[3]*r.v[3]);
	}
}

/// CPU OMP
int solve(int N,sparce_matrix const &A,float const *b,float *x,float thresh=1e-8,int limit=-1)
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
	return it;
}

#elif defined CPU_BARRIERS
std::string solver_name()
{
	return "Multiple CPU with Barriers";
}

void sparce_matrix::mpl(int from,int to,float const *vin,float *vout) const
{
	for(int i=from;i<to;i++) {
		row const &r = rows_[i];
		float v[4]= { 
			vin[r.c[0]],  
			vin[r.c[1]],  
			vin[r.c[2]],  
			vin[r.c[3]]
		};
		vout[i]=vin[i] + (v[0]*r.v[0]+v[1]*r.v[1]+v[2]*r.v[2]+v[3]*r.v[3]);
	}
}

float sum(float const *v,int n)
{
	float r=0;
	for(int i=0;i<n;i++)
		r+=v[i];
	return r;
}

class barrier {
	barrier(barrier const &other) = delete;
	void operator=(barrier const &other)=delete;
public:
	barrier(int n)
	{
		pthread_barrier_init(&b_,0,n);
	}
	void wait()
	{
		pthread_barrier_wait(&b_);
	}
	~barrier() 
	{
		pthread_barrier_destroy(&b_);
	}
private:
	pthread_barrier_t b_;
};

// CPU single thread
int solve(int N,sparce_matrix const &A,float const *b,float *x,float thresh=1e-8,int limit=-1)
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

	A.mpl(0,N,x,r);
	float grsold=0;
	for(int i=0;i<N;i++) {
		p[i] = r[i]=b[i]-r[i];
		grsold+=r[i]*r[i];
	}
	int it;
	
	float pAp_acc[16];
	float rsnew_acc[16];
	int threads = 2;
	barrier bar(threads);

	int iters;
	auto runner = [&](int from,int to,int id) {	
		float rsold = grsold;
		for(it=0;it<limit;it++) {
			A.mpl(from,to,p,Ap);
			bar.wait();

			float mypAp=0;
			for(int i=from;i<to;i++)
				mypAp+=p[i]*Ap[i];
			pAp_acc[id]=mypAp;
			bar.wait();

			float pAp = sum(pAp_acc,threads);
			float alpha = rsold / pAp;
			for(int i=from;i<to;i++)  {
				x[i]+=alpha*p[i];
				r[i]-=alpha*Ap[i];
			}
			float myrsnew=0;
			for(int i=from;i<to;i++) 
				myrsnew +=r[i]*r[i];
			rsnew_acc[id] = myrsnew;
			bar.wait();
			float rsnew = sum(rsnew_acc,threads);

			if(rsnew < th2) {
				if(id == 0)
					iters = it;
				return;
			}
			float factor = rsnew / rsold;
			for(int i=from;i<to;i++)
				p[i]=r[i]+factor*p[i];
			rsold = rsnew;
			bar.wait();
		}
	};
	std::thread tasks[16];
	int start = 0;
	int chunk = (N + threads-1)/threads;
	int id=0;
	for(start=0;start<N;start+=chunk) {
		int limit = std::min(start+chunk,N);
		tasks[id] = std::move(std::thread(runner,start,limit,id));
		id++;
	}
	for(int i=0;i<threads;i++) {
		tasks[i].join();
	}
	return iters;
}


#else

std::string solver_name()
{
	return "CPU Single Thread";
}

void sparce_matrix::mpl(int N,float const *vin,float *vout) const
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
		vout[i]=vin[i] + (v[0]*r.v[0]+v[1]*r.v[1]+v[2]*r.v[2]+v[3]*r.v[3]);
	}
}

// CPU single thread
int solve(int N,sparce_matrix const &A,float const *b,float *x,float thresh=1e-8,int limit=-1)
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
		rsnew=0;
		for(int i=0;i<N;i++)  {
			x[i]+=alpha*p[i];
			r[i]-=alpha*Ap[i];
		}
		for(int i=0;i<N;i++) {
			rsnew +=r[i]*r[i];
		}
		if(rsnew < th2)
			break;
		float factor = rsnew / rsold;
		for(int i=0;i<N;i++)
			p[i]=r[i]+factor*p[i];
		rsold = rsnew;


	}
	return it;
}

#endif

