#pragma once
#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <sstream>

#define USE_MT_SOLVER

#ifdef USE_MT_SOLVER
	#include <thread>

	#ifdef _WIN32
		#include <condition_variable>
		#include <mutex>
	#else
		#include <pthread.h>
	#endif

#endif

#ifndef _WIN32
#include <fstream>
#else
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

#endif


namespace cpu  {

#ifndef _WIN32

std::string cpu_name()
{
	std::ifstream proc("/proc/cpuinfo");
	if(!proc)
		return "Unknown CPU";
	std::string line;
	while(std::getline(proc,line)) {
		size_t p;
		if(line.find("model name")!=std::string::npos && (p=line.find(": "))!=std::string::npos) {
			return line.substr(p+2);
		}
	}
	return "Unknown CPU";
	
}

#else

std::string cpu_name()
{
	char name[256]={0};
	DWORD size=sizeof(name)-1;
	HKEY key;
	if(RegOpenKeyEx(HKEY_LOCAL_MACHINE,"HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0",0,KEY_READ,&key))
		return "Unknown CPU";
	if(RegQueryValueExA(key,"ProcessorNameString",0,0,(BYTE*)name,&size)) {
		RegCloseKey(key);
		return "Unknown CPU";
	}
	RegCloseKey(key);
	return name;
}

#endif


template<typename Type>
class kahan_sum {
public:
	kahan_sum() : sum_(), c_() {}
	operator Type () const
	{
		return sum_;
	}

	kahan_sum const &operator+=(Type inp)
	{
		Type y=inp - c_;
		Type t=sum_ + y;
		c_ = (t - sum_)  - y;
		sum_ = t;
		return *this;
	}
private:
	Type sum_;
	Type c_;
};


class sparce_matrix {
public:
	struct row {
		int c[4];
		float v[4];
	};
	void reserve(size_t n) 
	{
		rows_.resize(n);
		for(size_t i=0;i<n;i++) {
			for(int j=0;j<4;j++) {
				rows_[i].c[j]=i;
				rows_[i].v[j]=0;
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
				ind.v[i]=-0.25f;
				break;
			}
		}
		assert(i<4);
	}
	float mpl(int N,float const *vin,float *vout) const;
#ifdef USE_MT_SOLVER
	float mpl(int from,int to,float const *vin,float *vout) const;
#endif	

private:
	std::vector<row> rows_;
};


#if defined USE_MT_SOLVER

std::string solver_name()
{
	std::ostringstream ss;
	ss << "CPU MT solver running on " << cpu_name() << "; " << std::thread::hardware_concurrency() << " threads";
	return ss.str();
}

float sparce_matrix::mpl(int from,int to,float const *vin,float *vout) const
{
	kahan_sum<float> sum;
	for(int i=from;i<to;i++) {
		row r = rows_[i];
		float v[4]= { 
			vin[r.c[0]],  
			vin[r.c[1]],  
			vin[r.c[2]],  
			vin[r.c[3]]
		};
		float x=vin[i];
		float y=x+ (v[0]*r.v[0]+v[1]*r.v[1]+v[2]*r.v[2]+v[3]*r.v[3]);
		sum+=x*y;
		vout[i]=y;
	}
	return sum;
}

float sum(float const *v,int n)
{
	float r=0;
	for(int i=0;i<n;i++)
		r+=v[i];
	return r;
}

#ifndef _WIN32

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

#else

class barrier {
public:
	barrier(int count): threshold_(count), count_(count), generation_(0)
	{
		if (count == 0)
			throw std::invalid_argument("count cannot be zero.");
	}

	void wait()
	{
		std::unique_lock<std::mutex> guard(mutex_);		
		int gen = generation_;

		if (--count_ == 0)
		{
			generation_++;
			count_ = threshold_;
			cond_.notify_all();
			return;
		}

		while (gen == generation_)
			cond_.wait(guard);
	}

private:
	std::mutex mutex_;
	std::condition_variable cond_;
	int threshold_;
	int count_;
	int generation_;
};


#endif

// CPU several threads
std::pair<int,double> solve(int N,sparce_matrix const &A,float const *b,float *x,float thresh=1e-8,int limit=-1)
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
	{
		kahan_sum<float> rsold_s;
		for(int i=0;i<N;i++) {
			p[i] = r[i]=b[i]-r[i];
			rsold_s+=r[i]*r[i];
		}
		grsold = rsold_s;
	}
	
	int threads = std::thread::hardware_concurrency();

	std::vector<float> pAp_acc(threads);
	std::vector<float> rsnew_acc(threads);
	std::vector<std::thread> tasks(threads);
	barrier bar(threads);

	int iters;
	auto runner = [&](int from,int to,int id) {	
		float rsold = grsold;
		for(int it=0;it<limit;it++) {
			float mypAp = A.mpl(from,to,p,Ap);
			pAp_acc[id]=mypAp;
			bar.wait();

			float pAp = sum(pAp_acc.data(),threads);

			float alpha = rsold / pAp;
			kahan_sum<float> myrsnew;
			for(int i=from;i<to;i++) 
				x[i]+=alpha*p[i];
			for(int i=from;i<to;i++) {
				float tmp = r[i]-alpha*Ap[i];
				myrsnew += tmp*tmp;
				r[i]=tmp;
			}
			rsnew_acc[id] = myrsnew;
			bar.wait();
			float rsnew = sum(rsnew_acc.data(),threads);

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
	int start = 0;
	int chunk = (N + threads-1)/threads;
	int id=0;
	auto start_ts = std::chrono::high_resolution_clock::now();
	for(start=0;start<N;start+=chunk) {
		int limit = std::min(start+chunk,N);
		tasks[id] = std::move(std::thread(runner,start,limit,id));
		id++;
	}
	for(int i=0;i<threads;i++) {
		tasks[i].join();
	}
	auto end_ts = std::chrono::high_resolution_clock::now();
	double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end_ts-start_ts).count();
	return std::make_pair(iters,time);
}


#else

std::string solver_name()
{
	return "Single Thread solver running on " + cpu_name();
}

float sparce_matrix::mpl(int N,float const *vin,float *vout) const
{
	int size = rows_.size();
	kahan_sum<float> r;
	for(int i=0;i<size;i++) {
		row const &r = rows_[i];
		float v[4]= { 
			vin[r.c[0]],  
			vin[r.c[1]],  
			vin[r.c[2]],  
			vin[r.c[3]]
		};
		float a = vin[i];
		float b=a + (v[0]*r.v[0]+v[1]*r.v[1]+v[2]*r.v[2]+v[3]*r.v[3]);
		vout[i]=b;
		r+=a*b;
	}
	return r;
}

// CPU single thread
std::pair<int,double> solve(int N,sparce_matrix const &A,float const *b,float *x,float thresh=1e-8,int limit=-1)
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
	kahan_sum<float> rsold_sum;
	for(int i=0;i<N;i++) {
		p[i] = r[i]=b[i]-r[i];
		rsold_sum+=r[i]*r[i];
	}
	rsold = rsold_sum;
	int it;
	float rsnew = 0;	

	auto start_ts = std::chrono::high_resolution_clock::now();
	for(it=0;it<limit;it++) {
		float pAp = A.mpl(N,p,Ap);
		float alpha = rsold / pAp;
		for(int i=0;i<N;i++)  
			x[i]+=alpha*p[i];
		rsnew=0;
		kahan_sum<float> rsnew_sum;
		int id=0;
		int pos=0;
		for(int i=0;i<N;i++) {
			float tmp = r[i] - alpha*Ap[i];
			rsnew_sum += tmp*tmp;
		}
		rsnew = rsnew_sum;
		if(rsnew < th2)
			break;
		float factor = rsnew / rsold;
		for(int i=0;i<N;i++)
			p[i]=r[i]+factor*p[i];
		rsold = rsnew;
	}
	auto end_ts = std::chrono::high_resolution_clock::now();
	double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end_ts-start_ts).count();
	return std::make_pair(it,time);
}

#endif


class eq_solver {
public:
	void init_matrix(int N)
	{
		mat.reserve(N);
		N_ = N;
	}
	void add(int r,int c)
	{
		mat.add(r,c);
	}
	std::pair<int,double> solve(float const *b,float *x,float thresh,int limit)
	{
		return cpu::solve(N_,mat,b,x,thresh,limit);
	}
	std::string name()
	{
		return solver_name();
	}
	static int get_bytes_per_it()
	{
		return sizeof(sparce_matrix::row) + sizeof(float)* ( 2 + 3 + 3 + 3);
	}
	bool is_cpu()
	{
		return true;
	}

private:
	sparce_matrix mat;
	int N_;
};

} // namespace cpu
