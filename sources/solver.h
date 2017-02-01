#pragma once
#include <vector>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <thread>

#ifdef _WIN32
#  include <condition_variable>
#  include <mutex>
#else
#  include <pthread.h>
#endif

#ifndef _WIN32
#  include <fstream>
#else
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  include <windows.h>
#endif

#define OPTIMIZE_4BLK

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


class sparce_matrix {
public:
	union row {
		short c[4];
		unsigned long long v;
		bool operator==(row const &other) const 
		{
			return v==other.v;
		}
		bool operator!=(row const &other) const
		{
			return v!=other.v;
		}
		bool full() const
		{
			return c[0]!=0 && c[1]!=0 && c[2]!=0 && c[3]!=0;
		}
	};
	void reserve(int threads,size_t n) 
	{
		threads_ = threads;
		assert(n % (threads * 4)==0);
		step_ = n / threads_;
		rows_.resize(threads_);
		ratio_.resize(threads_,1.0);
		for(int i=0;i<threads_;i++) {
			rows_[i].resize(step_ + 1);
		}
		for(size_t i=0;i<rows_.size();i++) {
			for(int j=0;j<step_ + 1;j++) {
				for(int k=0;k<4;k++) {
					rows_[i][j].c[k]=0;
				}
			}
		}
	}
	void add(int r,int c)
	{
		row &ind = rows_[r / step_][r % step_];
		int diff = c - r;
		assert(abs(diff) <= 4096);
		if(diff < -1) {
			assert(ind.c[0]==0);
			ind.c[0] = -diff * 2 + 1;
		}
		else if(diff  == -1) {
			assert(ind.c[1]==0);
			ind.c[1] = 2 + 1;
		}
		else if(diff == 1) {
			assert(ind.c[2]==0);
			ind.c[2] = 2 + 1;
		}
		else if(diff > 1) {
			assert(ind.c[3]==0);
			ind.c[3] = diff * 2 + 1;
		}
		else {
			assert(!"Not get there!");
		}
	}
	void optimize()
	{
		#ifdef OPTIMIZE_4BLK
		compressed_.resize(threads_);
		for(int id=0;id<threads_;id++)
			optimize(id);
		#endif
	}
	int find_seq(int pos,std::vector<row> const &rows)
	{
		row r = rows[pos];
		if(!r.full())
			return pos;
		int e;
		for(e=pos+1;e<step_;e++) {
			if(rows[e]!=r)
				break;
		}
		e=e / 4 * 4;
		return e;
	}
	#ifdef OPTIMIZE_4BLK
	void optimize(int id)
	{
		std::vector<std::pair<int,int> > &cmp = compressed_[id];
		std::vector<row> const &rows = rows_[id];
		cmp.clear();
		cmp.reserve(step_);
		int start = 0;
		int bulk=0;
		while(start < step_) {
			int end = find_seq(start,rows);
			if(end == start) {
				start += 4;
				continue;
			}
			else if(end - start < 16 || start == 0 || end == step_) {
				start = end;
				continue;
			}
			else {
				cmp.push_back(std::make_pair(start,end));
				bulk+=end-start;
				start = end;
				continue;
			}
		}
		cmp.push_back(std::make_pair(step_,step_));
		if(step_ != 0) {
			ratio_[id] = double(step_ - bulk) / step_;
		}
	}
	#endif
	double mpl(int id,float const *vin,float *vout) const;
	double mpl(float const *vin,float *vout) const
	{
		double s=0;
		for(int i=0;i<int(rows_.size());i++)
			s+=mpl(i,vin,vout);
		return s;
	}
	double ratio() const
	{
		double s=0;
		for(size_t i=0;i<ratio_.size();i++)
			s+=ratio_[i];
		return s / ratio_.size();
	}
private:
	int threads_;
	int step_;
	std::vector<std::vector<row> > rows_;
	std::vector<double> ratio_;
	std::vector<std::vector<std::pair<int,int> > > compressed_;
};


std::string solver_name()
{
	std::ostringstream ss;
	ss << "CPU Solver running on " << cpu_name() << "; " << std::thread::hardware_concurrency() << " threads";
	return ss.str();
}

#ifdef OPTIMIZE_4BLK
double sparce_matrix::mpl(int id,float const *vin,float *vout) const
{
	std::vector<row> const &my_row = rows_[id];
	std::vector<std::pair<int,int> > const &cmp = compressed_[id];
	int cur_cmp = 0;
	int i=id * step_;
	int pos = 0;
	double sum = 0.0;
	while(pos < step_) {
		// mixed		
		std::pair<int,int> crange = cmp[cur_cmp++];

		for(;pos < crange.first;pos ++, i++) {
			row r = my_row[pos];
			float v[4]= { 
				vin[i - int(r.c[0] >>1)] * (r.c[0] & 1),  
				vin[i - int(r.c[1] & 1)] * (r.c[1] & 1),
				vin[i + int(r.c[2] & 1)] * (r.c[2] & 1),
				vin[i + int(r.c[3] >>1)] * (r.c[3] & 1),  
			};

			float x=vin[i];
			float y= x - (v[0]+v[1]+v[2]+v[3])/4;
			sum+=x*y;
			vout[i]=y;

		}

		// fixed 4 blocks
		// it is ok to pos=step_ - we have extra place
		row r = my_row[pos];
		int up = int(r.c[0]>>1);
		int dn = int(r.c[3]>>1);

		#ifdef USE_SIMD
		typedef float float4  __attribute__((vector_size(16)));
		typedef float float4u __attribute__((vector_size(16),aligned(4)));

		for(;pos < crange.second;pos+=4,i+=4) {
			float4 x=*(float4*)&vin[i];
			float4 v[4]= { 
				*(float4u*)&vin[i - up],  
				*(float4u*)&vin[i -  1],
				*(float4u*)&vin[i +  1],
				*(float4u*)&vin[i + dn], 
			};

			float4 y = x - 0.25f * (v[0]+v[1]+v[2]+v[3]);
			float4 sum4 = x*y;
			sum += (sum4[0]+sum4[1]) + (sum4[2]+sum4[3]);
			*(float4*)&vout[i] = y;
		}

		#else
		for(;pos < crange.second;pos++,i++) {
			float x=vin[i];
			float v[4]= { 
				vin[i - up],  
				vin[i -  1],
				vin[i +  1],
				vin[i + dn], 
			};

			float y = x - (v[0]+v[1]+v[2]+v[3]) / 4;
			sum+=x*y;
			vout[i]=y;
		}
		#endif

	}
	return sum;
}
#else
double sparce_matrix::mpl(int id,float const *vin,float *vout) const
{
	int from = id * step_;
	int to   = from + step_;
	double sum = 0.0;
	std::vector<row> const &my_row = rows_[id];
	for(int i=from,pos=0;i<to;i++,pos++) {
		row r = my_row[pos];
		float v[4]= { 
			vin[i - int(r.c[0]>>1)] * (r.c[0] & 1),  
			vin[i - int(r.c[1] & 1)] * (r.c[1] & 1),
			vin[i + int(r.c[2] & 1)] * (r.c[2] & 1),
			vin[i + int(r.c[3]>>1)] * (r.c[3] & 1),  
		};
		float x=vin[i];
		float y=x-(v[0]+v[1]+v[2]+v[3])/4;
		sum+=x*y;
		vout[i]=y;
	}
	return sum;
}
#endif

double sum(double const *v,int n)
{
	double r=0;
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

	A.mpl(x,r);
	double grsold=0;
	{
		double rsold_s=0.0;
		for(int i=0;i<N;i++) {
			p[i] = r[i]=b[i]-r[i];
			rsold_s+=r[i]*r[i];
		}
		grsold = rsold_s;
	}

	int threads = std::thread::hardware_concurrency();

	std::vector<double> pAp_acc(threads);
	std::vector<double> rsnew_acc(threads);
	std::vector<std::thread> tasks(threads);
	barrier bar(threads);

	
	#ifdef USE_SIMD
	typedef float float_type __attribute__((vector_size(16)));
	#define float_type_sum(x) ( ((x)[0] + (x)[1]) + ((x)[2] + (x)[3]) )
	static const int op_factor = 4;
	#else
	typedef float float_type;
	#define float_type_sum(x) (x)
	static const int op_factor = 1;
	#endif
	int iters;

	auto runner = [&](int from,int to,int id) {	
		double rsold = grsold;
		int items = (to - from) / op_factor;
		float_type *x_p  = (float_type*)(x+from);
		float_type *r_p  = (float_type*)(r+from);
		float_type *p_p  = (float_type*)(p+from);
		float_type *Ap_p = (float_type*)(Ap + from);
		for(int it=0;it<limit;it++) {
			double mypAp = A.mpl(id,p,Ap);
			pAp_acc[id]=mypAp;
			bar.wait();

			double pAp = sum(pAp_acc.data(),threads);

			float alpha = rsold / pAp;
			double myrsnew=0.0;
			for(int i=0;i<items;i++) 
				x_p[i]+=alpha*p_p[i];
			for(int i=0;i<items;i++) {
				float_type tmp = r_p[i]-alpha*Ap_p[i];
				float_type tmp2 = tmp*tmp;
				myrsnew += float_type_sum(tmp2);
				r_p[i]=tmp;
			}
			rsnew_acc[id] = myrsnew;
			bar.wait();
			double rsnew = sum(rsnew_acc.data(),threads);

			if(rsnew < th2) {
				if(id == 0)
					iters = it;
				return;
			}
			float factor = rsnew / rsold;
			for(int i=0;i<items;i++)
				p_p[i]=r_p[i]+factor*p_p[i];
			rsold = rsnew;
			bar.wait();
		}
	};
	int start = 0;
	int chunk = (N + threads-1)/threads;
	int id=0;
	auto start_ts = std::chrono::high_resolution_clock::now();
	if(grsold > 1e-5) {
		for(start=0;start<N;start+=chunk) {
			int limit = std::min(start+chunk,N);
			tasks[id] = std::move(std::thread(runner,start,limit,id));
			id++;
		}
		for(int i=0;i<threads;i++) {
			tasks[i].join();
		}
	}
	else {
		iters = 1;
	}
	auto end_ts = std::chrono::high_resolution_clock::now();
	double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end_ts-start_ts).count();
	return std::make_pair(iters,time);
}



class eq_solver {
public:
	void init_matrix(int N)
	{
		mat.reserve(std::thread::hardware_concurrency(),N);
		N_ = N;
	}
	void add(int r,int c)
	{
		mat.add(r,c);
	}
	std::pair<int,double> solve(float const *b,float *x,float thresh,int limit)
	{
		mat.optimize();
		return cpu::solve(N_,mat,b,x,thresh,limit);
	}
	std::string name()
	{
		return solver_name();
	}
	double get_bytes_per_it()
	{
		return sizeof(sparce_matrix::row) * mat.ratio() + sizeof(float)* (2 + 3 + 3 + 3);
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
