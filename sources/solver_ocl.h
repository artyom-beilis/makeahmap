#include "clcpp.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <assert.h>
#include <memory>
#include <chrono>

namespace gpu {

#define LOCAL_BLOCK	128
#define LOCAL_BLOCK_STR "128"

class eq_solver { 
public:
	struct row {
		cl_short4   c;
	};
	void init_matrix(size_t n)
	{
		rows_.resize(n);
		for(size_t i=0;i<n;i++) {
			for(int j=0;j<4;j++) {
				rows_[i].c.s[j]=0;
			}
		}
	}
	std::vector<row> rows_;
	void add(int r,int c)
	{
		row &ind = rows_[r];
		int i;
		assert(labs(r-c) < 32767);
		for(i=0;i<4;i++) {
			if(ind.c.s[i]==0) {
				ind.c.s[i]=c - r;
				break;
			}
		}
		assert(i<4);
	}
	context_with_program &context() { return ctx_; }
	eq_solver() : 
		mpl_(ctx_),
		mpl_prod_(ctx_),
		dot_prod_(ctx_),
		vsum_(ctx_),
		b_min_r_to_r_and_p_(ctx_),add_vectors_(ctx_),add_vectors_and_prod_(ctx_),add_vectors_to_(ctx_)
	{
		ctx_.load(
			"#define LOCAL_BLOCK " LOCAL_BLOCK_STR	"\n"
			R"xxx(
			typedef struct row {
				short4 c;
			} row;
			__kernel void matrix_mpl(
				int n,
				__global const row *rows,
				__global const float *vin,
				__global float *vout)
			{                                                                       
				int r = get_global_id(0);                                     
				if(r >= n) 
					return;
				__global const row *cr = rows + r;
				float4 cvA = (float4)( vin[cr->c.x + r],vin[cr->c.y + r],vin[cr->c.z + r ],vin[cr->c.w + r] );
				float4 cvB = (float4)( 
					(cr->c.x != 0 ? -0.25f : 0.0),
					(cr->c.y != 0 ? -0.25f : 0.0),
					(cr->c.z != 0 ? -0.25f : 0.0),
					(cr->c.w != 0 ? -0.25f : 0.0) 
				);
				vout[r] = vin[r] + dot(cvA,cvB);
			}

			void reduce(__local float *l_sum,__global float *r)
			{
				barrier(CLK_LOCAL_MEM_FENCE);
				int lid = get_local_id(0);  
				int wgid = get_group_id(0); 
				int dist = get_local_size(0);
				while(dist > 1) { 
					dist >>=1;
					if(lid < dist)  {
						l_sum[lid]+=l_sum[lid+dist];
					}
					barrier(CLK_LOCAL_MEM_FENCE);
				}
				if(lid == 0)
					r[wgid] = l_sum[0];
			}

			__kernel void matrix_mpl_and_prod(
				int n,
				__global const row *rows,
				__global const float *vin,
				__global float *vout,
				__global float *prod)
			{                          
				__local float l_sum[LOCAL_BLOCK];
				int r = get_global_id(0);
				int lid = get_local_id(0);  
				if(r < n) {
					__global const row *cr = rows + r;
					float4 cvA = (float4)( vin[cr->c.x + r],vin[cr->c.y + r],vin[cr->c.z + r ],vin[cr->c.w + r] );
					float4 cvB = (float4)( 
						(cr->c.x != 0 ? -0.25f : 0.0),
						(cr->c.y != 0 ? -0.25f : 0.0),
						(cr->c.z != 0 ? -0.25f : 0.0),
						(cr->c.w != 0 ? -0.25f : 0.0) 
					);
					float p = vin[r];
					float Ap = p + dot(cvA,cvB);
					vout[r] = Ap;
					l_sum[lid] = p*Ap;
				}
				else
					l_sum[lid] = 0.0f;
				reduce(l_sum,prod);
			}
			__kernel void dot_prod(
				int N,
				__global const float *A,
				__global const float *B,
				__global float *r)
			{ 
				 
				int gid = get_global_id(0);
				int lid = get_local_id(0);
				__local float l_sum[LOCAL_BLOCK];
				if(gid < N) 
					l_sum[lid] = A[gid]*B[gid]; 
				else
					l_sum[lid] = 0;
				reduce(l_sum,r);
			}

			__kernel void b_min_r_to_r_and_p(int N,__global const float *b,__global float *r,__global float *p)
			{
				int id = get_global_id(0);
				if(id < N) {
					float diff = b[id]-r[id];
					r[id]=diff;
					p[id]=diff;
				}
			}
			__kernel void vsum(int N,__global float const *A,__global float *R)
			{
				__local float l_sum[LOCAL_BLOCK];
				int gid = get_global_id(0);
				int lid = get_local_id(0);
				if(gid < N) 
					l_sum[lid] = A[gid];
				else
					l_sum[lid] = 0;
				reduce(l_sum,R);
			}
			__kernel void add_vectors(int N,__global float *x,__global float const *b,float factor)
			{
				int gid = get_global_id(0);
				if(gid < N)
					x[gid]+=b[gid]*factor;
			}
			__kernel void add_vectors_and_prod(int N,__global float *x,__global float const *b,float factor,__global float *prod)
			{
				int gid = get_global_id(0);
				int lid = get_local_id(0);
				__local float l_sum[LOCAL_BLOCK];
				if(gid < N) {
					float val = x[gid] + b[gid]*factor;
					x[gid]=val;
					l_sum[lid] = val*val;
				}
				else {
					l_sum[lid]=0;
				}
				reduce(l_sum,prod);
			}
			__kernel void add_vectors_to(int N,__global float *me,float factor,__global float const *other)
			{
				int gid = get_global_id(0);
				if(gid < N)
					me[gid]=me[gid]*factor + other[gid];
			}
		)xxx");

		mpl_.assign("matrix_mpl");
		mpl_.local(LOCAL_BLOCK);

		mpl_prod_.assign("matrix_mpl_and_prod");
		mpl_prod_.local(LOCAL_BLOCK);

		dot_prod_.assign("dot_prod");
		dot_prod_.local(LOCAL_BLOCK);

		vsum_.assign("vsum");
		vsum_.local(LOCAL_BLOCK);

		b_min_r_to_r_and_p_.assign("b_min_r_to_r_and_p");
		b_min_r_to_r_and_p_.local(LOCAL_BLOCK);
		add_vectors_.assign("add_vectors");
		add_vectors_.local(LOCAL_BLOCK);

		add_vectors_and_prod_.assign("add_vectors_and_prod");
		add_vectors_and_prod_.local(LOCAL_BLOCK);

		add_vectors_to_.assign("add_vectors_to");
		add_vectors_to_.local(LOCAL_BLOCK);
	}
	float post_reduce(memory_object<float> &rsrc,memory_object<float> &rtgt,int N)
	{
		memory_object<float> *src=&rsrc;
		memory_object<float> *tgt=&rtgt;
		int res_size = (N + LOCAL_BLOCK - 1) / LOCAL_BLOCK;
		while(res_size > 1) {
			vsum_(res_size,*src,*tgt);
			std::swap(src,tgt);
			res_size = (res_size + LOCAL_BLOCK - 1) / LOCAL_BLOCK;
		}
		float r=0.0;
		src->copy_to_host(&r,1);
		return r;
		
	}

	std::string name()
	{
		return "OpenCL solver running on " + context().name();
	}

	bool is_cpu()
	{
		return context().is_cpu();
	}

	std::pair<int,double>
	solve(float const *b_in,float *x0_in,float thresh=1e-6,int limit=-1)
	{
		int N = rows_.size();
		if(limit==-1)
			limit = rows_.size();
		memory_object<row> Mrows(ctx_,rows_.data(),rows_.size(),CL_MEM_READ_ONLY);
		memory_object<float> b(ctx_,b_in,N,CL_MEM_READ_ONLY);
		memory_object<float> r(ctx_,N);
		memory_object<float> Ap(ctx_,N);
		memory_object<float> p(ctx_,N);
		memory_object<float> x(ctx_,x0_in,N);
		size_t r1 = (N+LOCAL_BLOCK-1)/LOCAL_BLOCK;
		size_t r2 = (r1+LOCAL_BLOCK-1)/LOCAL_BLOCK;
		memory_object<float> rd1(ctx_,r1);
		memory_object<float> rd2(ctx_,r2);

		float th2=thresh*thresh;
		mpl_(N,Mrows,x,r);
		b_min_r_to_r_and_p_(N,b,r,p);

		dot_prod_(N,r,r,rd1);
		float rsold = post_reduce(rd1,rd2,N);
		int it=0;
		float rsnew = 0;
		auto start = std::chrono::high_resolution_clock::now();
		if(rsold >= 1e-5) {
			for(it=0;it<limit;it++) {
				mpl_prod_(N,Mrows,p,Ap,rd1);
				float pAp = post_reduce(rd1,rd2,N);
				float alpha = rsold / pAp;
				add_vectors_(N,x,p,alpha);
				add_vectors_and_prod_(N,r,Ap,-alpha,rd1);
				rsnew = post_reduce(rd1,rd2,N);
				if(rsnew < th2) 
					break;
				float factor = rsnew / rsold;
				add_vectors_to_(N,p,factor,r);
				rsold = rsnew;
			}
		}
		else {
			it = 1;
		}
		auto end = std::chrono::high_resolution_clock::now();
		double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end-start).count();
		x.copy_to_host(x0_in,N);
		return std::make_pair(it,time);
	}
	
	static int get_bytes_per_it()
	{
		return sizeof(eq_solver::row) + sizeof(float)* ( 2 + 3 + 3 + 3);
	}

	context_with_program ctx_;
	kernel mpl_,mpl_prod_;
	kernel dot_prod_;
	kernel vsum_;
	kernel b_min_r_to_r_and_p_,add_vectors_,add_vectors_and_prod_,add_vectors_to_;
};


} // gpu
