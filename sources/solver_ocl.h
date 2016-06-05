#include "clcpp.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <iostream>
#include <assert.h>
#include <memory>

#define SRC(x) #x

class eq_solver { 
public:
	struct row {
		cl_int4   c;
		cl_float4 m;
	};
	void init_matrix(size_t n)
	{
		rows_.resize(n);
		for(size_t i=0;i<n;i++) {
			for(int j=0;j<4;j++) {
				rows_[i].c.s[j]=i;
				rows_[i].m.s[j]=0.0f;
			}
		}
	}
	std::vector<row> rows_;
	int prev_r_;
	void add(int r,int c)
	{
		row &ind = rows_[r];
		int i;
		for(i=0;i<4;i++) {
			if(ind.c.s[i]==r) {
				ind.c.s[i]=c;
				ind.m.s[i]=-0.25f;
				break;
			}
		}
		assert(i<4);
	}
	context_with_program &context() { return ctx_; }
	eq_solver() : 
		mpl_(ctx_),
		dot_prod_(ctx_),
		vsum_(ctx_),
		b_min_r_to_r_and_p_(ctx_),add_vectors_(ctx_),add_vectors_to_(ctx_),
		last_res_size_(),last_local_size_()
	{
		ctx_.load(R"xxx(
			typedef struct row {
				int4 c;
				float4 m;
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
				float4 cv = (float4)( vin[cr->c.x],vin[cr->c.y],vin[cr->c.z],vin[cr->c.w] );
				vout[r] = vin[r] + dot(cr->m,cv);
			};
			__kernel void dot_prod(
				int N,
				__global const float *A,
				__global const float *B,
				__local float *l_sum,
				__global float *r)
			{  
				int gid = get_global_id(0);
				int lid = get_local_id(0);  
				int wgid = get_group_id(0); 
				int dist = get_local_size(0);
				if(gid < N) 
					l_sum[lid] = A[gid]*B[gid]; 
				else
					l_sum[lid] = 0;
				barrier(CLK_LOCAL_MEM_FENCE);
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

			__kernel void b_min_r_to_r_and_p(int N,__global const float *b,__global float *r,__global float *p)
			{
				int id = get_global_id(0);
				if(id < N) {
					float diff = b[id]-r[id];
					r[id]=diff;
					p[id]=diff;
				}
			}
			__kernel void vsum(int N,__global float const *A,__local float *l_sum,__global float *R)
			{
				int gid = get_global_id(0);
				int lid = get_local_id(0);  
				int wgid = get_group_id(0); 
				int dist = get_local_size(0);
				if(gid < N) 
					l_sum[lid] = A[gid];
				else
					l_sum[lid] = 0;
				barrier(CLK_LOCAL_MEM_FENCE);
				while(dist > 1) { 
					dist >>=1;
					if(lid < dist)  {
						l_sum[lid]+=l_sum[lid+dist];
					}
					barrier(CLK_LOCAL_MEM_FENCE);
				}
				if(lid == 0)
					R[wgid] = l_sum[0];
			}
			__kernel void add_vectors(int N,__global float *x,__global float const *b,float factor)
			{
				int gid = get_global_id(0);
				if(gid < N)
					x[gid]+=b[gid]*factor;
			}
			__kernel void add_vectors_to(int N,__global float *me,float factor,__global float const *other)
			{
				int gid = get_global_id(0);
				if(gid < N)
					me[gid]=me[gid]*factor + other[gid];
			}
		)xxx");

		mpl_.assign("matrix_mpl");
		dot_prod_.assign("dot_prod");
		vsum_.assign("vsum");
		b_min_r_to_r_and_p_.assign("b_min_r_to_r_and_p");
		add_vectors_.assign("add_vectors");
		add_vectors_to_.assign("add_vectors_to");
	}
	float dot_product(memory_object<float> &a,memory_object<float> &b,int N)
	{
		static size_t local = 0;
		if(local == 0)
			local = dot_prod_.get_local_workgroup();
		size_t res_size  = (N + local - 1)/local;
		if(last_res_size_ != res_size) {
			p1_.reset(new memory_object<float>(ctx_,res_size));
			last_res_size_ = res_size;
		}
		if(last_local_size_ != local) {
			p2_.reset(new memory_object<float>(ctx_,(res_size+local-1)/local));
			last_local_size_ = local;
		}

		dot_prod_.local(local);
		dot_prod_(N,a,b,local_memory(sizeof(float)*local),*p1_);

		memory_object<float> *src= p1_.get(), *dst = p2_.get();
		
		while(res_size > 1) {
			vsum_.local(local);
			vsum_(res_size,*src,local_memory(sizeof(float)*local),*dst);
			std::swap(src,dst);
			res_size = (res_size + local - 1) / local;
		}
		float s1=0.0;
		src->copy_to_host(&s1,1);
		return s1;
	}

	int solve(float const *b_in,float *x0_in,float thresh=1e-6,int limit=-1)
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

		float th2=thresh*thresh;
		mpl_(N,Mrows,x,r);
		b_min_r_to_r_and_p_(N,b,r,p);
		

		float rsold = dot_product(r,r,N);
		int it=0;
		float rsnew = 0;
		for(it=0;it<limit;it++) {
			mpl_(N,Mrows,p,Ap);
			float pAp = dot_product(p,Ap,N);
			float alpha = rsold / pAp;
			add_vectors_(N,x,p,alpha);
			add_vectors_(N,r,Ap,-alpha);
			rsnew = dot_product(r,r,N);
			if(rsnew < th2) 
				break;
			float factor = rsnew / rsold;
			add_vectors_to_(N,p,factor,r);
			rsold = rsnew;
		}
		x.copy_to_host(x0_in,N);
		return it;
	}

	context_with_program ctx_;
	kernel mpl_;
	kernel dot_prod_;
	kernel vsum_;
	kernel b_min_r_to_r_and_p_,add_vectors_,add_vectors_to_;
	std::unique_ptr<memory_object<float> > p1_,p2_;
	size_t last_res_size_,last_local_size_;
};

