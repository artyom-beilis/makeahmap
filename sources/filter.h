#pragma once
#include "clcpp.h"

void low_pass_filter(int radius,std::vector<std::vector<double> > const &ckernel,std::vector<std::vector<int> > const &safetypes,std::vector<std::vector<int> > &types)
{
	int rows = types.size();
	int cols = types[0].size();
	char const *ocl_kernel_code = 
	R"xxx(
	    __kernel void convolve_color(int rows,int cols,int radius,__global int const *safetypes,__global int *types,__global float const *ckernel)
	    {
		int r=get_global_id(0);
		int c=get_global_id(1);
		if(r > rows || c > cols)
			return;
                int ref_r = r+radius;
                int ref_c = c+radius;
		int safe_cols = cols + radius * 2;
		__global const int *sr_start = safetypes + safe_cols * (ref_r - 1) + ref_c - 1;
                __global const int *sr[3] = { sr_start, sr_start + safe_cols , sr_start + safe_cols*2 };
                int blue =      (sr[0][0] & 0x1) + 2*(sr[0][1] & 0x1) +   (sr[0][2] & 0x1)
                            + 2*(sr[1][0] & 0x1) + 4*(sr[1][1] & 0x1) + 2*(sr[1][2] & 0x1)
                            +   (sr[2][0] & 0x1) + 2*(sr[2][1] & 0x1) +   (sr[2][2] & 0x1);
                blue = (blue * 255 + 8)/ 16;
                if(blue == 255) {
                    types[r*cols + c] = sr[1][1];
		    return;
		}

                float red_sum    = 0;
                float red_weight = 0;
		int ks = 1 + radius * 2;
                for(int kr=0;kr<=ks;kr++) {
                	for(int kc=0;kc<=ks;kc++) {
	                        int type = safetypes[(r+kr)*safe_cols+c+kc];
				int red_type = type >> (16 + 4);
				float weight = ckernel[kr * ks + kc];
				int sum_factor  = (type & 1)  ^ 1; // blue == 0 for blue one of 0xFF or 0
				red_sum    += red_type * weight * sum_factor;
				red_weight += weight * sum_factor;
			}
                }
                int red = (int)(red_sum / red_weight + 0.5)*16;
                int green = (safetypes[ref_r * safe_cols + ref_c] >> 8) & 0xFF;
                types[r*cols + c] = (red << 16) + (green << 8) + blue;
            }
	)xxx";


	context_with_program ctx;
	std::cerr << " Using " << ctx.name() << std::endl;
	ctx.load(ocl_kernel_code);
	kernel filter(ctx);
	filter.assign("convolve_color");
	std::vector<int> local_safetypes((rows+radius*2)*(cols+radius*2));
	for(int r=0;r<rows+2*radius;r++)
		memcpy(&local_safetypes[r*(cols+radius*2)],&safetypes[r],sizeof(int)*(cols+radius*2));
	std::vector<float> local_kernel((2*radius+1)*(2*radius+1));
	int pos =0;
	for(int r=0;r<1+2*radius;r++) 
		for(int c=0;c<1+2*radius;c++)
			local_kernel[pos++]=ckernel[r][c];
		
	std::vector<int> local_types(rows*cols);
	auto start_ts = std::chrono::high_resolution_clock::now();
	memory_object<float> ocl_kernel(ctx,&local_kernel[0],local_kernel.size(),CL_MEM_READ_ONLY);
	memory_object<int> ocl_safetypes(ctx,&local_safetypes[0],local_safetypes.size(),CL_MEM_READ_ONLY);
	memory_object<int> ocl_types(ctx,rows*cols);

	filter.bind_all(rows,cols,radius,ocl_safetypes,ocl_types,ocl_kernel);
	filter.enqueue(rows,cols);
	ocl_types.copy_to_host(&local_types[0],local_types.size());
	auto end_ts = std::chrono::high_resolution_clock::now();
	double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end_ts-start_ts).count();
	std::cerr << "\n\nTime = " << time << "\n";
	for(int r=0;r<rows;r++)
		memcpy(&types[r][0],&local_types[r*cols],sizeof(int)*(cols));
}

