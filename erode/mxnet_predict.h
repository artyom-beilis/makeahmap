#pragma once

#include <mxnet/c_predict_api.h>
#include <string>
#include <stdexcept>
#include <fstream>
#include <vector>


class predictor {
public:
	static void init()
	{
	}
	predictor(std::string const &network,std::string const &trained,int batch,int ch,int height,int width) : 
		pred_hnd(0),
		batch_(batch),
		current_batch_(batch),
		chan_(ch),
		height_(height),
		width_(width),
		input_(batch*ch*height*width)
	{
		auto file_network = get_file(network);
		const char* input_key[1] = {"data"};
		const char** input_keys = input_key;
		auto file_param = get_file(trained);

    		const mx_uint input_shape_indptr[2] = { 0, 4 };
   		const mx_uint input_shape_data[4] = { batch,
                                        static_cast<mx_uint>(ch),
                                        static_cast<mx_uint>(height),
                                        static_cast<mx_uint>(width)};


		MXPredCreate(	file_network.data(),
				file_param.data(),
				file_param.size(),
				2, // gpu
				0, // dev_id
				1, // num_input_nodes
				input_keys,
		                input_shape_indptr,
                 		input_shape_data,
                 		&pred_hnd);
		if(!pred_hnd)
			throw std::runtime_error("Failed to create predictor: " + std::string(MXGetLastError()));
	}
	~predictor()
	{
		MXPredFree(pred_hnd);
	}

	int channels() 
	{
		return chan_;
	}
	int width() 
	{
		return width_;
	}
	int height() 
	{
		return height_;
	}
	void reshape(int batch)
	{
		if(batch > batch_)
			throw std::runtime_error("Batch size invaid");
		current_batch_ = batch;
	}
	float *input()
	{
		return input_.data();
	}
	void forward()
	{
		MXPredSetInput(pred_hnd, "data", input_.data(), input_.size());
    		MXPredForward(pred_hnd);
		
		mx_uint output_index = 0;

		mx_uint *shape = 0;
		mx_uint shape_len;

		// Get Output Result
		MXPredGetOutputShape(pred_hnd, output_index, &shape, &shape_len);

		size_t size = 1;
		for (mx_uint i = 0; i < shape_len; ++i) size *= shape[i];

		output_.resize(size);

		MXPredGetOutput(pred_hnd, output_index, &(output_[0]), size);
	}
	int classes()
	{
		return output_.size() / batch_;
        }
	int num()
	{
		return current_batch_;
	}
	
	float const *result()
	{
		return output_.data();
	}

private:
	std::vector<char> get_file(std::string const &name)
	{
		std::ifstream f(name.c_str(),std::ifstream::binary);
		if(!f)
			throw std::runtime_error("Failed to open file " + name);
		f.seekg(0,std::ifstream::end);
		size_t len = f.tellg();
		std::vector<char> r(len);
		f.seekg(0);
		f.read(r.data(),r.size());
		return r;
	}
    	PredictorHandle pred_hnd;
	int batch_,current_batch_,chan_,height_,width_;
	std::vector<float> input_,output_;
};
