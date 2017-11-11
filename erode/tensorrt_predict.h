#pragma once

#include <cuda_runtime_api.h>

#include "NvInfer.h"
#include "NvCaffeParser.h"

#include <iostream>
#include <vector>
#include <chrono>

using namespace nvinfer1;
using namespace nvcaffeparser1;


class predictor {
public:
	// Logger for GIE info/warning/errors
	class Logger : public ILogger			
	{
		void log(Severity severity, const char* msg) override
		{
			// suppress info-level messages
			if (severity != Severity::kINFO)
				std::cerr << msg << std::endl;
		}
	};


	~predictor()
	{
		context_->destroy();
		engine_->destroy();
		runtime_->destroy();
	}
	static void init()
	{
	}
	predictor(std::string const &net,std::string const &trained,int b,int c,int h,int w,char const *layer_name=nullptr) 
	{
    		IHostMemory *gieModelStream = nullptr;
		IBuilder* builder = createInferBuilder(logger_);
		INetworkDefinition* network = builder->createNetwork();
		ICaffeParser* parser = createCaffeParser();
		IBlobNameToTensor const *blobNameToTensor = parser->parse(net.c_str(),trained.c_str(),*network,DataType::kFLOAT);

		ITensor *output_blob = nullptr;
		if(layer_name) {
			output_blob = blobNameToTensor->find(layer_name);
		}
		else {
			output_blob = blobNameToTensor->find("prob");
			if(!output_blob)
				output_blob = blobNameToTensor->find("softmax");
			if(!output_blob)
				output_blob = blobNameToTensor->find("output");
		}
		if(!output_blob)
			throw std::runtime_error("Failed to find output layer");
		network->markOutput(*output_blob);
		nvinfer1::Dims output_dim = output_blob->getDimensions();
		classes_ = 1;
		for(int i=0;i<output_dim.nbDims;i++) {
			classes_ *= output_dim.d[i];
		}
		printf("\n");
		builder->setMaxBatchSize(b);
		builder->setMaxWorkspaceSize(1 << 20);
		ICudaEngine* engine = builder->buildCudaEngine(*network);
		if(!engine)
			throw std::runtime_error("Failed to build engine");
		network->destroy();
		parser->destroy();
		gieModelStream = engine->serialize();
		engine->destroy();
		builder->destroy();
		shutdownProtobufLibrary();
		runtime_ = createInferRuntime(logger_);
		engine_ = runtime_->deserializeCudaEngine(gieModelStream->data(), gieModelStream->size(), nullptr);
    		gieModelStream->destroy();
		context_ = engine_->createExecutionContext();
		batch_=curr_batch_=b;
		chan_=c;
		height_=h;
		width_=w;
		input_.resize(batch_*chan_*height_*width_);
		output_.resize(batch_*classes_);
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
			throw std::runtime_error("Invalid reshape size");
		curr_batch_ = batch;
	}
	float *input()
	{
		return input_.data();
	}
	void forward()
	{
    		auto start = std::chrono::high_resolution_clock::now();
		ICudaEngine const &engine = context_->getEngine();
		void *dev_buffers[2];
		int input_index = engine.getBindingIndex("data");
		int output_index = 1 - input_index;
		int input_size = curr_batch_ * chan_ * height_ * width_;
		int output_size = curr_batch_ * classes_;
		cudaMalloc(&dev_buffers[input_index],input_size * sizeof(float));
		cudaMalloc(&dev_buffers[output_index],output_size * sizeof(float));

		cudaStream_t stream;
		cudaStreamCreate(&stream);
		cudaMemcpyAsync(dev_buffers[input_index], input_.data(), input_size * sizeof(float), cudaMemcpyHostToDevice, stream);
		context_->enqueue(curr_batch_, dev_buffers, stream, nullptr);
		cudaMemcpyAsync(output_.data(), dev_buffers[output_index], output_size*sizeof(float), cudaMemcpyDeviceToHost, stream);
		cudaStreamSynchronize(stream);
		

		cudaStreamDestroy(stream);
		cudaFree(dev_buffers[input_index]);
		cudaFree(dev_buffers[output_index]);
	}
	int classes()
	{
		return classes_;
        }
	int num()
	{
		return curr_batch_;
	}
	
	float const *result()
	{
		return output_.data();
	}

private:
	Logger logger_;
	IRuntime* runtime_;
	ICudaEngine* engine_;
	IExecutionContext *context_;

	int batch_,curr_batch_,chan_,height_,width_,classes_;
	std::vector<float> input_,output_;

};
