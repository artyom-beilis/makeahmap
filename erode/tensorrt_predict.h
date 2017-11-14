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
		cudaFreeHost(input_);
		cudaFreeHost(output_);
		cudaFree(dev_buffers_[0]);
		cudaFree(dev_buffers_[1]);
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
		IBlobNameToTensor const *blobNameToTensor = 0;
		if(getenv("ENABLE_HALF_MODE")) {
			blobNameToTensor = parser->parse(net.c_str(),trained.c_str(),*network,DataType::kHALF);
			builder->setHalf2Mode(true);
		}
		else {
			blobNameToTensor = parser->parse(net.c_str(),trained.c_str(),*network,DataType::kFLOAT);
		}


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
		void *p1,*p2;
		cudaMallocHost(&p1,batch_*chan_*height_*width_*sizeof(float));
		cudaMallocHost(&p2,batch_*classes_*sizeof(float));
		input_ = static_cast<float*>(p1);
		output_ = static_cast<float*>(p2);

		{
			ICudaEngine const &engine = context_->getEngine();
			input_index_ = engine.getBindingIndex("data");
			output_index_ = 1 - input_index_;
			int input_size = batch_ * chan_ * height_ * width_;
			int output_size = batch_ * classes_;
			cudaMalloc(&dev_buffers_[input_index_],input_size * sizeof(float));
			cudaMalloc(&dev_buffers_[output_index_],output_size * sizeof(float));
		}
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
		return input_;
	}
	void forward()
	{
	
		int input_size = width_*height_*chan_*curr_batch_;
		int output_size = classes_ * curr_batch_;
		cudaStream_t stream;
		cudaStreamCreate(&stream);
    		auto start = std::chrono::high_resolution_clock::now();
		cudaMemcpyAsync(dev_buffers_[input_index_], input_, input_size * sizeof(float), cudaMemcpyHostToDevice, stream);
		context_->enqueue(curr_batch_, dev_buffers_, stream, nullptr);
		cudaMemcpyAsync(output_, dev_buffers_[output_index_], output_size*sizeof(float), cudaMemcpyDeviceToHost, stream);
    		auto e1 = std::chrono::high_resolution_clock::now();
		cudaStreamSynchronize(stream);
    		auto e2 = std::chrono::high_resolution_clock::now();
    		double t1 = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(e1-start).count();
    		double t2 = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(e2-start).count();
		//std::cerr << t1 << " " << t2 << std::endl;
		

		cudaStreamDestroy(stream);
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
		return output_;
	}

private:
	Logger logger_;
	IRuntime* runtime_;
	ICudaEngine* engine_;
	IExecutionContext *context_;

	int batch_,curr_batch_,chan_,height_,width_,classes_;
	int input_index_,output_index_;
	float *input_,*output_;
	void *dev_buffers_[2];

};
