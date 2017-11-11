#pragma once

#include <caffe/caffe.hpp>


class predictor {
public:
	static void init()
	{
		using namespace caffe;
		Caffe::set_mode(Caffe::GPU);
		Caffe::SetDevice(0);
	}
	predictor(std::string const &network,std::string const &trained,int b,int c,int h,int w) :
#ifdef USE_OCL
		net(network, caffe::TEST, caffe::Caffe::GetDefaultDevice())
#else
		net(network, caffe::TEST)
#endif
	{
    		net.CopyTrainedLayersFrom(trained);
		if(b < net.input_blobs()[0]->num() || c!=channels() || h!=height() || w != width())
			throw std::runtime_error("Invalid network topology");
	}
	int channels() 
	{
		return net.input_blobs()[0]->channels();
	}
	int width() 
	{
		return net.input_blobs()[0]->width();
	}
	int height() 
	{
		return net.input_blobs()[0]->height();
	}
	void reshape(int batch)
	{
		int c = channels();
		int h = height();
		int w = width();
		net.input_blobs()[0]->Reshape(batch,c,h,w);
		net.Reshape();
	}
	float *input()
	{
    		caffe::Blob<float> *input  = net.input_blobs()[0];
    		return input->mutable_cpu_data();
	}
	void forward()
	{
		net.Reshape();
		net.Forward();
	}
	int classes()
	{
	    caffe::Blob<float> *output = net.output_blobs()[0];
    	    int classes = output->width()  * output->height() * output->channels();
	    return classes;
        }
	int num()
	{
	    caffe::Blob<float> *output = net.output_blobs()[0];
	    return output->num();
	}
	
	float const *result()
	{
    	    	caffe::Blob<float> *output = net.output_blobs()[0];
    	    	return output->cpu_data();
	}

private:
	 caffe::Net<float> net;
};
