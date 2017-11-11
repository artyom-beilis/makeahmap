#ifdef USE_MXNET
#include "mxnet_predict.h"
#elif defined USE_TENSORRT
#include "tensorrt_predict.h" 
#else
#include "caffe_predict.h"
#endif

#include <fstream>
#include <sstream>
#include <chrono>
#include <map>
#include <iostream>
#include <string.h>

std::vector<float> read_pgm(std::string const &name)
{
	std::ifstream f(name.c_str(),std::fstream::binary);
	int width,height,range;
	std::string header,line,meta;
	if(!std::getline(f,header))
		throw std::runtime_error("Failed to get header of " + name);
	if(header!="P5" && header!="P6")
		throw std::runtime_error(name + " is not pgm/ppm");
	for(;;) {
		if(!std::getline(f,line))
			throw std::runtime_error("Unexpeced EOF in " + name);
		if(!line.empty() && line[0]=='#')
			continue;
		meta += "\n" + line;
		std::istringstream ss(meta);
		ss >> width >> height >> range;
		if(ss)
			break;
	}
    int chan = header == "P5" ? 1 : 3;
	std::vector<float> res(height*width*chan);
	std::vector<char> row(width*chan);
	int p=0;
	for(int i=0;i<height;i++) {
		f.read(&row[0],width*chan);
		for(int j=0;j<width;j++,p++) {
            for(int c=0;c<chan;c++)
    			res[c*width*height + p] = (unsigned char)(row[j*chan + c]) / 255.0f;
		}
	}
	return res;

}


std::map<std::string,std::vector<float> > load_universal_model(std::string const &name)
{
    std::map<std::string,std::vector<float> > result;
    std::ifstream f(name.c_str());
    if(!f)
        throw std::runtime_error("Falied to open " + name);
    std::string line;
    while(std::getline(f,line)) {
        std::istringstream ss(line);
        std::string name;
        int len;
        ss >> name >> len;
        if(!ss)
            break;
        std::vector<float> &tmp = result[name];
        tmp.resize(len);
        for(int i=0;i<len;i++)
            ss >> tmp[i];
        if(!ss)
            throw std::runtime_error("Failed to read line " + name);
    }
    return result;
}

#if 0
void update_network(Net<float> &net,std::map<std::string,std::vector<float> > const &data)
{
    std::vector<boost::shared_ptr<Layer<float> > > const &layers = net.layers();
    for(size_t i=0;i<layers.size();i++) {
        std::string name = net.layer_names()[i];
        std::vector<boost::shared_ptr<Blob<float> > > const &blobs = layers[i]->blobs();
        for(size_t j=0;j<blobs.size();j++) {
            Blob<float> &blob = *blobs[j];
            size_t size=blob.count();
            std::string dname = name;
            if(j==0)
                dname += "_weight";
            else
                dname += "_bias";
            auto p = data.find(dname);
            if(p==data.end())
                throw std::runtime_error("No layer data for " + dname);
            std::vector<float> const &values = p->second;
            if(values.size() != size) {
                throw std::runtime_error("Size of blob " + dname + "=" + std::to_string(size)+ " differs from " + std::to_string(values.size()));
            }
            memcpy(blob.mutable_cpu_data(),values.data(),size*sizeof(float));
            std::cout << "Loaded " << dname << std::endl;
        }
    }
}
#endif

int main(int argc,char **argv)
{
    if(argc < 3+4) {
	    std::cerr << argv[0] << " model batch channels height width trainedmodel [img1 img2 ... ]\n";
	    return 1;
    }
    std::vector<std::string> files;
    int batch = atoi(argv[1]);
    int chan = atoi(argv[2]);
    int height = atoi(argv[3]);
    int width = atoi(argv[4]);
    std::string model_file  = argv[1+4];
    std::string trained_file = argv[2+4];
    for(int i=3+4;i<argc;i++)
	    files.push_back(argv[i]);
    
    predictor::init();

    predictor net(model_file,trained_file,batch,chan,height,width);
    
    net.reshape(files.size());
    float* input_data = net.input();

    for(size_t i=0;i<files.size();i++) {
	    std::vector<float> f=read_pgm(files[i]);
            size_t exp_size = chan*width*height;
	    if(f.size()!=exp_size)
		    throw std::runtime_error("Invalid file size for " + files[i]);
	    memcpy(input_data + i * exp_size,f.data(),exp_size*sizeof(float));
    }

    auto start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<1000;i++) {
	    net.reshape(files.size());
	    net.input();
	    net.forward();
	    net.result();
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end-start).count();

    float const *output_data = net.result();
    int classes = net.classes();
    int total = net.num();
    for(int i=0;i<total;i++) {
    	float mv = 0;
    	int mp = 0;
    	for(int j=0;j<classes;j++) {
    	    if(output_data[i*classes+j] > mv) {
    		    mp = j;
    		    mv = output_data[i*classes+j];
    	    }
	}
		
        printf("RESULT: %2d: %2d %8.6f      :",i,mp,mv);
	for(int j=0;j<classes;j++) {
		printf("%5.3f ",output_data[i*classes+j]);
	}
	printf("\n");
    }
    std::cout << "Time " << time << " ms total " << (time / files.size()) << " per icon " << std::endl;

}
