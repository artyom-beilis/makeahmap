#include <caffe/caffe.hpp>

#include <fstream>
#include <chrono>

using namespace caffe;


static const double factor = 1.0f;

std::vector<float> read_pgm(std::string const &name)
{
	std::ifstream f(name.c_str(),std::fstream::binary);
	int width,height,range;
	std::string header,line,meta;
	if(!std::getline(f,header))
		throw std::runtime_error("Failed to get header of " + name);
	if(header!="P5")
		throw std::runtime_error(name + " is not pgm");
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
	std::vector<float> res(height*width);
	std::vector<char> row(width);
	int p=0;
	for(int i=0;i<height;i++) {
		f.read(&row[0],width);
		for(int j=0;j<width;j++) {
			res[p++] = 1 - (unsigned char)(row[j]) / 255.0f * factor;
		}
	}
	return res;

}



int main(int argc,char **argv)
{
    if(argc < 3) {
	    std::cerr << argv[0] << " model trainedmodel [img1 img2 ... ]\n";
	    return 1;
    }
    std::vector<std::string> files;
    std::string model_file  = argv[1];
    std::string trained_file = argv[2];
    for(int i=3;i<argc;i++)
	    files.push_back(argv[i]);
    
    static const int dim = 28;
  

    Caffe::set_mode(Caffe::GPU);
    #ifdef USE_OCL
    Caffe::SetDevice(0);
    Net<float> net(model_file, TEST, Caffe::GetDefaultDevice());
    #else
    Caffe::SetDevice(0);
    Net<float> net(model_file, TEST);
    #endif
    net.CopyTrainedLayersFrom(trained_file);
    Blob<float> *input  = net.input_blobs()[0];
    
    
    float* input_data = input->mutable_cpu_data();
    input->Reshape(files.size(),1,dim,dim);
    for(size_t i=0;i<files.size();i++) {
	    std::vector<float> f=read_pgm(files[i]);
	    if(f.size()!=dim*dim)
		    throw std::runtime_error("Invalid file size for " + files[i]);
	    memcpy(input_data + i * dim*dim,f.data(),dim*dim*sizeof(float));
    }

    auto start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<1000;i++) {
	    input->Reshape(files.size() /*- i%2*/,1,dim,dim);
	    input->mutable_cpu_data(); 
	    net.Reshape();
	    net.Forward();
    	    Blob<float> *output = net.output_blobs()[0];
    	    output->cpu_data();
    }
    auto end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end-start).count();

    Blob<float> *output = net.output_blobs()[0];
    float const *output_data = output->cpu_data();
    int total = output->width()  * output->height() * output->num() * output->channels();
    for(int i=0;i<total/10;i++) {
	float mv = 0;
	int mp = 0;
	for(int j=0;j<10;j++) {
	    if(output_data[i*10+j] > mv) {
		    mp = j;
		    mv = output_data[i*10+j];
	    }
	}
		
        printf("%2d: %2d %8.6f\n",i,mp,mv);
    }
    std::cout << "Time " << time << " ms total " << (time / files.size()) << " per icon " << std::endl;

}
