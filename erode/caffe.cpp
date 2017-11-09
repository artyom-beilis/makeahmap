#include <caffe/caffe.hpp>

#include <fstream>
#include <chrono>

using namespace caffe;


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
    
    Caffe::set_mode(Caffe::GPU);
    #ifdef USE_OCL
    Caffe::SetDevice(0);
    Net<float> net(model_file, TEST, Caffe::GetDefaultDevice());
    #else
    Caffe::SetDevice(0);
    Net<float> net(model_file, TEST);
    #endif
    //net.CopyTrainedLayersFrom(trained_file);
    update_network(net,load_universal_model(trained_file));

    Blob<float> *input  = net.input_blobs()[0];

    int chan = input->channels();
    int width = input->width();
    int height = input->height();
    
    
    float* input_data = input->mutable_cpu_data();
    input->Reshape(files.size(),chan,height,width);

    for(size_t i=0;i<files.size();i++) {
	    std::vector<float> f=read_pgm(files[i]);
        size_t exp_size = chan*width*height;
	    if(f.size()!=exp_size)
		    throw std::runtime_error("Invalid file size for " + files[i]);
	    memcpy(input_data + i * exp_size,f.data(),exp_size*sizeof(float));
    }

    auto start = std::chrono::high_resolution_clock::now();
    for(int i=0;i<1000;i++) {
	    input->Reshape(files.size() /*- i%2*/,chan,height,width);
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
    int classes = output->width()  * output->height() * output->channels();
    int total = output->num();
    for(int i=0;i<total;i++) {
    	float mv = 0;
    	int mp = 0;
    	for(int j=0;j<classes;j++) {
    	    if(output_data[i*classes+j] > mv) {
    		    mp = j;
    		    mv = output_data[i*classes+j];
    	    }
	}
		
        printf("%2d: %2d %8.6f\n",i,mp,mv);
    }
    /*std::vector<boost::shared_ptr<Layer<float> > > const &layers = net.layers();
    for(size_t i=0;i<layers.size();i++) {
        std::cout << "Layer " << i << " " << net.layer_names()[i] << std::endl;
        std::vector<boost::shared_ptr<Blob<float> > > const &blobs = layers[i]->blobs();
        for(size_t j=0;j<blobs.size();j++) {
            std::cout << "   " << j << ":" << blobs[j]->shape_string() << std::endl;
        }
    }*/

    std::cout << "Time " << time << " ms total " << (time / files.size()) << " per icon " << std::endl;

}
