#include <caffe/caffe.hpp>

#include <fstream>

using namespace caffe;

int main(int argc,char **argv)
{
    std::string model_file  = argv[1];
    std::string trained_file = argv[2];
    static const int dim = 28;
    //Caffe::set_mode(Caffe::GPU);
    //Caffe::SetDevice(0);
    Net<float> net(model_file, TEST, Caffe::GetDefaultDevice());
    net.CopyTrainedLayersFrom(trained_file);
    Blob<float> *input  = net.input_blobs()[0];

    input->Reshape(3,1,dim,dim);

    float* input_data = input->mutable_cpu_data();
    for(int r=0;r<dim;r++) {
        for(int c=0;c<dim;c++) {
            float dx = float(r - dim/2) / dim * 2;
            float dy = float(c - dim/2) / dim * 2;
            if(-0.6 < dx && dx < 0.6 && -0.1 < dy && dy < 0.1)
                *input_data++ = 0.0f;
            else
                *input_data++ = 1.0f;
            std::cout << input_data[-1];
        }
        std::cout <<"\n";
    }

    for(int r=0;r<dim;r++) {
        for(int c=0;c<dim;c++) {
            float dx = float(r - dim/2) / dim * 2;
            float dy = float(c - dim/2) / dim * 2;
            float r2 = dx*dx + dy*dy;
            if(0.6 < r2 && r2 < 0.8)
                *input_data++ = 0.0f;
            else
                *input_data++ = 1.0f;
            std::cout << input_data[-1];
        }
        std::cout <<"\n";
    }

    std::ifstream f("5.data");
    char buf[dim*dim];
    f.read(buf,sizeof(buf));
    for(int i=0;i<dim*dim;i++) {
        *input_data++ = (unsigned char)buf[i] / 255.0f;
        std::cout << (input_data[-1] > 0.5);
        if(i%dim==dim-1)
            std::cout << "\n";
    }

    net.Reshape();
    net.Forward();

    Blob<float> *output = net.output_blobs()[0];
    float const *output_data = output->cpu_data();
    int total = output->width()  * output->height() * output->num() * output->channels();
    for(int i=0;i<total;i++)
        printf("%2d %8.6f\n",i,output_data[i]);

}
