#define ENABLE_PROF
//#define PRINT_DEVICE
#include <clcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>

#ifdef USE_HALF

#include <immintrin.h>

inline unsigned short float2half(float f)
{
    return _cvtss_sh(f,0);
}

inline float half2float(unsigned short f)
{
    return _cvtsh_ss(f);
}

#endif

std::string get_program(std::string const &file_name)
{
    std::ifstream f(file_name.c_str());
    if(!f)
        throw std::runtime_error("Failed to read file " + file_name);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

int main()
{
    try {
        context_with_program ctx;
        ctx.load(get_program("prg.cl"));

        kernel test(ctx,"test");
        static const int N = 1024;
        float (*A)[N]=(float (*)[N])malloc(sizeof(float)*N*N);
        float (*B)[N]=(float (*)[N])malloc(sizeof(float)*N*N);
        float (*C)[N]=(float (*)[N])malloc(sizeof(float)*N*N);
        for(int i=0;i<N;i++) {
            for(int j=0;j<N;j++) {
                A[i][j]=2*float(i)/N - 0.5;
                B[i][j]=2*float(j)/N - 0.5;
            }
        }
        memory_object<float> clC(ctx,N*N);
        memory_object<float> clA(ctx,A[0],N*N);
        memory_object<float> clB(ctx,B[0],N*N);

        for(int i=0;i<100;i++)
            test(nd_range(N,N/4),N,clA,clB,clC);

        clC.copy_to_host(C[0],N*N);

        /*for(int i=0;i<N;i++) {
            for(int j=0;j<N;j++) {
                std::cout << C[i][j] << " " ;
            }
            std::cout << "\n";
        }*/

    }
    catch(std::exception const &e) {
        std::cerr << e.what() << std::endl;
    }
}


