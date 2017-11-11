MXNET_ROOT=$HOME/Packages/mxnet/mxnet/
CUDA_ROOT=/home/artik/Packages/cuda/cuda_8/
TENSORRT_ROOT=$HOME/Packages/cuda/TensorRT-2.1.2/
function bld()
{
	ROOT=$1
	DEFS=$2
	g++ $DEFS -O2 -g -Wall -std=c++11 caffe.cpp -I $ROOT/include/ -L $ROOT/lib/ -Wl,-rpath=$ROOT/lib/ -lcaffe -lboost_system -lglog -o $3 || exit 1
}


g++ -DUSE_TENSORRT -O2 -g -Wall -std=c++11 -I $TENSORRT_ROOT/include -L $TENSORRT_ROOT/lib -Wl,-rpath=$TENSORRT_ROOT/lib caffe.cpp -lnvinfer -lnvcaffe_parser -lcudart -o mnist_tensorrt || exit 1
g++ -DUSE_MXNET -O2 -g -Wall -std=c++11 -I $MXNET_ROOT/include -L $MXNET_ROOT/lib -Wl,-rpath=$MXNET_ROOT/lib:$CUDA_ROOT/lib64 caffe.cpp  -lmxnet -o mnist_mxnet|| exit  1
bld /opt/caffe/caffe-ocl -DUSE_OCL mnist_ocl
bld /opt/caffe/caffe-cuda-cudnn/ "" mnist_cudnn
bld /opt/caffe/caffe-cuda/ "" mnist_cuda

