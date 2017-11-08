function bld()
{
	ROOT=$1
	DEFS=$2
	g++ $DEFS -O2 -g -Wall -std=c++11 caffe.cpp -I $ROOT/include/ -L $ROOT/lib/ -Wl,-rpath=$ROOT/lib/ -lcaffe -lboost_system -lglog -o $3
}

bld /opt/caffe/caffe-ocl -DUSE_OCL mnist_ocl
bld /opt/caffe/caffe-cuda-cudnn/ "" mnist_cudnn
bld /opt/caffe/caffe-cuda/ "" mnist_cuda

