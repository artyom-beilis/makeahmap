#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <chrono>

int main()
{
	size_t size = 1024*1024*1024;
	void *ptr = malloc(size);
	assert(ptr);
	for(int i=0;i<3;i++)
		memset(ptr,i,size);
	int points=30;
	{
		auto start = std::chrono::high_resolution_clock::now();
		for(int i=0;i<points;i++) {
			memset(ptr,i+1,size);
		}
		auto end = std::chrono::high_resolution_clock::now();
		double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end-start).count();
		std::cout << "Write:     "<< std::setw(8) << (size * points / time / (1024.0*1024*1024)) <<" GB/s"<< std::endl;
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		int *ip = (int*)(ptr);
		size_t isize = size/sizeof(*ip);
		int sum=0;
		for(int i=0;i<points;i++) {
			for(size_t j=0;j<isize;j++)
				sum+=ip[j];
		}
		auto end = std::chrono::high_resolution_clock::now();
		double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end-start).count();
		volatile int vsum = sum;
		std::cout << "Read:      "<< std::setw(8) << (size * points / time / (1024.0*1024*1024)) <<" GB/s "<< std::endl;
	}
	{
		auto start = std::chrono::high_resolution_clock::now();
		int *ip = (int*)(ptr);
		size_t isize = size/sizeof(*ip);
		int sum=0;
		for(int i=0;i<points/2;i++) {
			memset(ptr,i+5,size);
			for(size_t j=0;j<isize;j++)
				sum+=ip[j];
		}
		auto end = std::chrono::high_resolution_clock::now();
		double time = std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1> > >(end-start).count();
		volatile int vsum = sum;
		std::cout << "Read/Write:"<< std::setw(8) << (size * points / time / (1024.0*1024*1024)) <<" GB/s "<< std::endl;
	}
	free(ptr);
	return 0;
}
