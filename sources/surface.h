#pragma once
#include <vector>
#include <memory>
class surface_solver {
	surface_solver(surface_solver const &)=delete;
	void operator=(surface_solver const &)=delete;
public:
	surface_solver();
	~surface_solver();
	std::string name();
	std::pair<int,double> run(std::vector<std::vector<char> > const &bmask,std::vector<std::vector<float> > &bvalues,float thresh);
private:
	struct data;
	std::unique_ptr<data> impl;
		
};

