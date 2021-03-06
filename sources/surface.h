#pragma once
#include <vector>
#include <memory>


struct surface_solver_options {
	bool allow_cpu;
	bool force_cpu;
	float threshold;
	int initial_grid;
	int platform_id;
	surface_solver_options() :
		allow_cpu(false),
		force_cpu(false),
		threshold(0.5),
		initial_grid(256),
		platform_id(0)
	{
	}
};

class surface_solver_base {
	surface_solver_base(surface_solver_base const &)=delete;
	void operator=(surface_solver_base const &)=delete;
public:
	
	struct stats {
		double bandwidth;
		int iterations;
		double time;
	};
	
	virtual bool is_cpu() = 0;
	virtual std::string name() = 0;
	virtual stats run(std::vector<std::vector<char> > const &bmask,std::vector<std::vector<float> > &bvalues,float thresh) = 0;
	surface_solver_base() {}
	virtual ~surface_solver_base() {}
};

std::unique_ptr<surface_solver_base> get_solver(surface_solver_options const &opt);

