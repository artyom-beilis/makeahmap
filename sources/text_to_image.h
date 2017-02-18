#pragma once
#include "image.h"
#include <vector>

class render_text {
public:
	render_text(double lat1,double lat2,double lon1,double lon2,int size) :
		lat1_(lat1),
		lat2_(lat2),
		lon1_(lon1),
		lon2_(lon2),
		height_(size),
		width_(size)
		//raster_(size,size)
	{
	}
	void load(std::string const &csv);
	void optimize() {}
	void render(std::vector<std::vector<bool> > &grid);
	
private:
	
	struct label {
		double lat,lon;
		std::string src_text;
		std::string text;
		int height,width;
		int row,col;
		int dx,dy;
	};
	double lat1_,lat2_,lon1_,lon2_;
	int height_;
	int width_;
	
	bool is_utf8_;
	void get_metrix(label &l);
	void handle_line(std::string const &line,std::string const &csv,int no);
	std::vector<std::string> split(std::string const &s);
	std::vector<label> labels_;
};