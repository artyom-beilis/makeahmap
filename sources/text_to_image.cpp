#include "text_to_image.h"
#include "font.h"
#include <fstream>
#include <stdexcept>
#include <math.h>
#include <iostream>

void render_text::load(std::string const &csv)
{
	labels_.clear();
	std::ifstream f(csv.c_str());
	if(!f) 
		throw std::runtime_error("Failed to open " + csv);
	char bom[4]={};
	f.read(bom,3);
	if(bom==std::string("\xEF\xBB\xBF")) {
		is_utf8_ = true;
	}
	else
	{
		f.seekg(0);
		is_utf8_ = false;
	}
	int line_no = 0;
	std::string line;
	while(std::getline(f,line)) {
		line_no++;
		handle_line(line,csv,line_no);
	}
}

void render_text::get_metrix(render_text::label &l)
{
	l.text = l.src_text;
	l.height = FONT_HEIGHT;
	l.width  = get_print_str_len(l.text.c_str());
}

std::vector<std::string> render_text::split(std::string const &s)
{
	std::vector<std::string> r(1);
	bool inside_string=false;
	char prev = 0;
	for(size_t i=0;i<s.size();i++) {
		switch(s[i]) {
		case '"':
			inside_string = !inside_string;
			if(prev == '"')
				r.back() += '"';
			break;
		case ',':
			if(inside_string) {
				r.back() += ',';
			}
			else {
				r.push_back(std::string());
			}
			break;
		default:
			r.back()+= s[i];
		}
		prev = s[i];
	}
	return r;
}

void render_text::render(std::vector<std::vector<bool> > &grid)
{
	for(size_t i=0;i<labels_.size();i++) {
		label &l=labels_[i];
		print_str(l.text.c_str(),l.row + l.dy,l.col + l.dx,1,grid);
		grid[l.row][l.col]=true;
	}
}

void render_text::handle_line(std::string const &line,std::string const &file,int no)
{
	std::vector<std::string> fields = split(line);
	if(no == 1)
		return;
	if(fields.size() == 1 && fields[0].find_first_not_of(" \t\r\n")==std::string::npos)
		return;
	if(fields.size() < 3) {
		throw std::runtime_error("Expected at least 3 fields in line [" + line + "]" + std::to_string(no));
	}
	label l;
	l.src_text = fields[0];
	l.lat = atof(fields[1].c_str());
	l.lon = atof(fields[2].c_str());
	get_metrix(l);
	l.row = round(height_ * (l.lat - lat2_) / (lat1_ - lat2_));
	l.col = round(width_  * (l.lon - lon1_) / (lon2_ - lon1_));
	if(l.row < 0 || l.col < 0 || l.row >= height_ || l.col >= width_) {
		std::cerr << "   WARNING: "  << ("Coordinates " + std::to_string(l.lat) + ", "  + std::to_string(l.lon) + " for " + l.src_text + " outsize the lat/lon box in line " + std::to_string(no) + " file " + file) << std::endl;
	}
	l.dx = 0;
	l.dy = 0;
	labels_.push_back(l);
}


