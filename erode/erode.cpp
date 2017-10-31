#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>

typedef std::vector<std::vector<float > > image_type;

image_type vals(float v,int w,int h)
{
	return image_type(h,std::vector<float>(w,v));
}

image_type read_pgm(std::string const &name)
{
	std::ifstream f(name.c_str(),std::fstream::binary);
	int width,height,range;
	std::string header,line,meta;
	if(!std::getline(f,header))
		throw std::runtime_error("Failed to get header of " + name);
	if(header!="P5")
		throw std::runtime_error(name + " is not pgm");
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
	image_type res(height,std::vector<float>(width));
	std::vector<char> row(width);
	for(int i=0;i<height;i++) {
		f.read(&row[0],width);
		for(int j=0;j<width;j++) {
			res[i][j] = (unsigned char)(row[j]) / 255.0f;
		}
	}
	return res;

}

void write_pgm(std::string const &name,image_type const &img)
{
	int height = img.size();
	int width = img.at(0).size();
	std::vector<char> buf(width);
	std::ofstream f(name.c_str(),std::ofstream::binary);
	f<< "P5\n" << width << " " << height << " 255\n";
	for(int i=0;i<height;i++) {
		for(int j=0;j<width;j++)
			buf[j] = int(std::min(1.0f,std::max(0.0f,img[i][j])) * 255 + 0.5);
		f.write(&buf[0],width);
	}
}

std::pair<int,int> get_dir(image_type const &alt,int r,int c)
{
	float mv=alt[r][c];
	int min_r = r;
	int min_c = c;
	for(int dr = -1;dr <= 1;dr++) {
		for(int dc = -1;dc<=1;dc++) {
			if(alt[r+dr][c+dc] < mv) {
				mv = alt[r+dr][c+dc];
				min_r = r+dr;
				min_c = c+dc;
			}
		}
	}
	return std::make_pair(min_r,min_c);
}


void rain(image_type &alt,int w,int h,int iterations)
{
	float const rain = 0.01f;
	image_type water = vals(rain,w,h);
	image_type d_water = vals(0,w,h);
	image_type dissolved = vals(0,w,h);
	image_type errode = vals(0,w,h);
	image_type d_dissolved = dissolved;
	float const delta = 0.01f;
	for(int count = 0;count < iterations;count++) {
		float evaporated = 0;
		float evaporated_in_ponds = 0;

		for(int r=1;r<h-1;r++) {
			for(int c=1;c<w-1;c++) {
				if(alt[r][c] < delta) {
					d_water[r][c]     = -water[r][c];
					evaporated += water[r][c];
					d_dissolved[r][c] = -dissolved[r][c];
					errode[r][c] = 0;
				}
				else {
					auto move_to = get_dir(alt,r,c);
					int rnew = move_to.first;
					int cnew = move_to.second;
					if(move_to.first == r && move_to.second == c) {
						evaporated_in_ponds += water[r][c];
						errode[r][c] -= dissolved[r][c];
						d_water[r][c] = -water[r][c];
					}
					else {
						float errode_one = delta * water[r][c];
						errode[r][c] += errode_one;
						d_dissolved[r   ][c   ] -= dissolved[r][c];
						d_dissolved[rnew][cnew] += errode_one + dissolved[r][c];
						d_water[r   ][c   ] -= water[r][c];
						d_water[rnew][cnew] += water[r][c];
					}
				}
			}
		}
		float remainder = 0.0;
		for(int r=0;r<h;r++) {
			for(int c=0;c<w;c++) {
				water[r][c]     += d_water[r][c] + rain * rand() / RAND_MAX;
				alt[r][c]       -= errode[r][c];
				dissolved[r][c] += d_dissolved[r][c];
				d_water[r][c] = d_dissolved[r][c] = errode[r][c] = 0.0f;
				
				remainder += water[r][c];
			}
		}
		for(int i=0;i<h;i++) {
			water[i][0]=0;
			water[i][w-1]=0;
			water[0][i]=0;
			water[h-1][i]=0;
		}
		std::cout << remainder << " ev " << evaporated << " evp " << evaporated_in_ponds<< std::endl;
		if(remainder < delta)
			break;
	}
}


int main()
{
	image_type in = read_pgm("test.pgm");
	rain(in,in.at(0).size(),in.size(),1000);
	write_pgm("res.pgm",in);
}
