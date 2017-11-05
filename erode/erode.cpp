#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <math.h>
#include <assert.h>

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

std::pair<float,float> get_dir(image_type const &alt,int r,int c)
{
	float mv=alt[r][c];
	bool found = false;
	for(int dr = -1;!found && dr <= 1;dr++) {
		for(int dc = -1;!found && dc<=1;dc++) {
			if(alt[r+dr][c+dc] < mv) {
				found = true;
			}
		}
	}
	if(!found) {
		return std::make_pair(0.0f,0.0f);
	}
	float dx =   (2*alt[r][c-1] + alt[r-1][c-1] + alt[r+1][c-1])
		   - (2*alt[r][c+1] + alt[r-1][c+1] + alt[r+1][c+1]);
	float dy =   (2*alt[r-1][c] + alt[r-1][c-1] + alt[r-1][c+1])
		   - (2*alt[r+1][c] + alt[r+1][c-1] + alt[r+1][c+1]);

	float norm = sqrtf(dx*dx+dy*dy);
	if(norm < 1e-5) {
		float angle = 3.1415926f * 2.0f * rand() / RAND_MAX;
		return std::make_pair(sinf(angle),cosf(angle));
	}
	else
		return std::make_pair(dx / norm,dy / norm);
}


void rain(image_type &alt,int w,int h,int iterations)
{
	float const rain = 0.01f;
	float const delta = 0.001f;
	image_type water = vals(rain,w,h);
	image_type d_water = vals(0,w,h);
	image_type dissolved = vals(0,w,h);
	image_type errode = vals(0,w,h);
	image_type d_dissolved = dissolved;
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
					float dx = move_to.first;
					float dy = move_to.second;
					if(dx == 0.0f && dy==0.0f) {
						evaporated_in_ponds += water[r][c];
						errode[r][c] -= dissolved[r][c];
						//errode[r][c]  = 0;
						d_water[r][c] = -water[r][c];
						d_dissolved[r][c] = -dissolved[r][c];
					}
					else {
						float wx = dx*dx;
						float wy = dy*dy;
						int rnew = dy < 0 ? r + 1 : r - 1;
						int cnew = dx < 0 ? c + 1 : c - 1;
						
						float d_watr     = water[r][c];
						float errode_one = delta * water[r][c];
						float diss       = dissolved[r][c];

						errode     [r][c] += errode_one;
						d_dissolved[r][c] -= diss;
						d_water    [r][c] -= d_watr;

						float d_diss = errode_one + diss;

						d_dissolved[rnew][c] += d_diss * wy;
						d_water    [rnew][c] += d_watr * wy;

						d_dissolved[r][cnew] += d_diss * wx;
						d_water    [r][cnew] += d_watr * wx;
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
		//std::cout << remainder << " ev " << evaporated << " evp " << evaporated_in_ponds<< std::endl;
		if(remainder < delta)
			break;
	}
}

void smooth(image_type &in,int times)
{
	image_type tmp = in;
	int w = tmp.at(0).size();
	int h = tmp.size();
	for(int i=0;i<times;i++) {
		for(int r=1;r<h-1;r++) {
			for(int c=1;c<w-1;c++) {
				tmp[r][c] = 
					( 1*in[r-1][c-1] + 2*in[r-1][c  ] + 1*in[r-1][c+1] +
					  2*in[r  ][c-1] + 4*in[r  ][c  ] + 2*in[r  ][c+1] +
					  1*in[r+1][c-1] + 2*in[r+1][c  ] + 1*in[r+1][c+1] ) / 16;
			}
		}
		tmp.swap(in);
	}
}


bool check(float &a,float d)
{
	a += d;
	if(a < 0) {
		a = 0.0f;
		return false;
	}
	return true;
}

bool update(image_type &alt,float r,float c,float delta)
{
	int ir = r;
	int ic = c;
	float wr = (c-ic);
	float wl = 1.0f - wr;
	float wb = (r-ir);
	float wt = 1.0f - wb;

	assert(ir >= 0);
	assert(ic >= 0);
	assert(ir+1 < int(alt.size()));
	assert(ic+1 < int(alt.at(0).size()));

	bool b1 = check(alt[ir  ][ic  ],delta * wl*wt);
	bool b2 = check(alt[ir  ][ic+1],delta * wr*wt);
	bool b3 = check(alt[ir+1][ic  ],delta * wl*wb);
	bool b4 = check(alt[ir+1][ic+1],delta * wr*wb);
	return b1 && b2 && b3 && b4;
}

void drop(image_type &alt)
{
	int w = alt.at(0).size();
	int h = alt.size();
	float r = rand() % h;
	float c = rand() % w;
	float dissolved = 0;
	static const float remove = 0.001;
	while(r >= 1 && r<h-2 && c>=1 && c<w-2) {
		auto dir = get_dir(alt,int(r+0.5),int(c+0.5));
		float dx = dir.first;
		float dy = dir.second;
		if(dx == 0.0f && dy==0.0f) {
			update(alt,r,c,dissolved);
			return;
		}
		if(!update(alt,r,c,-remove)) {
			return;
		}
		r = r+dy;
		c = c+dx;
		dissolved += remove;
	}
}


int main(int argc,char **argv)
{
	int iterations = 100;
	char const *name="test.pgm";
	if(argc>=2)
		iterations = atoi(argv[1]);
	if(argc>=3)
		name = argv[2];
	image_type in = read_pgm(name);
	//rain(in,in.at(0).size(),in.size(),iterations);
	for(int it=0;it < iterations;it++) {
		for(size_t i=0;i<in.size()*in.size();i++)
			drop(in);
		if(it%5 == 4) 
			smooth(in,1);
		std::cout << it << std::endl;
	}
	write_pgm("res.pgm",in);
}
