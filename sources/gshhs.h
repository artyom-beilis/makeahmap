/*
 * Copyright (c) 2013 Artyom Beilis
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#pragma once 
#include "surface.h"
#include <chrono>
#include "bmp.h"
#include "fileio.h"
#include "downloader.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include <string.h>
#include <stdexcept>
#include <math.h>
#include <functional>
#include <assert.h>
#include <map>
#include <queue>

extern std::string output_dir;

struct water_properties {
    double lat_shift;
    double lon_shift;
    int max_level;
    int default_width;
    std::map<int,int> level_to_width;
    water_properties() :
        lat_shift(0.0),
        lon_shift(0.0),
        max_level(3), // major rivers only
        default_width(1900),
        level_to_width({ { 1, 2500 }, { 2, 2200  }, { 3 , 1900 }  })
    {
    }
};


inline void switch_endian(int &x)
{
	union {
		char buf[4];
		int val;
	} tmp;
	tmp.val = x;
	std::swap(tmp.buf[0],tmp.buf[3]);
	std::swap(tmp.buf[1],tmp.buf[2]);
	x=tmp.val;
}



// global data

template<typename T>
T range_limit(T x,T min,T max)
{
	if(x < min)
		return min;
	if(x > max)
		return max;
	return x;
}

class water_generator {
	struct GSHHS {  /* Global Self-consistent Hierarchical High-resolution Shorelines */
		int id;         /* Unique polygon id number, starting at 0 */
		int n;          /* Number of points in this polygon */
		int flag;       /* = level + version << 8 + greenwich << 16 + source << 24 + river << 25 */
		/* flag contains 5 items, as follows:
		 * low byte:    level = flag & 255: Values: 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake
		 * 2nd byte:    version = (flag >> 8) & 255: Values: Should be 7 for GSHHS release 7 (i.e., version 2.0)
		 * 3rd byte:    greenwich = (flag >> 16) & 1: Values: Greenwich is 1 if Greenwich is crossed
		 * 4th byte:    source = (flag >> 24) & 1: Values: 0 = CIA WDBII, 1 = WVS
		 * 4th byte:    river = (flag >> 25) & 1: Values: 0 = not set, 1 = river-lake and level = 2
		 */
		int west, east, south, north;   /* min/max extent in micro-degrees */
		int area;       /* Area of polygon in 1/10 km^2 */
		int area_full;  /* Area of original full-resolution polygon in 1/10 km^2 */
		int container;  /* Id of container polygon that encloses this polygon (-1 if none) */
		int ancestor;   /* Id of ancestor polygon in the full resolution set that was the source of this polygon (-1 if none) */

		void endian()
		{
			switch_endian(id);
			switch_endian(n);
			switch_endian(flag);
			switch_endian(west);
			switch_endian(east);
			switch_endian(south);
			switch_endian(north);
			switch_endian(area);
			switch_endian(area_full);
			switch_endian(container);
			switch_endian(ancestor);
		}
	};

	struct point {
		int x,y;
		point(int xi=0,int yi=0) : x(xi),y(yi) {}
		void endian()
		{
			switch_endian(x);
			switch_endian(y);
		}
	};
    

public:
	
	static const int sea_mark = 0;
	static const int land_mark = 1;
	static const int lake_mark = 2;
	static const int river_mark = 3;
    int pixel(int r,int c)
    {
        r=(water_size - 1) - r ;
        return internal_pixel(r,c);
    }
    
	bool mark_lake_river_as_river;

	int water_size;
	std::vector<unsigned char> watermap;
	std::vector<int> land_ids;
    std::vector<std::vector<unsigned char> > water_types;
	
	water_generator(double flat1,double flat2,double flon1,double flon2,int water_size_in)
	{
		mark_lake_river_as_river = false;
		lat1=int(round(1e6*flat1));
		lat2=int(round(1e6*flat2));
		lon1=int(round(1e6*flon1));
		lon2=int(round(1e6*flon2));
		water_size = water_size_in;
		watermap.resize(water_size*water_size / 4,0);
		land_ids.resize(water_size*water_size,-1);
		

		if(lat1 > lat2)
			std::swap(lat1,lat2);
		if(lon1 > lon2)
			std::swap(lon1,lon2);

		col_2_lon = double(lon2 - lon1) / water_size;
		lon_2_col = 1.0/col_2_lon;
		row_2_lat = double(lat2 - lat1) / water_size;
		lat_2_row = 1.0/ row_2_lat;
	}
	~water_generator()
	{
		close();
	}

	struct segment {
        double r0,r1,c0,c1;
		double length;
		double n_x,n_y;
        double width;
		void calc_params()
		{
			length = sqrt((r0 - r1)*(r0 - r1) + (c0-c1)*(c0-c1));
			if(length < 0.000001) {
				n_x = 1;
				n_y = 0;
			}
			else {
				n_x = (c1 - c0) / length;
				n_y = (r1 - r0) / length;
			}
			
		}
		double distance(double r,double c) const
		{
			double x0=c,y0=r;
			double p_y=r0,p_x=c0;
			double dx = x0 - p_x;
			double dy = y0 - p_y;
			double proj = dx * n_x + dy * n_y;
			if(proj < 0) 
				return sqrt((dx*dx) + (dy*dy));
			if(proj > length) {
				double dr = r1-r;
				double dc = c1-c;
				return sqrt(dr*dr + dc*dc);
			}
			double ref_x = p_x + n_x * proj - x0;
			double ref_y = p_y + n_y * proj - y0;
			double proj_norm = ref_x * n_y - ref_y * n_x;
			return fabs(proj_norm);
		}
	};

	void load_land(std::string file_name,std::vector<std::vector<int16_t> > &elev,double area_limit_sq_m,int alt_limit)
	{
		open(file_name);
		int area_limit_points = static_cast<int>(area_limit_sq_m * 8 * 8);
		while(get()) {
			int level = (hdr.flag & 255);
			int river = (hdr.flag >> 25) & 1; 

			int value;
			if((level & 1) == 1)
				value = land_mark;
			else
				value = lake_mark;
            
			if(river && mark_lake_river_as_river)
				value = river_mark;
			
			int patch_rows = rmax - rmin+1;
			
			// prepare index
			int total_lines = 0;
			std::vector<std::vector<double> > intersection_points(patch_rows);
			
			assert(poly.front().x == poly.back().x);
			assert(poly.front().y == poly.back().y);

			std::list<segment> cur_segments;
			for(int i = 0;i < points - 1;i++) {
				point start = poly[i];
				point end = poly[i+1];
					
				
				int bot = std::min(start.y,end.y);
				int top = std::max(start.y,end.y);
				
				
				segment s;
                s.width = 0;
				s.r0 = lat_2_row * (start.y - lat1);
				s.r1 = lat_2_row * (end.y   - lat1);
				s.c0 = lon_2_col * (start.x - lon1);
				s.c1 = lon_2_col * (end.x   - lon1);
				s.calc_params();
				cur_segments.push_back(s);

				int row_start  = int(floor(lat_2_row * (bot - lat1))) - 1;
				int row_end   = int(ceil (lat_2_row * (top - lat1))) + 1;	
				
				row_start = range_limit(row_start,rmin,rmax);
				row_end   = range_limit(row_end  ,rmin,rmax);
				
				for(int r=row_start;r<=row_end;r++) {
					double y = (r * row_2_lat) + lat1;
					double x=0;
					if(get_intersection_point(r,y,start,end,x)) {
						intersection_points[r-rmin].push_back(x);
						total_lines ++;
					}
				}
			}
			for(int i=0;i<patch_rows;i++) {
				if(intersection_points[i].size() % 2 != 0) {
					printf("Row case = %d size=%d\n",i+rmin,int(intersection_points[i].size()));
					for(unsigned j=0;j<intersection_points[i].size();j++)
						printf("   x=%f\n",intersection_points[i][j]);
				}
				assert(intersection_points[i].size() % 2 == 0);
				std::sort(intersection_points[i].begin(),intersection_points[i].end());
			}


			// make sets of lines
			std::vector<std::pair<int,std::pair<int,int> > > lines;
			lines.reserve(total_lines);
			for(int r=rmin,row_index=0;r<=rmax;r++,row_index++) {
				std::vector<double> &index = intersection_points[row_index];
				for(size_t i=0;i<index.size();i+=2) {
					int x1 = index[i];
					int x2 = index[i+1];
					if(x1 > x2)
						std::swap(x1,x2);
					int c1 = int(ceil(lon_2_col * (x1 - lon1)));
					int c2 = int(floor(lon_2_col * (x2 - lon1)));
					c1=std::max(0,std::min(c1,water_size-1));
					c2=std::max(0,std::min(c2,water_size-1));
					
					assert(c2-c1>=-1);
					if(c1<=c2)
						lines.push_back(std::make_pair(r,std::make_pair(c1,c2)));
					assert(0  <= c1);
					assert(c2 < water_size);
				}
			}

			bool write_down = true;
			double average = 0.0;
			int total = 0;
			for(size_t i=0;i<lines.size();i++) {
				int r=water_size - 1 - lines[i].first;
				int c1 = lines[i].second.first;
				int c2 = lines[i].second.second;
				for(int c=c1;c<=c2;c++) {
					average+= elev[r][c];
				}
				total += c2 - c1 + 1;
			}
			if(total > 0) {
				average /= total;
				if(value == lake_mark && average > alt_limit)
					write_down = false;
				if(total < area_limit_points) 
					write_down = false;
			}
			else
				write_down = false;
			assert(total >= 0);
			if(write_down) {
				segments.splice(segments.end(),cur_segments);
				for(size_t i=0;i<lines.size();i++) {
					int r=lines[i].first;
					int c1 = lines[i].second.first;
					int c2 = lines[i].second.second;
					for(int c=c1;c<=c2;c++) {
						mark(value,r,c,-1,hdr.id);
					}
				}
			}

		}
		close();
		double max_len = 0;
		for(auto s : segments) {
			max_len = std::max(max_len,s.length);
		}
		std::cout << "-- Total segments:" << segments.size() << ", maximal segment lenght:" << max_len << std::endl;
	}
#if 1
	void write_segments_and_alt(std::vector<std::vector<int16_t> > &elev)
	{
		static const int factor = 8;
		int map_size = water_size * factor;
		std::vector<std::vector<unsigned char> > map(map_size,std::vector<unsigned char>(map_size,0));
		for(int r=0;r<water_size-1;r++) {
			for(int c=0;c<water_size-1;c++) {
				for(int dr=0;dr<factor;dr++) {
					for(int dc=0;dc<factor;dc++) {
						int v00 = elev[water_size-r-1][c];
						int v01 = elev[water_size-r-1][c+1];
						int v10 = elev[water_size-r-2][c];
						int v11 = elev[water_size-r-2][c+1];
						int v0 = ((factor - dc) * v00 + dc * v01 + factor/2) / factor;
						int v1 = ((factor - dc) * v10 + dc * v11 + factor/2) / factor;
						int v  = ((factor - dr) * v0  + dr * v1  + factor/2) / factor;
						map[r*factor+dr][c*factor+dc]=v + 128;
					}
				}
				if(internal_pixel(r,c) == land_mark)
					map[r*factor][c*factor]=255;
				else
					map[r*factor][c*factor]=0;
			}
		}
		for(segment const &s : segments) {
			double y0 = s.r0 * factor;
			double y1 = s.r1 * factor;
			double x0 = s.c0 * factor;
			double x1 = s.c1 * factor;
			if(fabs(x0-x1) > fabs(y0-y1)) {
				double A = (y1-y0) / (x1-x0);
				double B = y0 - x0 * A;
				if(x0 > x1)
					std::swap(x0,x1);
				int xs=round(x0);
				int xe=round(x1);
				assert(xs<=xe);
				for(int x=xs;x<=xe;x++) {
					int y=round(A * x + B);
					if(x>=0 && x<map_size && y>=0 && y<map_size)
						map[y][x]=255;
				}
			}
			else {
				double A = (x1-x0) / (y1-y0);
				double B = x0 - y0*A;
				if(y0 > y1)
					std::swap(y0,y1);
				for(int y=round(y0);y<=round(y1);y++) {
					int x=round(A * y + B);
					if(x>=0 && x<map_size && y>=0 && y<map_size)
						map[y][x]=255;
				}
				
			}
		}
		std::ofstream f((output_dir + "/elevation_mapping.pgm").c_str(),std::ios::binary);
		f<<"P5 " << map_size << " " << map_size << " 255\n";
		for(int r=map_size-1;r>=0;r--)
			f.write(reinterpret_cast<char *>(&map[r][0]),map_size);
	}


    void update_elevations(std::vector<std::vector<int16_t> > &elev,double slope,int alt_limit,surface_solver_options const &opt)
    {
        std::vector<std::vector<int16_t> > calc_alt(water_size,std::vector<int16_t>(water_size));
	    std::vector<std::vector<float> > prev_values;
        std::cout << "-- Converting water bodies to elevations map\n";
        std::cout << "--- Iceles and Lakes\n";
        calc_elevation_boundaries(calc_alt,slope,alt_limit);
        std::cout << "--- Rivers\n";
        calc_river_boundaries(calc_alt,slope,alt_limit);
        save_water_elev(calc_alt,alt_limit);
        std::cout << "-- Calculating elevation modifications to bring water to 0 feet\n";
        std::unique_ptr<surface_solver_base> solver = get_solver(opt);
        std::cout << "   Using:" << solver->name() << "\n\n";
        std::cout << "    " << std::setw(9) << "Grid" << "|" << std::setw(10) << "Vertices"<< "|" << std::setw(9) << "Time (s)" << "|" << std::setw(10) << "Iterations" << "| Bandwidth GB/s" << std::endl;
        double total = 0;
        for(int N=opt.initial_grid;N<=water_size;N*=2) {
            std::vector<std::vector<char> > bmask(N,std::vector<char>(N));
            std::vector<std::vector<float> > bvalues(N,std::vector<float>(N));
            int factor = water_size / N;
            int vertices = (N-1)*4;
            for(int rc=0;rc<N;rc++) {
                bmask[rc][0]=bmask[rc][N-1]=bmask[0][rc]=bmask[N-1][rc]=1;
            }
            for(int r=1;r<N-1;r++) {
                for(int c=1;c<N-1;c++) {
                    if(calc_alt[r*factor][c*factor]<alt_limit) {
                        bmask[r][c]=1;
                        bvalues[r][c]=elev[r*factor][c*factor];
                    }
                    else {
                        vertices++;
                        bmask[r][c]=0;
                        float val = 0;
                        if(!prev_values.empty()) {
                            if((r&1)==0 && (c&1)==0)
                                val=prev_values[r/2][c/2];
                            else if((r&1) == 0) 
                                val=(prev_values[r/2][c/2] + prev_values[r/2][c/2+1])/2;
                            else if((c&1) == 0) 
                                val=(prev_values[r/2][c/2] + prev_values[r/2+1][c/2])/2;
                            else
                                val=(prev_values[r/2][c/2] + prev_values[r/2][c/2+1] +
                                        prev_values[r/2][c/2] + prev_values[r/2+1][c/2])/4;
                        }
                        bvalues[r][c]=val;
                    }
                }
            }
            std::cout<< "    " << std::setw(9) << (std::to_string(N) + "x" + std::to_string(N))  << "|" << std::setw(10) << vertices << "|" << std::flush;
            auto stat = solver->run(bmask,bvalues,opt.threshold);
            total+=stat.time;
            std::cout << std::setw(9) << stat.time << "|"<< std::setw(10) << stat.iterations << "|" << std::setw(9) <<  stat.bandwidth << std::endl;
            prev_values.swap(bvalues);
        }
        std::cout << "\n    Total optimization time " << total << " s"<< std::endl;
        float minv=0;
        float maxv=0;
        for(int r=0;r<water_size;r++) {
            for(int c=0;c<water_size;c++) {
                minv=std::min(minv,prev_values[r][c]);
                maxv=std::max(maxv,prev_values[r][c]);
            }
        }
        std::cout << "    Maximal correction  range max=" << maxv << " min="<<minv << std::endl;
        std::vector<std::vector<int16_t> > ndiff(water_size,std::vector<int16_t>(water_size));
        std::ofstream mem((output_dir + "/sea_level.pgm").c_str(),std::ios::binary);
        std::ofstream neg((output_dir + "/flatten_terrain.pgm").c_str(),std::ios::binary);
		mem<<"P5 " << water_size << " " << water_size << " 255\n";
		neg<<"P5 " << water_size << " " << water_size << " 255\n";
        int16_t max_diff = 0;
        for(int r=0;r<water_size;r++) {
            for(int c=0;c<water_size;c++) {
                int orig_alt = elev[r][c];
                if(calc_alt[r][c]<alt_limit) {
                    elev[r][c]=calc_alt[r][c];
                }
                else {
                    elev[r][c]=elev[r][c]-prev_values[r][c]+alt_limit;
                    if(elev[r][c]<1) {
                        if(orig_alt >= 0) {
                            ndiff[r][c] = alt_limit - elev[r][c];
                            max_diff=std::max(ndiff[r][c],max_diff);
                        }
                        elev[r][c]=1;
                    }
                }
                float v=(prev_values[r][c]-minv) / (maxv-minv) * 254;
                unsigned char color=static_cast<unsigned char>(v);
                mem<<color;
            }
        }
        std::cout << "    Maximal land modificationt to below 0 alt " << max_diff << std::endl;
        for(int r=0;r<water_size;r++)
            for(int c=0;c<water_size;c++)
                neg << static_cast<unsigned char>(ndiff[r][c]*255/max_diff);
    }
    struct bounding_box {
        int rmin,rmax;
        int cmin,cmax;
    };
    bounding_box get_bbox_from_segment(segment const &s,int max_dist)
    {
        bounding_box bb;
        bb.rmin=std::max(0,int(std::min(s.r0,s.r1) - max_dist));
        bb.cmin=std::max(0,int(std::min(s.c0,s.c1) - max_dist));
        bb.rmax=std::min(water_size-1,int(std::max(s.r0,s.r1) + max_dist));
        bb.cmax=std::min(water_size-1,int(std::max(s.c0,s.c1) + max_dist));
        return bb;
    }
    void calc_river_boundaries(std::vector<std::vector<int16_t> > &elev,double slope,int alt_limit)
    {
        double factor = 660 * slope;
        int max_dist = int((alt_limit / slope / 600) + 1);
        for(segment const &s : rivers_) {
            float mdepth = - s.width * slope / 2;
            bounding_box bb=get_bbox_from_segment(s,max_dist);
			for(int r=bb.rmin;r<=bb.rmax;r++) {
				for(int c=bb.cmin;c<=bb.cmax;c++) {
					int depth = int(round(s.distance(r,c) * factor + mdepth));
                    depth = std::max(-alt_limit,std::min(depth,alt_limit));
                    elev[water_size - 1 - r][c] = std::min<int16_t>(elev[water_size - 1 - r][c],depth);
				}
			}
        }
    }
    void save_water_elev(std::vector<std::vector<int16_t> > &elev,int alt_limit)
    {
        std::ofstream f((output_dir + "/water_bodies.pgm").c_str(),std::ios::binary);
        f<<"P5 " << water_size << " " << water_size << " 255\n";
        for(int r=0;r<water_size;r++) {
            for(int c=0;c<water_size;c++) {
                float alt = elev[r][c];
                int ialt = alt/alt_limit * 127 + 127;
                unsigned char calt = std::max(0,std::min(255,ialt));
                f<<calt;
            }
        }
    }
    void calc_elevation_boundaries(std::vector<std::vector<int16_t> > &elev,double slope,int alt_limit)
    {
		double alt_for_cell = slope * 660;
		double max_dist = alt_limit / alt_for_cell + 1;
		std::vector<std::vector<double> > dist(	water_size,std::vector<double>(water_size,max_dist));
		for(segment const &s: segments) {
            bounding_box bb=get_bbox_from_segment(s,max_dist);
			for(int r=bb.rmin;r<=bb.rmax;r++) {
				for(int c=bb.cmin;c<=bb.cmax;c++) {
					double d = s.distance(r,c);
					if(d < dist[r][c])
						dist[r][c] = d;
				}
			}
		}
		for(int r=0;r<water_size;r++) {
			int elev_r = water_size - r - 1;
			for(int c=0;c<water_size;c++) {
				int pos = water_size * r + c;
				int type = (watermap[pos / 4] >> ((pos % 4)*2)) & 0x3;
				int sig = -1;
				if(type == land_mark)  {
					sig = 1;
				}
				int idist = int(round(dist[r][c] * alt_for_cell * sig));
				idist = std::max(-alt_limit,std::min(alt_limit,idist));
				elev[elev_r][c] = idist;
			}
		}
    	    	//write_segments_and_alt(elev);
		return;
    }

#else
	
    void update_elevations(std::vector<std::vector<int16_t> > &elev,double slope)
    {
        water_types.clear();
        water_types.resize(water_size,std::vector<unsigned char>(water_size,0));
		std::set<int> review;
        int pos;
		
		for(int r=0;r<water_size;r++) {
            pos = ((water_size - 1)- r)  * water_size;
			for(int c=0;c<water_size;c++,pos++)  {
                int type = (watermap[pos / 4] >> ((pos % 4)*2)) & 0x3;
                int type_c = c;
                int type_r = (water_size - 1 - r);
                water_types[r][c] |= 1u << (type*2);
				if(r==0 || c==0 || r==water_size-1 || c==water_size-1)
					type = sea_mark;
				switch(type) {
                case land_mark:
					for(int dr=-1;dr<=1;dr++) {
						for(int dc=-1;dc<=1;dc++) {
							int mr = std::max(0,std::min(dr + r,water_size-1));
							int mc = std::max(0,std::min(dc + c,water_size-1));
							int type_mr = (water_size - 1 - mr);
							int tmp_pos = type_mr * water_size + mc;
							int tmp_value = internal_pixel(type_mr,mc);
							if(tmp_value == land_mark && land_ids[pos] != land_ids[tmp_pos] && land_ids[tmp_pos]!=-1) {
								elev[mr][mc] = -10;
								review.insert(mr * water_size + mc);
							}
						}
					}
                    break;
                case sea_mark:
                case lake_mark:
                case river_mark:
					{
						elev[r][c] = -10;
						review.insert(r*water_size + c);
					}
                    break;
				}
			}
		}
        int max_fix = 0;
        static const int horizontal_feet = 660;
        int direct_feet_limit = static_cast<int>(slope * horizontal_feet);
        int cross_feet_limit = static_cast<int>(slope * 1.41421 * horizontal_feet);
		
		std::queue<int> review_queue;
		for(std::set<int>::iterator rp = review.begin();rp!=review.end();rp++)
			review_queue.push(*rp);
        while(!review_queue.empty()) {
            int pos = review_queue.front();
            review_queue.pop();
            int c_r = pos / water_size;
            int c_c = pos % water_size;
            int alt = elev[c_r][c_c];
            for(int dr = -1;dr<=1;dr++) {
                for(int dc=-1;dc<=1;dc++) {
                    if(dr==0 && dc==0)
                        continue;
                    int r=c_r + dr;
                    int c=c_c + dc;
                    if(r < 0 || r>=water_size || c< 0 || c>=water_size)
                        continue;
                    int limit = dr!=0 && dc!=0 ? cross_feet_limit : direct_feet_limit;
                    if(elev[r][c] > alt + limit) {
                        int diff = elev[r][c] - (alt + limit);
                        max_fix = std::max(max_fix,diff);
                        elev[r][c] = alt + limit;
                        review_queue.push(r * water_size + c);
                    }
                }
            }
        }
        std::cout << "  Maximal correction is " << max_fix << std::endl;
    }
#endif	
    
    typedef std::function<void(int/*id*/,int /*segment*/,int/*r*/,int/*c*/)> point_callback_type;
    
    struct callback_guard {
        callback_guard(point_callback_type const &cb,water_generator &gen) : 
            gen_(&gen)
        {
            gen_->cb_ = cb;
        }
        ~callback_guard()
        {
            gen_->cb_ = point_callback_type();
        }
        water_generator *gen_;
    };
    
    
    int alt(std::vector<std::vector<int16_t> > const &elev,double r,double c)
    {
        if(r<=0 || c<=0 || r>=water_size-1 || c>=water_size-1)
            return 0;
        return elev[water_size - int(round(r)) - 1][int(round(c))];
    }

	void load_rivers(std::string file_name,std::vector<std::vector<int16_t> > const &elev,water_properties const &prop,int alt_limit)
	{
        int lat_shift = int(round(prop.lat_shift*1e6));
        int lon_shift = int(round(prop.lon_shift*1e6));
        double max_segment_length = 0.0;
		open(file_name);
		while(get()) {
			int level = (hdr.flag & 255);
            if(prop.max_level != -1 && level > prop.max_level)
                continue;
            int end_point = 0;
			for(int i=0;i<points;i++) {
                poly[i].y+=lat_shift;
                poly[i].x+=lon_shift;
            }
            double width = prop.default_width;
            auto specific_width = prop.level_to_width.find(level);
            if(specific_width != prop.level_to_width.end())
                width = specific_width->second;
			for(int i=points-2;i>=end_point;i--) {
				point start = poly[i+1];
				point end = poly[i];
				if(!(between(start.x,lon1,lon2) && between(start.y,lat1,lat2))
				   && !(between(end.x,lon1,lon2) && between(end.y,lat1,lat2)))
				{
					continue;
				}

						
				segment s;
				s.r0 = lat_2_row * (start.y - lat1);
				s.r1 = lat_2_row * (end.y   - lat1);
				s.c0 = lon_2_col * (start.x - lon1);
				s.c1 = lon_2_col * (end.x   - lon1);
				s.calc_params();
                
                int segment_alt = std::max(alt(elev,s.r0,s.c0),alt(elev,s.r1,s.c1));
                if(segment_alt > alt_limit)
                    break;
                s.width = width;
                rivers_.push_back(s);
                max_segment_length = std::max(s.length,max_segment_length);
			}
		}
        std::cout<< "--- Total river segments: " << rivers_.size() << "; max river segment length " << max_segment_length << std::endl;
        close();
	}

    
    void make_border()
    {
        int border_size = 34 * 2;
        int river_remove_limit = 16;
        for(int d=0;d<border_size;d++) {
            int r1=d;
            int r2=(water_size-1)-d;
            for(int c=0;c<water_size;c++) {
                if(d < river_remove_limit) {
                    mark(sea_mark,r1,c);
                    mark(sea_mark,r2,c);
                    mark(sea_mark,c,r1);
                    mark(sea_mark,c,r2);
                }
                else {
                    to_land_if_river(r1,c);
                    to_land_if_river(r2,c);
                    to_land_if_river(c,r1);
                    to_land_if_river(c,r2);
                }
            }
        }
    }
    
private:
   
    point_callback_type cb_;
    std::list<segment> rivers_;

    int internal_pixel(int r,int c)
    {
        int pos = r * water_size + c;
        int rpos = pos / 4;
        int shift = (pos % 4) * 2;
        int mark = (watermap[rpos] >> shift) & 0x3;
        return mark;
    }
    
    void to_land_if_river(int r,int c)
    {
        if(internal_pixel(r,c)==river_mark) {
            mark(land_mark,r,c);
            assert(internal_pixel(r,c)==land_mark);
        }
    }

  	void mark(int type,int r,int c,int segment=-1,int id=-1,int level=-1)
	{
        if(cb_) {
            cb_(hdr.id,segment,(water_size-1)-r,c);
            return;
        }
		int p = r * water_size + c;
        int preal = p / 4;
        int off = (p % 4) * 2;
        unsigned char current = watermap[preal];

        if(type == river_mark) {
            int value = (current >> off) & 0x3;
	    	if(value != land_mark)
    			return;
        }

        current &= ~(3u << off);
        current |=  (type << off);
		watermap[preal] = current;
		land_ids[p] = id;
	}



	bool between(int x,int min,int max)
	{
		return (min<=x && x<=max);
	}


	void draw_point(point p,double radius,int segment)
	{
        int iradius = floor(radius);
        double r_center_f = (p.y - lat1)*lat_2_row;
        double c_center_f = (p.x - lon1)*lon_2_col;
		int r_center = int(round(r_center_f));
		int c_center = int(round(c_center_f));
		for(int r = r_center-iradius;r<=r_center+iradius;r++) {
			for(int c = c_center-iradius;c<=c_center + iradius;c++) {
                double dr = r - r_center_f;
                double dc = c - c_center_f;
                double R2 = dr*dr + dc * dc;
				if(R2 > radius*radius)
					continue;
				if(c < 0 || c >= water_size || r<0 || r>=water_size)
					continue;
				mark(river_mark,r,c,segment);
			}
		}
	}

	bool get_intersection_point(int r,double y,point p1,point p2,double &x)
	{
		y+=0.005;
		int ymin = std::min(p1.y,p2.y);
		int ymax = std::max(p1.y,p2.y);
		if(!(ymin < y && y< ymax))
			return false;
		double y0 = p1.y - y;
		double y1 = p2.y - y;
		double x0 = p1.x;
		double x1 = p2.x;

		// ensure we never exactly on the line/vertix
		//y0 = y0+0.05;
		//y1 = y1+0.05;
		
		double x_fp = (x0 * double (y1-y0) - y0 * double(x1-x0)) / (y1 - y0);

		// fit into integer range (otherwise use some infinity)
		static const double min_range = -370e6;
		static const double max_range =  370e6;

		if(x_fp < min_range)
			x_fp = min_range;
		else if(x_fp > max_range)
			x_fp = max_range;
		
		x = x_fp;
		return true;
	}

	
	void open(std::string file_name)
	{
        downloader::manager::instance().check(file_name,"gshhs");
		f.open(file_name);
	}
	void close()
	{
        f.close();
	}

    int fix_longitude_sign(int x)
    {
        if(x>180*1000*1000)
        {
			x -= 360 * 1000 * 1000;
        }
        return x;
    }

    void skip()
    {
        f.skip(points * sizeof(point));
    }
	bool get()
	{
		for(;;) {
			if(!f.read(&hdr,sizeof(hdr)))
				return false;
			hdr.endian();
			points = hdr.n;
			greenwich = (hdr.flag >> 16) & 1;
            if(hdr.west > hdr.east) { // ignore polygons passing 180E - 180W
                skip();
				continue;
            }
            if(points <=2) {
                skip();
                continue;
            }
			if(hdr.north < lat1 || hdr.south > lat2) {
                skip();
                continue;
            }
            if(greenwich) {
                if(hdr.west > lon2 || hdr.east < lon1) {
                    skip();
                    continue;   
                }
            }
            else {
                int west = fix_longitude_sign(hdr.west);
                int east = fix_longitude_sign(hdr.east);
                if(west > lon2 || east < lon1) {
                    skip();
                    continue;   
                }
            }

			poly.resize(points);
			if(!f.read(&poly[0],sizeof(point)*points)) {
				throw std::runtime_error("Failed to read file - unexpected EOF");
			}
			for(int i=0;i<points;i++) {
				poly[i].endian();
			}
    		if(greenwich) {
				for(int i=0;i<points;i++) {
					if(poly[i].x > hdr.east)
						poly[i].x -= 360 * 1000 * 1000;
				}
			}
            else {
   				for(int i=0;i<points;i++) {
                    poly[i].x = fix_longitude_sign(poly[i].x);
				}
                hdr.west = fix_longitude_sign(hdr.west);
                hdr.east = fix_longitude_sign(hdr.east);
            }
            
			west_col = int(floor(lon_2_col * (hdr.west - lon1)))-1;
			east_col = int(ceil(lon_2_col * (hdr.east - lon1)))-1;
			north_row = int(ceil(lat_2_row * (hdr.north - lat1)))+1;
			south_row = int(floor(lat_2_row * (hdr.south - lat1)))+1;
			rmin = range_limit(south_row,0,water_size - 1);
			rmax = range_limit(north_row,0,water_size - 1);
			cmin = range_limit(west_col,0,water_size - 1);
			cmax = range_limit(east_col,0,water_size - 1);
			return true;
		}
	}

	int points,greenwich;

	int west_col,east_col,north_row,south_row,cmin,cmax,rmin,rmax;

	int lat1,lat2,lon1,lon2;
	double col_2_lon,lon_2_col,row_2_lat,lat_2_row;
	
	std::vector<point> poly;
	
	std::list<segment> segments;
	
	GSHHS hdr;

	fileio f;
};



// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

