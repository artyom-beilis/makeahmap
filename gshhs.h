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

#include "bmp.h"
#include "fileio.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <stdint.h>
#include <string.h>
#include <stdexcept>
#include <math.h>
#include <functional>
#include <assert.h>
#include <map>


struct water_properties {
    double lat_shift;
    double lon_shift;
    int max_level;
    int default_width;
    std::map<int,int> level_to_width;
    std::map<int,int> start_points;
    int colors[4];
    int depth_range[4];
    water_properties() :
        lat_shift(0.0),
        lon_shift(0.0),
        max_level(-1),
        default_width(1),
        level_to_width({ { 1, 6 }, { 2, 5 }, { 3 , 4 } , { 4 , 3 }, { 5, 3} , { 5 , 2 } }),
        colors({0,0,0,0}),
        depth_range({0,0,0,0})
    {
        lat_shift = 0;
        lon_shift = 0;
        max_level = -1;
        default_width = 2;
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

	void load_land(std::string file_name)
	{
		open(file_name);
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
			std::vector<std::vector<int> > intersection_points(patch_rows);
			
			assert(poly.front().x == poly.back().x);
			assert(poly.front().y == poly.back().y);

			for(int i = 0;i < points - 1;i++) {
				point start = poly[i];
				point end = poly[i+1];
				
				int bot = std::min(start.y,end.y);
				int top = std::max(start.y,end.y);

				int row_start  = int(floor(lat_2_row * (bot - lat1))) - 1;
				int row_end   = int(ceil (lat_2_row * (top - lat1))) + 1;	
				
				row_start = range_limit(row_start,rmin,rmax);
				row_end   = range_limit(row_end  ,rmin,rmax);
				
				for(int r=row_start;r<=row_end;r++) {
					int y = int(round(r * row_2_lat)) + lat1;
					int x=0;
					if(get_intersection_point(y,start,end,x)) {
						intersection_points[r-rmin].push_back(x);
					}
				}
			}

			for(int i=0;i<patch_rows;i++) {
				assert(intersection_points[i].size() % 2 == 0);
				std::sort(intersection_points[i].begin(),intersection_points[i].end());
			}

			// fill the lines
			for(int r=rmin,row_index=0;r<=rmax;r++,row_index++) {
				std::vector<int> &index = intersection_points[row_index];
				for(size_t i=0;i<index.size();i+=2) {
					int x1 = index[i];
					int x2 = index[i+1];
					int c1 = int(round(lon_2_col * (x1 - lon1)));
					int c2 = int(round(lon_2_col * (x2 - lon1)));
					c1=std::max(0,c1);
					c2=std::min(water_size-1,c2);
					for(int c=c1;c<=c2;c++) {
						mark(value,r,c);
					}
				}
			}

		}
		close();
	}
	
    void save_waterc_map(std::string out,water_properties const &prop)
	{
		bmp::header hdr(16384,16384);
        std::vector<unsigned char> zeros(16384,0);
        int padding=(16384 - water_size)/2;
        outfile f(out);
        f.write(&hdr,sizeof(hdr));

		std::vector<unsigned char> data(water_size);
        

		int pos = 0;
        for(int i=0;i<padding;i++)
            f.write(&zeros[0],16384);
		for(int r=0;r<water_size;r++) {
			for(int c=0;c<water_size;c++,pos++)  {
                int type = (watermap[pos / 4] >> ((pos % 4)*2)) & 0x3;
                unsigned char color = prop.colors[type];
                data[c]=color;
			}
            f.write(&zeros[0],padding);
			f.write(&data[0],data.size());
            f.write(&zeros[0],padding);
		}
        for(int i=0;i<padding;i++)
            f.write(&zeros[0],16384);
        f.close();
	}

	void save_waterd_map(std::string out,water_properties const &prop)
	{
        water_types.clear();
        
        water_types.resize(water_size/4,std::vector<unsigned char>(water_size/4,0));
        
		bmp::header hdr(16384,16384);
        std::vector<unsigned char> zeros(16384,0);
        int padding=(16384 - water_size)/2;
        outfile f(out);
		f.write(&hdr,sizeof(hdr));
        
		std::vector<unsigned char> data(water_size);

        std::vector<int> lut[4];
        for(int i=0;i<4;i++) {
            int l = prop.depth_range[i];
            int l2 = l*l;
            lut[i].resize(l2+1,0);
            for(int j=0;j<l2;j++) {
                int v=int(round(sqrt(l2-j) * 254 / l));
                if(v<0)
                    v=0;
                if(v>254)
                    v=254;
                lut[i][j]=v;
            }
        }

		int pos = 0;
        for(int i=0;i<padding;i++)
            f.write(&zeros[0],16384);
        int prev_percent = -1;
        
		for(int r=0;r<water_size;r++) {
            int percent = (r+1) * 100 / water_size;
            if(percent != prev_percent) {
                if(prev_percent!=-1)
                    std::cout << "\b\b\b\b\b";
                prev_percent = percent;
                std::cout << std::setw(3) << percent << "% " << std::flush;
            }

            bool prev_water_only = false;
            int prev_type = -1;
			for(int c=0;c<water_size;c++,pos++)  {
                int type = internal_pixel(r,c);
                if(type == land_mark) {
                    prev_type = land_mark;
                    prev_water_only = false;
                    data[c] = 255;
                    continue;
                }
                int depth_dist = prop.depth_range[type];
                if(depth_dist == 0) {
                    prev_water_only = false;
                    prev_type = type;
                    data[c] = 0;
                    continue;
                }
                int closest_dist2 = depth_dist * depth_dist;
                
                int r_start = std::max(r-depth_dist,0);
                int r_end  = std::min(r+depth_dist,water_size-1);
                int c_start = std::max(c-depth_dist,0);
                int c_end  = std::min(c+depth_dist,water_size-1);

                if(prev_water_only && prev_type == type) {
                    c_start = c_end;
                }
                bool water_only = true;
                for(int tr = r_start,dr=r_start - r;tr<=r_end;tr++,dr++) {
                    int r2 = dr*dr;
                    for(int tc = c_start,dc = c_start - c;tc<=c_end;tc++,dc++) {
                        int type = internal_pixel(tr,tc);
                        if(type == land_mark) {
                            int c2 = dc*dc;
                            int D = r2 + c2;
                            if(D < closest_dist2)
                                closest_dist2 = D;
                            water_only = false;
                        }
                    }
                }
                prev_water_only = water_only;
                prev_type = type;
                
                data[c] = lut[type][closest_dist2];
			}
            for(int c=0;c<water_size;c++) {
                int type = internal_pixel(r,c);
                int type_c = c/4;
                int type_r = (water_size - 1 - r)/4;
                water_types[type_r][type_c] |= 1u << (type*2);
                if(type != land_mark && data[c]!=0)
                    water_types[type_r][type_c] |= 2u << (type*2);
            }
            f.write(&zeros[0],padding);
			f.write(&data[0],data.size());
            f.write(&zeros[0],padding);
		}
        for(int i=0;i<padding;i++)
            f.write(&zeros[0],16384);
        f.close();
	}
    
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
    
	void load_rivers(std::string file_name,water_properties const &prop)
	{
        int lat_shift = int(round(prop.lat_shift*1e6));
        int lon_shift = int(round(prop.lon_shift*1e6));
		open(file_name);
		while(get()) {
			int level = (hdr.flag & 255);
            if(prop.max_level != -1 && level > prop.max_level)
                continue;
            int end_point = 0;
            auto black_list = prop.start_points.find(hdr.id);
            if(black_list!=prop.start_points.end())
                end_point = black_list->second;
			for(int i=0;i<points;i++) {
                poly[i].y+=lat_shift;
                poly[i].x+=lon_shift;
            }
            double width = prop.default_width;
            auto specific_width = prop.level_to_width.find(level);
            if(specific_width != prop.level_to_width.end())
                width = specific_width->second;
            double radius = width / 2;
			for(int i=points-2;i>=end_point;i--) {
				point start = poly[i+1];
				point end = poly[i];
				if(!(between(start.x,lon1,lon2) && between(start.y,lat1,lat2))
				   && !(between(end.x,lon1,lon2) && between(end.y,lat1,lat2)))
				{
					continue;
				}

				draw_point(start,radius,i);
				draw_point(end,radius,i);

				if(labs(start.x - end.x) > labs(start.y - end.y)) {

					int c_start = int(round((start.x - lon1) * lon_2_col));
					int c_end   = int(round((end.x   - lon1) * lon_2_col));
					
					c_start = range_limit(c_start,cmin,cmax);
					c_end   = range_limit(c_end,cmin,cmax);

					if(c_start > c_end)
						std::swap(c_start,c_end);

					double a=double(end.y - start.y) / double(end.x-start.x);
					double b=start.y - a * start.x;
					
					for(int c=c_start;c<=c_end;c++) {
						double x = col_2_lon * c + lon1;
						double y = a*x + b;
                        draw_point(point(int(round(x)),int(round(y))),radius,i);
                        /*
						int r = int(round((y - lat1) * lat_2_row - half_width));
						for(int rr=r;rr<=r+integer_width;rr++) {
							if(rr >= 0 && rr<water_size) {
								mark(river_mark,rr,c);
							}
						}*/
					}
				}
				else if(start.y != end.y) {
					
					int r_start = int(round((start.y - lat1) * lat_2_row));
					int r_end   = int(round((end.y   - lat1) * lat_2_row));
					
					r_start = range_limit(r_start,rmin,rmax);
					r_end   = range_limit(r_end,rmin,rmax);

					double a=double(end.x - start.x) / double(end.y - start.y);
					double b=start.x - a * start.y;

					if(r_start > r_end)
						std::swap(r_start,r_end);
					
					for(int r=r_start;r<=r_end;r++) {
						double y = row_2_lat * r + lat1;
						double x = a*y + b;
                        
                        draw_point(point(int(round(x)),int(round(y))),radius,i);
                        /*
						int c = int(round((x - lon1) * lon_2_col - half_width));
						if(c<0 || c>=water_size)
							continue;
						
						for(int cc=c-width;cc<=c+width;cc++) {
							if(cc >= 0 && cc<water_size) {
								mark(river_mark,r,cc);
							}
						}*/
					}
				}

			}
		}
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

  	void mark(int type,int r,int c,int segment=-1)
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

	bool get_intersection_point(int y,point p1,point p2,int &x)
	{
		int y0 = p1.y - y;
		int y1 = p2.y - y;
		int x0 = p1.x;
		int x1 = p2.x;

		// ensure we never exactly on the line/vertix
		y0 = y0*2+1;
		y1 = y1*2+1;

		// line does not cross
		if((y0 > 0 && y1 >0) || (y0 < 0 && y1 < 0))
			return false;
		
		double x_fp = (x0 * double (y1-y0) - y0 * double(x1-x0)) / (y1 - y0);

		// fit into integer range (otherwise use some infinity)
		static const double min_range = -370e6;
		static const double max_range =  370e6;

		if(x_fp < min_range)
			x_fp = min_range;
		else if(x_fp > max_range)
			x_fp = max_range;
		
		x = int(round(x_fp));
		return true;
	}

	
	void open(std::string file_name)
	{
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
            
			west_col = int(round(lon_2_col * (hdr.west - lon1)));
			east_col = int(round(lon_2_col * (hdr.east - lon1)));
			north_row = int(round(lat_2_row * (hdr.north - lat1)));
			south_row = int(round(lat_2_row * (hdr.south - lat1)));
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
	GSHHS hdr;

	fileio f;
};



// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

