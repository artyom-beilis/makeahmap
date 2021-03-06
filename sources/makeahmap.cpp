/*
 * Copyright (c) 2013-2106 Artyom Beilis
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


#include <tiffio.h>
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
#include <iomanip>
#include <limits>
#include <string>
#include <thread>
#include <streambuf>
#include <set>
#include <functional>
#include "bmp.h"
#include "gshhs.h"
#include "dem.h"
#include "downloader.h"
#include "surface.h"
#include "util.h"
#include "getversion.h"
#include "font.h"
#include "image.h"
#include "text_to_image.h"

std::vector<std::vector<int16_t> > elevations;
int map_size = 512;
int cbm_size = -1;
int river_correction_limit = -1;
int altitude_offset = 0;
bool remove_entire_river = false;
dem::db_properties db_type;

std::vector<std::vector<unsigned char> > bare_types;
//std::vector<std::vector<int> > types;
image<int> types;

static const int lower_lat = -65;

double lat1 = -1000,lat2 = -1000;
double lon1 = -1000,lon2 = -1000;

water_properties rv_prop;

double river_north_shift = 0.0;
double river_east_shift = 0.0;
bool fix_river_slopes = true;
int water_alt_limit = 30000; //higher than everest
double lake_or_island_min_size = 0.5; // sq miles
double water_to_land_slope = 0.1; 
int water_to_land_range = 100; // feet
int splattype_filter_radius = 7;
int check_updates=1;
bool disable_sse=false;

std::string map_name = "map";
std::string tiff_file="./data/globcover/globcover.tif";
std::string shores = "./data/gshhs/gshhs_f.b";
std::string rivers = "./data/gshhs/wdb_rivers_f.b";
std::string grid_dir = "./images";
std::string tile_set = "europe";
std::string tile_set_path = "./resources/groundmapping";
std::string custom_mapping;
std::string output_dir = "./output";
std::string download_sources = "download_sources.txt";
std::string cbm_labels;
std::string temp_dir = "./temp";
std::string cbm_grid_type = "game";
bool auto_download_enabled = true;
bool disable_ssl_check = false;
char const *real_file=0;
char const *color_file=0;
enum altitude_handling_type {
    shade_hills,
    dark_hills,
    untouched_hills
} altitude_handling = shade_hills;
double shading_factor = 1.0;
double map_scale;

surface_solver_options solver_options;

class parsing_error : public std::runtime_error {
public:
    parsing_error(std::string const &v) : std::runtime_error(v) {}
};


class color {
public:
    unsigned char r,g,b;
    color(int ir=0,int ig=0,int ib=0) : r(ir), g(ig), b(ib) {}
    color(std::string const &color) :r(0),g(0),b(0)
    {
        bool ok=false;
        do {
            if(color.empty())
                break;
            if(color[0]=='h') {
                if(color.size()==4) {
                    r=from_hex(color[1])*17;
                    g=from_hex(color[2])*17;
                    b=from_hex(color[3])*17;
                    ok=true;
                }
                else if(color.size()==7) {
                    r=from_hex(color[1])*16+from_hex(color[2]);
                    g=from_hex(color[3])*16+from_hex(color[4]);
                    b=from_hex(color[5])*16+from_hex(color[6]);
                    ok=true;
                }
                else 
                    break;
            }
            else if(color[0]=='(' && color.back()==')') {
                size_t p1=color.find_first_of(',');
                if(p1==std::string::npos) 
                    break;
                size_t p2=color.find_first_of(',',p1+1);
                if(p2==std::string::npos)
                    break;
                std::string sr=color.substr(1,p1-1);
                std::string sg=color.substr(p1+1,p2-(p1+1));
                std::string sb=color.substr(p2+1,color.size()-1-(p2+1));
                r=atoi(sr.c_str());
                g=atoi(sg.c_str());
                b=atoi(sb.c_str());
                ok=true;
            }
        } while(0);
        if(!ok)
            throw parsing_error("Invalid color specification `"+color+"'");
    }
private:
    static int from_hex(char c)
    {
        if('0' <= c && c<= '9')
            return c-'0';
        if('a' <= c && c<= 'f')
            return c-'a'+10;
        if('A' <= c && c<= 'F')
            return c-'A'+10;
        throw parsing_error(std::string("Invalid hexadecimal digit `") + c + "'" );
    }
};

color map_colors[5] = {
        color( 178, 210, 174 ), //   0- 49 - farm 
        color( 165, 200, 163 ), //  50- 99 - forest
        color( 220, 220, 240 ), // 100-149 - snow
        color( 216, 225, 185 ), // 150-199 - grass, grocky grass, sandy grass, swamp
        color( 246, 241, 220 )  // 200-249 - beach, river bad, rock
};

color water_color(167,201,253);
color beach_color(125,151,190);
color grid_color(64,64,64); 



// tiles
// ntt0000  Deep water ntt0000 will always be deep water  (This is a greyscale BMP designed to vary the water over distances)
// ntt0001  Grass
// ntt0002  Forest 1 (evergreen)
// ntt0003  Forest 2 (deciduous)
// ntt0004  Farm 1
// ntt0005  Farm 2
// ntt0006  Rock
// ntt0007  Swamp
// ntt0008  Rocky Grass
// ntt0009  Sandy Grass
// ntt0010 A  Beach
// ntt0011 B  River bed /
// ntt0012 C  Snow / Coral

static const unsigned char beach_type = 0xAA;
/*

AH2 types 

static const int map_type_in[][3] = {
 { 11,  0x44 }, //Post-flooding or irrigated croplands (or aquatic)
 { 14,  0x0080FF }, //Rainfed croplands
 { 20,  0x0040FF }, //Mosaic cropland (50-70%) / vegetation (grassland/shrubland/forest) (20-50%)
 { 30,  0x42 }, //Mosaic vegetation (grassland/shrubland/forest) (50-70%) / cropland (20-50%) 
 { 40,  0x21 }, //Closed to open (>15%) broadleaved evergreen or semi-deciduous forest (>5m)
 { 50,  0xA00000 }, //Closed (>40%) broadleaved deciduous forest (>5m)
 { 60,  0xA00000 }, //Open (15-40%) broadleaved deciduous forest/woodland (>5m)
 { 70,  0x33 }, //Closed (>40%) needleleaved evergreen forest (>5m)
 { 90,  0x31 }, //Open (15-40%) needleleaved deciduous or evergreen forest (>5m)
 { 100, 0x32 }, //Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)
 { 110, 0x11 }, //Mosaic forest or shrubland (50-70%) / grassland (20-50%)
 { 120, 0x12 }, //Mosaic grassland (50-70%) / forest or shrubland (20-50%) 
 { 130, 0x21 }, //Closed to open (>15%) (broadleaved or needleleaved, evergreen or deciduous) shrubland (<5m)
 { 140, 0x21 }, //Closed to open (>15%) herbaceous vegetation (grassland, savannas or lichens/mosses)
 { 150, 0x11 }, //Sparse (<15%) vegetation
 { 160, 0x24 }, //Closed to open (>15%) broadleaved forest regularly flooded (semi-permanently or temporarily) - Fresh or brackish water
 { 170, 0x25 }, //Closed (>40%) broadleaved forest or shrubland permanently flooded - Saline or brackish water
 { 180, 0x24 }, //Closed to open (>15%) grassland or woody vegetation on regularly flooded or waterlogged soil - Fresh, brackish or saline water
 { 190, 0x88 },    // Artificial surfaces and associated areas (Urban areas >50%)
 { 200, 0xAA },    // Bare areas
 { 210, 0xAA },    // Waterbodies
 { 220, 0xCC },    // Permanent snow and ice
 { 230, 0xAA },    // No data (burnt areas, clouds,…)
 { -1, 0}
};
*/
static const int map_type_in[][2] = {
 { 11,  16 }, //Post-flooding or irrigated croplands (or aquatic)
 { 14,  18 }, //Rainfed croplands
 { 20,  17 }, //Mosaic cropland (50-70%) / vegetation (grassland/shrubland/forest) (20-50%)
 { 30,  17 }, //Mosaic vegetation (grassland/shrubland/forest) (50-70%) / cropland (20-50%) 
 { 40,  9 }, //Closed to open (>15%) broadleaved evergreen or semi-deciduous forest (>5m)
 { 50,  10 }, //Closed (>40%) broadleaved deciduous forest (>5m)
 { 60,  8 }, //Open (15-40%) broadleaved deciduous forest/woodland (>5m)
 { 70,  10 }, //Closed (>40%) needleleaved evergreen forest (>5m)
 { 90,  8 }, //Open (15-40%) needleleaved deciduous or evergreen forest (>5m)
 { 100, 9 }, //Closed to open (>15%) mixed broadleaved and needleleaved forest (>5m)
 { 110, 17 }, //Mosaic forest or shrubland (50-70%) / grassland (20-50%)
 { 120, 16 }, //Mosaic grassland (50-70%) / forest or shrubland (20-50%) 
 { 130, 5 }, //Closed to open (>15%) (broadleaved or needleleaved, evergreen or deciduous) shrubland (<5m)
 { 140, 5 }, //Closed to open (>15%) herbaceous vegetation (grassland, savannas or lichens/mosses)
 { 150, 2 }, //Sparse (<15%) vegetation
 { 160, 8 }, //Closed to open (>15%) broadleaved forest regularly flooded (semi-permanently or temporarily) - Fresh or brackish water
 { 170, 8 }, //Closed (>40%) broadleaved forest or shrubland permanently flooded - Saline or brackish water
 { 180, 7 }, //Closed to open (>15%) grassland or woody vegetation on regularly flooded or waterlogged soil - Fresh, brackish or saline water
 { 190, 19 },    // type 19 Artificial surfaces and associated areas (Urban areas >50%)
 { 200, 0 },    // type 0 Bare areas
 { 210, 1 },    // type 1 Waterbodies
 { 220, 15 },    // type 15 Permanent snow and ice
 { 230, 0 },    // type 0 No data (burnt areas, clouds,…)
 { -1, 0}
};



bool is_water(int c) { return c==210 || c==230; }

static int map_type[256];

static const int map_in[][4] = {
    {11,170,240,240},
    {14,255,255,100},
    {20,220,240,100},
    {30,205,205,102},
    {40,0,100,0},
    {50,0,160,0},
    {60,170,200,0},
    {70,0,60,0},
    {90,40,100,0},
    {100,120,130,0},
    {110,140,160,0},
    {120,190,150,0},
    {130,150,100,0},
    {140,255,180,50},
    {150,255,235,175},
    {160,0,120,90},
    {170,0,150,120},
    {180,0,220,130},
    {190,195,20,0},
    {200,255,245,215},
    {210,0,70,200},
    {220,255,255,255},
    {230,0,0,0},
    {-1,0,0,0}
};

static unsigned char tmap[256][3];

unsigned tile_to_color(unsigned type)
{
	if(type < 16u) 
		type = type << (4 + 16);
	else 
		type = (((type - 16u) * 0x40u) << 8) + 0xFF;
	return type;
}
void prepare_map()
{
    for(int i=0;map_in[i][0]!=-1;i++) {
        for(int j=0;j<3;j++) 
            tmap[map_in[i][0]][j]=map_in[i][j+1];
    }
	for(int i=0;map_type_in[i][0]!=-1;i++) {
        map_type[map_type_in[i][0]] = tile_to_color(map_type_in[i][1]);
    }
}


void load_custom_type_mapping()
{
    prepare_map();
	std::string path;
	std::string tag;
    if(custom_mapping.empty()) {
		path=tile_set_path + "/" + tile_set + ".csv";
		tag = tile_set + ".csv";
	}
	else {
		path = custom_mapping;
		tag = custom_mapping;
	}
    std::cout <<"- Loading ground type mapping " << tag << std::endl;
    std::ifstream csv(path.c_str());
    if(!csv) {
        throw std::runtime_error("Failed to open file " + path);
    }
    std::string line;
    int lineno=0;
    while(std::getline(csv,line)) {
        lineno++;
        std::istringstream ss(line);
        ss.unsetf(std::ios_base::basefield);
        unsigned gcover=0,ah=0;
        char c;

        ss >> std::dec >> gcover >> c >> ah;
        if(!ss || c!=',' || gcover >255u || ah > 19u) {
            std::ostringstream tmp;
            tmp << "Invalid line " << lineno << " in file " << path;
            throw std::runtime_error(tmp.str());
        }
        map_type[gcover]=tile_to_color(ah);
    }
}

bool parse(std::string const &in,std::string &key,std::string &value,std::string &index)
{
    key.erase();
    value.erase();
    index.erase();

    size_t key_start = in.find_first_not_of(" \t\r\n");
    if(key_start==std::string::npos || in[key_start]=='#')
        return false;
    size_t key_end = in.find_first_of(" \r\n\t=",key_start);
    if(key_end==std::string::npos)
        throw parsing_error("Key missing value");
    size_t p = in.find_first_not_of(" \r\n\t",key_end);
    if(p==std::string::npos)
        throw parsing_error("Key missing value");
    if(in[p]=='=') {
        p = in.find_first_not_of(" \r\n\t",p+1);
        if(p==std::string::npos)
            throw parsing_error("Key missing value");
    }
    if(in[p]=='#')
            throw parsing_error("Key missing value");
    
    size_t value_start,value_end;
    char c;
    if((c=in[p])=='"' || c=='\'') {
        p++;
        value_start=p;
        p=in.find_first_of(c,p);
        if(p==std::string::npos)
            throw parsing_error("Incomplete string");
        value_end = p;
        p++;
    }
    else {
        value_start=p;
        p=in.find_first_of(" \t\r\n#",p);
        if(p==std::string::npos)
            p=in.size();
        if(value_start==p)
            throw parsing_error("Empty value");
        value_end=p;
    }
    p=in.find_first_not_of(" \t\r\n",p);
    if(p!=std::string::npos && in[p]!='#') 
        throw parsing_error("Unexpected content after a value");
    value = in.substr(value_start,value_end-value_start);
    key = in.substr(key_start,key_end-key_start);
    p=key.find_first_of("[]");
    if(p!=std::string::npos) {
        if( p==0 ||
            key[p]!='[' 
            || key.back()!=']' 
            || (index = key.substr(p+1,key.size() - p - 2)).find_first_of("[]")!=std::string::npos
            || index.empty())
        {
            throw parsing_error("Key with index format should be `key[index]' ");
        }
        key = key.substr(0,p) + "[]";
    }
    return true;
}

int get_type_index(std::string const &index)
{
    if(index=="sea")
        return water_generator::sea_mark;
    if(index=="land")
        return water_generator::land_mark;
    if(index=="lake")
        return water_generator::lake_mark;
    if(index=="river")
        return water_generator::river_mark;
    return -1;
}


bool yesno(std::string const &key,std::string const &value)
{
    if(value == "yes")
        return true;
    if(value == "no")
        return false;
    throw parsing_error("Invalid key value "  + key + " expected to be yes or no");
}


class logger : public std::streambuf {
public:
    logger(std::ostream &out,std::ofstream &log)
    {
        out_ = &out;
        buf_ = out.rdbuf(this);
        log_ = &log;
    }
    void close()
    {
        if(buf_ != 0) {
            out_->rdbuf(buf_);
            buf_ = 0;
        }
    }
    ~logger()
    {
        close();
    }
    virtual std::streamsize xputn(char const *s,std::streamsize n)
    {
        log_->write(s,n);
        log_->flush();
        std::streamsize r = buf_->sputn(s,n);
        if(buf_->pubsync()!=0)
            return -1;
        return r;
    }
    virtual int overflow(int c)
    {
        if(c!=EOF) {
            char ch=c;
            log_->write(&ch,1);
            log_->flush();
        }
        if(c!=EOF)
            buf_->sputc(c);
        if(buf_->pubsync()!=0)
            return -1;
        if(c!=EOF)
            return c;
        return 0;
    }
private:
    std::streambuf *buf_;
    std::ostream *out_;
    std::ofstream *log_;
};


void load_profile(std::string file_name,std::ofstream &log)
{
    std::ifstream in(file_name.c_str());
    if(!in) {
        throw std::runtime_error("Failed to open file " + file_name);
    }

	char bom[4]={};
	in.read(bom,3);
	if(bom!=std::string("\xEF\xBB\xBF")) {
		in.seekg(0);
	}

    std::string dem_prefix;
    bool via_scale = false;
    bool via_coord = false;
    double scale = -1;
    double lat_c=-1000,lon_c=-1000;
    int line = 0;
    log<< "- Loading Profile from " << file_name << std::endl;
    try {
        while(!in.eof()) {
            std::string s;
            std::getline(in,s);
            line++;
            std::string key,value,index;
            try {
                if(!parse(s,key,value,index))
                    continue;
            }
            catch(...) {
                log << "Exception in line " << line << " :" << s << std::endl;
                throw;
            }
            if(index.empty())
                log << "  " << key << "=" << value << std::endl;
            else
                log << "  " << key << " with " << index << " =" << value << std::endl;

            if(key == "map_size") {
                map_size = atoi(value.c_str());
                switch(map_size) {
                case 64:
                case 128:
                case 256:
                case 512:
                    break;
                default:
                    throw parsing_error("Invalid map size " + value);
                }
            }
            else if(key == "cbm_size") {
                cbm_size = atoi(value.c_str());
            }
			else if(key == "cbm_labels") {
				cbm_labels = value;
			}
            else if(key == "cbm_grid") {
                if(value == "game" || value=="geo" || value=="geo_dd" || value=="geo_dm") {
                    if(value == "geo")
                        cbm_grid_type = "geo_dd";
                    else
                        cbm_grid_type = value;
                }
                else
                    throw parsing_error("Invalid value " + value + "  for cbm_grid, expected one of : game, geo, geo_dd, geo_dm");
            }
			else if(key == "splattype_filter_radius") {
				splattype_filter_radius = atoi(value.c_str());
			}
            else if(key=="type_mapping") {
                custom_mapping = value;
            }
            else if(key=="tile_set") {
				auto valid = util::dir(tile_set_path,".csv");
				if(valid.find(value + ".csv") == valid.end()) {
					std::ostringstream ss;
					ss << "Invalid tile set: " << value << " expected one of:";
					for(std::string const &name : valid)
						ss << name.substr(0,name.size()-4) << " ";
					throw parsing_error(ss.str());
				}
				tile_set = value;
            }
            else if(key=="dem") {
                if(value == "srtm3") {
                    db_type = dem::srtm3();
                }
                else if(value == "srtm30") {
                    db_type = dem::srtm30();
                }
                else if(value=="gtopo30") {
                    db_type = dem::gtopo30();
                }
                else {
                    throw parsing_error("Invalid database type name " + value + " valid are srtm3, srtm30 or gtopo30");
                }
                db_type.altitude_offset = altitude_offset;
            }
            else if(key == "altitude_offset") {
                altitude_offset = db_type.altitude_offset = atoi(value.c_str());
            }
            else if(key == "scale") {
                via_scale = true;
                scale = atof(value.c_str());
            }
            else if(key == "map_name") {
                map_name = value;
            }
			else if(key == "check_updates") {
				check_updates = atoi(value.c_str());
			}
			else if(key == "disable_sse") {
				disable_sse = yesno(key,value);
			}
            else if(key == "shores") {
                shores = value;
            }
            else if(key == "rivers") {
                rivers = value;
            }
            else if(key == "download_sources") {
                download_sources = value;
            }
            else if(key == "download") {
                if(value == "yes")
                    auto_download_enabled = true;
                else if(value == "nossl") {
                    auto_download_enabled = true;
                    disable_ssl_check=true;
                }
                else if(value == "no")
                    auto_download_enabled = false;
                else
                    throw parsing_error("Key download need values yes or no");
            }
            else if(key == "temporary_dir") {
                temp_dir = value;
            }
            else if(key == "grid_dir") {
                grid_dir = value;
            }
			else if(key == "water_alt_limit") {
				water_alt_limit = atoi(value.c_str());
			}
			else if(key == "lake_or_island_min_size") {
				lake_or_island_min_size = atof(value.c_str());
			}
			else if(key == "water_to_land_slope") {
				water_to_land_slope = atof(value.c_str());
			}
   			else if(key == "water_to_land_range") {
				water_to_land_range = atoi(value.c_str());
			}
         else if(key == "fix_river_slopes") {
                if(value == "yes")
                    fix_river_slopes = true;
                else if(value == "no")
                    fix_river_slopes = false;
                else
                    throw parsing_error("Invalid value of fix_river_slopes option, should be yes or no");
            }
            else if(key == "river_width") {
                rv_prop.default_width = atoi(value.c_str());
                rv_prop.level_to_width.clear();
            }
            else if(key == "river_width[]") {
                int level = atoi(index.c_str());
                int width = atoi(value.c_str());
                rv_prop.level_to_width[level]=width;
            }
            else if(key == "river_level") {
                rv_prop.max_level = atoi(value.c_str());
            }
            else if(key == "river_north_shift") {
                river_north_shift = atof(value.c_str());
            }
            else if(key == "river_east_shift") {
                river_east_shift = atof(value.c_str());
            }
            else if(key == "globcover_tiff_path") {
                tiff_file = value;
            }
            else if(key == "dem_path") {
                dem_prefix = value;
            }
            else if(key == "output_dir") {
                output_dir = value;
            }
            else if(key == "elevation_marking") {
                if(value == "shade")
                    altitude_handling = shade_hills;
                else if(value == "dark")
                    altitude_handling = dark_hills;
                else if(value == "disable")
                    altitude_handling = untouched_hills;
                else
                    throw parsing_error("elevation_marking should be one of shade, dark, disable");
            }
            else if(key == "elevation_shading_factor") {
                shading_factor = atof(value.c_str());
            }
            else if(key == "solver_allow_cpu") {
                solver_options.allow_cpu=yesno(key,value);
            }
            else if(key == "solver_force_cpu") {
                solver_options.force_cpu=yesno(key,value);
            }
            else if(key == "solver_threshold") {
                solver_options.threshold = atof(value.c_str());
            }
            else if(key == "solver_platform_id") {
                solver_options.platform_id = atoi(value.c_str());
            }
            else if(key == "solver_initial_grid") {
                // hidden
                solver_options.initial_grid = atoi(value.c_str());
            }
            else if(key == "map_color[]") {
                color c(value);
                if(index=="beach")
                    beach_color=c;
                else if(index=="water")
                    water_color=c;
                else if(index=="grid")
                    grid_color=c;
                else if(index=="farm")
                    map_colors[0]=c;
                else if(index=="forest")
                    map_colors[1]=c;
                else if(index=="snow")
                    map_colors[2]=c;
                else if(index=="grass")
                    map_colors[3]=c;
                else if(index=="desert")
                    map_colors[4]=c;
        else
                    throw parsing_error("Invalid map_color[] index `" + index + "'");
            }
            else if(key == "lat1" || key == "lat2" || key=="lat") {
                double v = atof(value.c_str());
                if(!(-60 < v && v < 90))
                    throw parsing_error("Invalid latitude " + value + " should be in rage [-60,90]");
                if(key == "lat1") {
                    lat1 = v;
                    via_coord = true;
                }
                else if(key == "lat2") {
                    lat2 = v;
                    via_coord = true;
                }
                else {
                    lat_c = v;
                    via_scale = true;
               }
            }
            else if(key == "lon1" || key == "lon2" || key=="lon") {
                double v = atof(value.c_str());
                if(!(-180 < v && v < 180))
                    throw parsing_error("Invalid longitude " + value + " should be in rage [-180,180]");
                if(key == "lon1") {
                    lon1 = v;
                    via_coord = true;
                }
                else if(key == "lon2") {
                    lon2 = v;
                    via_coord = true;
                }
                else if(key == "lon") {
                    lon_c = v;
                    via_scale = true;
                }
            }
            else {
                throw parsing_error("Unknown key `" + key + "' with value `"+value+"'");
            }
        }
    }
    catch(parsing_error const &e) {
        std::ostringstream ss;
        ss << "Error in " << file_name  << ", line " << line << ": " << e.what();
        throw std::runtime_error(ss.str());
    }
    if(via_scale == via_coord) {
        throw std::runtime_error("Configuration error: either lat/lon/scale or lat1/lat2/lon1/lon2 should be specified");
    }
    if(via_coord) {
        if(lat1 == -1000 || lat2 == -1000 || lon1 == -1000 || lon2 == -1000) {
            throw std::runtime_error("Configuration error: the latitude or longitude range and not fully defined");
        }
        lat_c = (lat1+lat2)/2;
        lon_c = (lon1+lon2)/2;
        double dist_miles = (std::max(lat1,lat2)-std::min(lat1,lat2)) * 1.852 * 60 / 1.60934;
        map_scale = dist_miles / map_size;
    }
    else { // via scale
        if(lat_c == -1000 || lon_c == -1000 || scale == -1) {
            throw std::runtime_error("Configuration error: the latitude, longitude or scane are not fully defined");
        }
        map_scale = scale;
        double nmiles = scale * map_size * 1.60934 / 1.852;
        lat1=lat_c - (nmiles / 2 / 60 );
        lat2=lat_c + (nmiles / 2 / 60 );
        double diff = (nmiles / 2 / 60 ) / cos(lat_c / 180 * 3.14159);
        lon1=lon_c - diff;
        lon2=lon_c + diff;
        if(lat1 < -60 || 90 < lat2 )
            throw std::runtime_error("Configuration error: the latitude should be below N90 and above S60");
            
    }
    if(db_type.rows == 0) {
        throw std::runtime_error("Configuration error: undefined dem - should be one of srtm30, srtm3, gtopo30");
    }
    if(db_type.code == "srtm3" && (lat1 > 60 || lat2 > 60)) {
        throw std::runtime_error("Latitude above N60 is not supported with srtm3, please use srtm30");
    }
    if(!dem_prefix.empty()) {
        db_type.directory = dem_prefix;
    }
    double miles_to_lat =  1.60934 / 1.852 / 60;
    double miles_to_lon =  miles_to_lat / cos(lat_c / 180 * 3.14159);
    rv_prop.lat_shift = miles_to_lat * river_north_shift;
    rv_prop.lon_shift = miles_to_lon * river_east_shift;
    
}


void resample_type()
{
    // resample
    int bmp_size = map_size * 8;
    
    types=std::move(image<int>(bmp_size,bmp_size,8,8));
    
    int bh = bare_types.size();
    int bw = bare_types[0].size();
    double r_factor = double(bh) / bmp_size;
    double c_factor = double(bw) / bmp_size;
    for(int r=0;r<bmp_size;r++) {
        for(int c=0;c<bmp_size;c++) {
            int br = int(round(r * r_factor));
            int bc = int(round(c * c_factor));
            if(br >= bh)
                br = bh - 1;
            if(bc >= bw)
                bc = bw - 1;
            types[r][c] = bare_types[br][bc];
        }
    }
}

void write_reference_bmp()
{
    std::string fname = output_dir + "/ground_coverage.bmp";
    outfile f(fname);
    int rows = types.height();
    int cols = types.width();
    bmp::header hdr(rows,cols);

    f.write(&hdr,hdr.offset);
    
	std::vector<bmp::rgbq> c(cols);
    for(int i=types.height()-1;i>=0;i--) {
		for(int j=0;j<cols;j++) {
			unsigned char color = types[i][j];
            f.write(&color,1);
		}
    }
    f.close();
}

void recolor()
{
    int bmp_size = map_size * 8;
    for(int r=0;r<bmp_size;r++) {
        for(int c=0;c<bmp_size;c++) {
			types[r][c]=map_type[types[r][c]];
        }
    }
}

std::vector<int> gconv(std::vector<int> const &in)
{
	std::vector<int> r(in.size() + 2,0);
	for(unsigned i=0;i<in.size();i++) {
		r[i+0] += in[i];
		r[i+1] += in[i] * 2;
		r[i+2] += in[i];
	}
	return r;
}

std::vector<float> prod(std::vector<int> const &in)
{
	int N=in.size();
	std::vector<float> kernel(N*N);
	int sum = 0;
    int pos = 0;
	for(int r=0;r<N;r++) {
		for(int c=0;c<N;c++) {
			int p=in[r]*in[c];
            kernel[pos++] = p;
			sum += p;
		}
	}
    pos = 0;
	for(int r=0;r<N;r++) {
		for(int c=0;c<N;c++) {
			kernel[pos++] /= sum;
		}
	}
	return kernel;
}

std::vector<float> get_kernel(int radius)
{
	std::vector<int> sec(1,1);
	for(int i=0;i<radius;i++)
		sec = gconv(sec);
	return prod(sec);
}

image<int> expanded_types(int radius)
{
	int rows = types.height();
	int cols = types.width();
    image<int> safetypes(rows + radius*2,cols + radius*2,8,8);
    for(int r=0;r<rows;r++) {
        for(int c=0;c<cols;c++) {
            safetypes[r+radius][c+radius]=types[r][c];
        }
    }
    for(int r=0;r<rows;r++) {
        for(int dc=0;dc<radius;dc++) {
            safetypes[r+radius][dc]=types[r][0];
            safetypes[r+radius][radius + cols + dc]=types[r][cols-1];
        }
    }
    for(int dr=0;dr<radius;dr++) {
        for(int c=0;c<cols+radius*2;c++) {
            safetypes[dr][c]=safetypes[radius][c];
            safetypes[radius + rows + dr][c]=safetypes[rows + radius - 1][c];
        }
    }
    return std::move(safetypes);
}

float low_pass_filter_splattype(int radius)
{
	if(radius <= 0)
		return 0;
	auto start_ts = std::chrono::high_resolution_clock::now();


	std::vector<float> kernel = get_kernel(radius);
    std::vector<unsigned> ikernel(kernel.size());
    for(size_t i=0;i<ikernel.size();i++) {
        ikernel[i] = kernel[i] / 16 * UINT_MAX;
    }
    int eradius = (radius + 3)/4*4;

    image<int> safetypes = std::move(expanded_types(eradius));

	int rows = types.height();
	int cols = types.width();
	
	std::function<void(int,int)> runner;
	
    #ifdef USE_SIMD
	
	if(!disable_sse 
		&& __builtin_cpu_supports("sse3")
		&& __builtin_cpu_supports("sse4.1") 
		&& __builtin_cpu_supports("sse4.2"))
	{

		typedef unsigned int int4 __attribute__ ((vector_size(16)));
		typedef unsigned int int4u __attribute__ ((vector_size(16),aligned(4)));
		
		runner = [&](int from,int to) __attribute__((force_align_arg_pointer,target("sse3"),target("sse4.1"),target("sse4.2"))) {	
			for(int r=from;r<to;r++) {
				for(int c=0;c<cols;c+=4) {
					int ref_r = r+eradius;
					int ref_c = c+eradius;
					int4 sr[3][3] = {
						{ *(int4u*)&safetypes[ref_r-1][ref_c-1],*(int4u*)&safetypes[ref_r-1][ref_c+0],*(int4u*)&safetypes[ref_r-1][ref_c+1] },
						{ *(int4u*)&safetypes[ref_r  ][ref_c-1],*(int4u*)&safetypes[ref_r  ][ref_c+0],*(int4u*)&safetypes[ref_r  ][ref_c+1] },
						{ *(int4u*)&safetypes[ref_r+1][ref_c-1],*(int4u*)&safetypes[ref_r+1][ref_c+0],*(int4u*)&safetypes[ref_r+1][ref_c+1] },
					};
					int4 blue =
									(sr[0][0] & 0x1) + 2*(sr[0][1] & 0x1) +   (sr[0][2] & 0x1)
								+ 2*(sr[1][0] & 0x1) + 4*(sr[1][1] & 0x1) + 2*(sr[1][2] & 0x1)
								+   (sr[2][0] & 0x1) + 2*(sr[2][1] & 0x1) +   (sr[2][2] & 0x1);
					
					int4 center = sr[1][1];
					
					blue = (blue * 255 + 8)/ 16;

					if(blue[0] == 255 && blue[1] == 255 && blue[2]==255 && blue[3] == 255) {
						continue;
					}

					int4 red_sum    = {0,0,0,0};
					int4 red_weight = {0,0,0,0};
					int kp = 0;
					for(int dr=-radius,kr=0;dr<=radius;dr++,kr++) {
						for(int dc=-radius,kc=0;dc<=radius;dc++,kc++) {
							int4 type = *(int4u*)&safetypes[ref_r+dr][ref_c+dc];
							unsigned int kv = ikernel[kp++];
							int4 weight = (type & 1) == 0 ? int4{kv,kv,kv,kv} : int4{0,0,0,0};
							red_weight += weight;
							red_sum    += weight * (type>>(16+4));
						}
					}
					int4 safe_red_weight = red_weight == 0 ? int4{1,1,1,1} : red_weight;
					int4 red = red_weight == 0 ? int4{0,0,0,0} : ((red_sum + red_weight/2) / safe_red_weight) * 16;
					int4 green = (center >> 8) & 0xFF;
					
					*(int4u*)&types[r][c] = blue == 255 ? center : ((red << 16) + (green << 8) + blue);
				}
			}
		}; // runner*/
		
		std::cout << " using sse4.2 " << std::flush;
    }
	else 
    #endif
	{
		
		runner = [&](int from,int to) {	
			for(int r=from;r<to;r++) {
				for(int c=0;c<cols;c++) {
					int ref_r = r+eradius;
					int ref_c = c+eradius;
					int *sr[3] = {
						&safetypes[ref_r-1][ref_c-1],
						&safetypes[ref_r][ref_c-1],
						&safetypes[ref_r+1][ref_c-1]
					};
					unsigned blue =
									(sr[0][0] & 0x1) + 2*(sr[0][1] & 0x1) +   (sr[0][2] & 0x1)
								+ 2*(sr[1][0] & 0x1) + 4*(sr[1][1] & 0x1) + 2*(sr[1][2] & 0x1)
								+   (sr[2][0] & 0x1) + 2*(sr[2][1] & 0x1) +   (sr[2][2] & 0x1);
					blue = (blue * 255 + 8)/ 16;
					if(blue == 255)
						continue;

					unsigned red_sum    = 0;
					unsigned red_weight = 0;
					int kp = 0;
					for(int dr=-radius,kr=0;dr<=radius;dr++,kr++) {
						for(int dc=-radius,kc=0;dc<=radius;dc++,kc++) {
							unsigned type = safetypes[ref_r+dr][ref_c+dc];
							unsigned sum_factor  = (type & 1)  == 0; // blue == 0 for blue one of 0xFF or 0
							unsigned weight = sum_factor * ikernel[kp++];
							red_weight += weight;
							red_sum    += weight * (type >> (16 + 4));
						}
					}
					
					unsigned red = ((red_sum + red_weight/2)/red_weight)*16;
					unsigned green = (safetypes[ref_r][ref_c] >> 8) & 0xFF;
					types[r][c] = (red << 16) + (green << 8) + blue;
				}
			}
		};
	}
   

    // split for multiple threads
	int threads = std::thread::hardware_concurrency();
	std::vector<std::thread> tasks(threads);
	int chunk = (rows + threads-1)/threads;
    for(int start=0,id=0;start<rows;start+=chunk,id++) {
        int limit = std::min(start + chunk,rows);
        tasks[id]=std::move(std::thread(runner,start,limit));
    }

    for(int i=0;i<threads;i++) {
        tasks[i].join();
    }

    auto end_ts = std::chrono::high_resolution_clock::now();
	float time = std::chrono::duration_cast<std::chrono::duration<float, std::ratio<1> > >(end_ts-start_ts).count();
    return time;
}


void write_gndtype()
{
    std::string fname = output_dir + "/splattype.bmp";
    outfile f(fname);
    int rows = types.height();
    int cols = types.width();
    bmp::header hdr(4096,4096,32);
    f.write(&hdr,hdr.offset);
    std::vector<bmp::rgbq> c(cols);
    int padding = (4096 - rows)/2;
    
    std::vector<unsigned char> zeros(4096*4,0);
    for(int i=0;i<padding;i++)
        f.write(&zeros[0],zeros.size());
    for(int i=rows-1;i>=0;i--) {
        f.write(&zeros[0],padding*4);
		for(int j=0;j<cols;j++) {
			int color = types[i][j];
			c[j].r = (color >> 16)& 0xFF;
			c[j].g = (color >> 8)& 0xFF;
			c[j].b = (color >> 0)& 0xFF;
			c[j].res = 0;
		}
        f.write(&c[0],cols*4);
        f.write(&zeros[0],padding*4);
    }
    for(int i=0;i<padding;i++)
        f.write(&zeros[0],zeros.size());
    f.close();
}

unsigned altitude_to_bmp(std::string file,std::vector<std::vector<int16_t> > const &elev)
{
    int max = 0;
    int matrix_size = map_size * 8;
    for(int i=0;i<matrix_size;i++) {
        for(int j=0;j<matrix_size;j++) {
            if(elev[i][j] > max)
                max = elev[i][j];
        }
    }
    outfile f(file);
    bmp::header hdr(map_size * 8,map_size * 8);
    f.write(&hdr,sizeof(hdr));
    std::vector<uint8_t> row(matrix_size,0);
    unsigned div = max;
    if(div == 0)
        div = 1;
    for(int i=matrix_size-1;i>=0;i--) {
        for(int j=0;j<matrix_size;j++) {
            int v=elev[i][j];
            if(v<0)
                v=0;
            row[j] = v * 255u / div;
        }
        f.write(&row[0],matrix_size);
    }
    f.close();
    return max;
}

void save_elevations_file()
{
    std::string elev_file = output_dir +"/" + map_name + ".elv";
    outfile f(elev_file);
    int matrix_size = map_size * 2;
    int padding = (1024 - matrix_size) / 2;
    std::vector<int16_t> zero_row(1024,0);
    for(int i=0;i<padding;i++) {
        f.write(&zero_row[0],2*1024);
    }
    for(int i=matrix_size-1;i>=0;i--) {
        f.write(&zero_row[0],2*padding);
        f.write(&elevations[i][0],2*matrix_size);
        f.write(&zero_row[0],2*padding);
    }
    for(int i=0;i<padding;i++) {
        f.write(&zero_row[0],2*1024);
    }
    f.close();
}

void save_height_file()
{
    std::string elev_file = output_dir +"/" + map_name + ".raw";
    outfile f(elev_file);
	int size = elevations.size();
    static const int raw_size = 4096;
	int padding = (raw_size - size)/2;
	std::vector<int16_t> zero(raw_size, int16_t(-water_to_land_range));
    char const *zero_ptr = reinterpret_cast<char*>(&zero[0]);
	if(padding > 0) {
		for(int i=0;i<padding;i++)
			f.write(zero_ptr,raw_size*2);
	}
    for(int i=elevations.size()-1;i>=0;i--) {
		if(padding > 0)
			f.write(zero_ptr,padding*2);
        f.write(&elevations[i][0],elevations[i].size()*2);
		if(padding > 0)
			f.write(zero_ptr,padding*2);
    }
	if(padding > 0) {
		for(int i=0;i<padding;i++)
			f.write(zero_ptr,raw_size*2);
	}
    f.close();
}


int row_from_lat(double lat,int height)
{
    return int((90 - lat)/(90 - lower_lat) * height);
}

int col_from_lon(double lon,int width)
{
    return int((lon + 180) / 360 * width);
}

void load_globcover_data()
{
    TIFFSetWarningHandler(0);
    downloader::manager::instance().check(tiff_file,"globcover");
    TIFF *in = TIFFOpen(tiff_file.c_str(),"r");
    if(!in) 
        throw std::runtime_error("Failed to open file " + tiff_file);
    uint32 imw,imh;
    TIFFGetField(in, TIFFTAG_IMAGELENGTH, &imh);
    TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &imw);
    size_t len = TIFFScanlineSize(in);

    if(len < imw) {
        TIFFClose(in);
        std::ostringstream ss;
        ss << "Internal error line length = " << len << " < im-width" << imw;
        throw std::runtime_error(ss.str());
    }
    int width = imw;
    int height = imh;

    int r0=row_from_lat(lat1,height);
    int r1=row_from_lat(lat2,height);
    if(r1 < r0)
        std::swap(r1,r0);
    int c0=col_from_lon(lon1,width);
    int c1=col_from_lon(lon2,width);
    if(c1 < c0)
        std::swap(c1,c0);

    std::vector<unsigned char> buf(len);        

    int h = r1-r0 + 1;
    int w = c1-c0 + 1;

    bare_types.resize(h,std::vector<unsigned char>(w,0));

    for(int r=r0,inr=0;r<=r1;r++,inr++) {
        TIFFReadScanline(in,&buf[0],r);
        if(0<= c0 && c1 <= width)
            memcpy(&bare_types[inr][0],&buf[c0],w);
        else if(c0 < 0) {
            memcpy(&bare_types[inr][0],&buf[width + c0],-c0);
            memcpy(&bare_types[inr][-c0],&buf[0],w + c0);
        }
        else { // if c1 > width
            memcpy(&bare_types[inr][0],&buf[c0],(width - c0));
            memcpy(&bare_types[inr][width - c0],&buf[0],c1 - width + 1);
        }
    }
    TIFFClose(in);

}

void land_border()
{
    int border_size = 2;
    for(int d=0;d<border_size;d++) {
        int r1=d;
        int r2=(map_size*2 - 1) - d;
        for(int c=0;c<map_size * 2;c++) {
            elevations[r1][c]=0;
            elevations[r2][c]=0;
            elevations[c][r1]=0;
            elevations[c][r2]=0;
        }
    }
}

int get_min_elev(int vr,int vc)
{
    return  std::min(
                std::min(elevations[vr][vc],elevations[vr][vc+1]),
                std::min(elevations[vr+1][vc],elevations[vr+1][vc+1])
            );
}

int get_max_elev(int vr,int vc)
{
    return  std::max(
                std::max(elevations[vr][vc],elevations[vr][vc+1]),
                std::max(elevations[vr+1][vc],elevations[vr+1][vc+1])
            );
}

int water_to_elev_row(int wr)
{
    int vr = (wr - 14) / 16;
    return vr;
}

int water_to_elev_col(int wc)
{
    int vc = (wc - 2) / 16;
    return vc;
}


int elev_to_water_row(int vr)
{
    int wr = vr * 16 + 14;
    return wr;
}

int elev_to_water_col(int vc)
{
    int wc = vc * 16 + 2;
    return wc;
}

struct ground_surface {
    int r,c;
    bool has_lake;
    bool was_adjusted;
    ground_surface(int rin=0,int cin=0,bool lake=false) : 
        r(rin),
        c(cin),
        has_lake(lake),
        was_adjusted(false)
    {
    }
};

typedef std::vector<ground_surface> surface_set;

surface_set get_water_set(water_generator &gen)
{
    surface_set water_elevations;
    int max_water = map_size * 32;
    int elev_size = map_size * 2;
    for(int vr = 0; vr < elev_size - 1;vr++) {
        for(int vc=0;vc < elev_size - 1; vc++) {
            if( get_max_elev(vr,vc) == 0)
                continue;
            bool has_water = false;
            bool has_lake = false;
            for(int dr=0;dr<16 && !has_lake;dr++) {
                for(int dc=0;dc<16 && !has_lake;dc++) {
                    int wr = elev_to_water_row(vr) + dr;
                    int wc = elev_to_water_col(vc) + dc;
                    if(wr < 0 || wr >= max_water || dc <0 || dc >= max_water)
                        continue;
                    int type = gen.pixel(wr,wc);
                    if(type!= water_generator::land_mark)
                        has_water = true;
                    if(type== water_generator::lake_mark)
                        has_lake = true;
                }
            }
            if(has_water)
                water_elevations.push_back(ground_surface(vr,vc,has_lake));
        }
    }
    return water_elevations;
}
/*
bool run_correction_iteration(surface_set &surfaces)
{
    // The maximal allowed slope is 120 feet per 5280 feet (i.e. 1 mile). Elevations grid is
    // 0.5 mile so maximal allowed slope is 60 feet per grid. However the higher slope can be
    // achieved if vertices placed in triangle as LHL such that the shortest distance is at central point and it is
    // sqrt(0.5^2 + 0.5^2) ~= 707. So we need to have maximal allowed drop of 60 * 0.707 feed per 4 vertices which is
    // 42 feet
    static const int max_drop = 42; 
    bool ok;
    ok = true;
    for(size_t i=0;i<surfaces.size();i++) {
        ground_surface &s = surfaces[i];
        int min = get_min_elev(s.r,s.c);
        int max = get_max_elev(s.r,s.c);
        int local_max_drop = max_drop - 1;
        if(s.has_lake)
            local_max_drop = 0;
        if(max - min > local_max_drop) {
            int16_t limit = min + local_max_drop;
            for(int dr=0;dr<2;dr++)
                for(int dc=0;dc<2;dc++)
                    elevations[s.r+dr][s.c+dc] = std::min(elevations[s.r+dr][s.c+dc],limit);
            s.was_adjusted = true;
            ok = false;
        }
    }
    return ok;
}
*/
/*
void fix_river_elevations(water_generator &gen)
{
    std::cout << "- Fixing river and lake slopes... " << std::flush;

    // all surfaces containing water above altitude 0
    surface_set surfaces = get_water_set(gen);
    
    std::vector<std::vector<int16_t> > saved_elevations(elevations);
    // run corrections
    int iterations = 0;
    while(!run_correction_iteration(surfaces))
    {
        iterations ++;
    }
    int adjusted = 0;
    for(size_t i=0;i<surfaces.size();i++)
        if(surfaces[i].was_adjusted)
            adjusted++;
    //save file
    int vertices = 0;
    for(int r=0;r<map_size*2;r++) {
        for(int c=0;c<map_size*2;c++) {
            assert(saved_elevations[r][c] >= elevations[r][c]);
            int16_t diff = saved_elevations[r][c] - elevations[r][c];
            if(diff!=0) {
                vertices++;
            }
            saved_elevations[r][c] = diff;
            
        }
    }
    unsigned max = altitude_to_bmp(output_dir+"/water_altitude_correction.bmp", saved_elevations);
    std::cout   << " completed" << std::endl 
                << "   Total iteration:         " << std::setw(6) << iterations  << std::endl
                << "   Maximal alt. correction :" << std::setw(6) << max << " feet"  << std::endl
                << "   Total surfaces adjusted :" << std::setw(6) << adjusted  << " of " << surfaces.size() << " containing non-0 alt water"<< std::endl
                << "   Total vertices adjusted :" << std::setw(6) << vertices  << std::endl;
}

void fix_sea_elevations(water_generator &gen)
{
    int max_water = map_size * 32;
    int elev_size = map_size * 2;
    for(int vr = 0; vr < elev_size - 1;vr++) {
        for(int vc=0;vc < elev_size - 1; vc++) {
            for(int dr=0;dr<16;dr++) {
                for(int dc=0;dc<16;dc++) {
                    int wr = elev_to_water_row(vr) + dr;
                    int wc = elev_to_water_col(vc) + dc;
                    if(wr < 0 || wr >= max_water || dc <0 || dc >= max_water)
                        continue;
                    if(gen.pixel(wr,wc)==water_generator::sea_mark) {
                        elevations[vr][vc] = 0;
                        elevations[vr][vc+1] = 0;
                        elevations[vr+1][vc] = 0;
                        elevations[vr+1][vc+1] = 0;
                        goto internal_loop_exit;
                    }
                }
            }
            internal_loop_exit:
                ;
        }
    }
}
*/
void make_map_color_index(bmp::header &hdr)
{
    
    for(int c=0;c<5;c++) {
        for(int i=0;i<50;i++) {
            double f=i/50.0 + 0.5;
            bmp::rgbq *p = &hdr.ih.colors[c*50 + i ];
            p->r = std::min(std::max(0,int(floor(map_colors[c].r * f))),255);
            p->g = std::min(std::max(0,int(floor(map_colors[c].g * f))),255);
            p->b = std::min(std::max(0,int(floor(map_colors[c].b * f))),255);
        }
    }
    // water
    {
        bmp::rgbq *p = &hdr.ih.colors[255];
        p->r=water_color.r;
        p->g=water_color.g;
        p->b=water_color.b;
    }
    // beach
    {
        bmp::rgbq *p = &hdr.ih.colors[254];
        p->r=beach_color.r;
        p->g=beach_color.g;
        p->b=beach_color.b;
    }
    // grid
    {
        bmp::rgbq *p = &hdr.ih.colors[253];
        p->r=grid_color.r;
        p->g=grid_color.g;
        p->b=grid_color.b;
    }
    // text
    {
        bmp::rgbq *p = &hdr.ih.colors[252];
        p->r=0;
        p->g=0;
        p->b=0;
    }
}

int get_color_from_type(int type,double brightness_factor)
{
	
	if((type & 0xFF) >= 127) {
		type = 16+((type >> 8) / 0x40 & 0x3);
	}
	else {
		type = type >> (8+8+4);
	}
	assert(0<=type && type<=19);

    // 255 - water
    // 254 - beach
    // 253 - grid line
    // 252 - TBD
    // 251 - TBD
    // 250 - TBD
    // 200-249 - beach, river bad, rock
    // 150-199 - grass, grocky grass, sandy grass, swamp
    // 100-149 - snow
    //  50- 99 - forest
    //   0- 49 - farm 
    
    int basic = 4 - (type / 4);
/*    switch(type & 0xF) {
    case 4:
    case 5:
        basic = 0; // farm1,2
        break;
    case 2:
    case 3:
        basic = 1; // forest1,2 
        break;
    case 12:
        basic = 2; // snow
        break;
    case 1:
    case 7:
    case 8:
    case 9:
        basic = 3; // grass, grocky grass, sandy grass, swamp
        break;
    case 6:
    case 10:
    case 11:
        basic = 4; // beach, river bad, rock
        break;
    case 0:
        return 255;
    default:
        throw std::runtime_error("Internal error invalid ground type");
    }*/
    int brightness=std::max(0,std::min(49,int(round(brightness_factor*49))));
    int color = 50 * basic + brightness;
    return color;
}

double get_scaled_elevation(double r,double c)
{
    int size = map_size * 8;
    r*=size;
    c*=size;
    int rl = std::max(0,std::min(int(floor(r)),size-1));
    int cl = std::max(0,std::min(int(floor(c)),size-1));
    int rh = std::min(rl+1,size-1);
    int ch = std::min(cl+1,size-1);
    if(r < rl)
        r=rl;
    if(c < cl)
        c=cl;
    double rwh = r - rl;
    double cwh = c - cl;
    double rwl = 1-rwh;
    double cwl = 1-cwh;
    double elev =   rwl*(cwl*elevations[rl][cl] + cwh * elevations[rl][ch]) 
                  + rwh*(cwl*elevations[rh][cl] + cwh * elevations[rh][ch]);
    return elev;
}

std::vector<std::vector<bool> > load_grid(size_t tsize)
{
    std::ostringstream ss;
    ss << grid_dir << "/grid_" << tsize <<"x"<<tsize<<"px_" << map_size<<"miles.tif";
    std::string file_name = ss.str();
    std::vector<std::vector<bool> > grid(tsize,std::vector<bool>(tsize,false));
    TIFF *in = 0;
    try {
        std::vector<unsigned char> buf(tsize);
        in = TIFFOpen(file_name.c_str(),"r");
        if(!in) {
            throw std::runtime_error("Failed to open file " + file_name);
        }
        uint32 imw,imh;
        TIFFGetField(in, TIFFTAG_IMAGELENGTH, &imh);
        TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &imw);
        size_t len = TIFFScanlineSize(in);
        if(imw != tsize || imh!=tsize || len!=tsize)
            throw std::runtime_error("Invalid tiff format for " + file_name);
        for(unsigned i=0;i<tsize;i++) {
            TIFFReadScanline(in,&buf[0],i);
            for(unsigned j=0;j<tsize;j++) {
                bool is_grid = buf[j]<128;
                grid[i][j]=is_grid;
            }
        }
        TIFFClose(in);
        in=0;
    }
    catch(...) {
        if(in)
                TIFFClose(in);
        throw;
    }
    return grid;
}

std::string to_deg(double v,char Pos,char Neg,int acc)
{
    std::ostringstream ss;
    int total_seconds = int(round((v)*3600));
	bool neg = total_seconds < 0;
	total_seconds = labs(total_seconds);
	
    int degree  = total_seconds / 3600;
    int minutes = (total_seconds % 3600) / 60;
    int seconds = total_seconds % 60;
    ss << abs(degree) << "\xB0";
    if(acc >= 1 && minutes!=0)
        ss << minutes << "'";
    if(acc >= 2 && seconds!=0)
        ss << seconds << '"';
    if(neg)
        ss << Neg;
    else
        ss << Pos;
    return ss.str();
}

std::string get_coord(double lat,double lon,int acc_lan,int acc_lon)
{
    std::ostringstream ss;
    if(lon > 180)
        lon -= 360;
    else if(lon < -180)
        lon += 360;
    if(cbm_grid_type == "geo_dd")
        ss <<std::fixed  << std::setprecision(acc_lan) << lat << ", " << std::setprecision(acc_lon) << lon;
    else
        ss << to_deg(lat,'N','S',acc_lan) << ", " << to_deg(lon,'E','W',acc_lon);
    return ss.str();
}


int get_coord_spacing(int acc_lat,int acc_lon)
{
    double lats[4]={lat1,lat1,lat2,lat2};
    double lons[4]={lon1,lon2,lon1,lon2};
    int pixels = 1;
    for(int i=0;i<4;i++) {
        std::string coord = get_coord(lats[i],lons[i],acc_lat,acc_lon);
        pixels =std::max(get_print_str_len(coord.c_str()),pixels);
    }
    return pixels;
}


struct tick_grid {
    int digits;
    double tick;
    double tick2;
    double tick3;
    double start;
    double end;
};

tick_grid get_optimal_grid_step(double v1,double v2,int image_size,int pixels_needed = 0)
{
    bool decimal_grid = cbm_grid_type == "geo_dd";

    double vmin = std::min(v1,v2);
    double vmax = std::max(v1,v2);
    double total_steps = double(image_size) / pixels_needed;
    double range = (vmax - vmin) / total_steps;
    int    power = ceil(log(range)/log(10));
    double tick = pow(10.0,power);
    double res = range / tick;
    tick_grid r;
    r.digits = - power;
    if(res < 0.2) {
        tick /= 5;
        r.digits ++;
    }
    else if(res < 0.5) {
        tick /= 2;
        r.digits ++;
    }
    if(r.digits < 0)
        r.digits = 0;

    if(decimal_grid || tick >= 5) {
        r.tick = tick;
        r.tick2 = tick / 5;
        r.tick3 = tick / 10;
    }
    else {
        if(tick < 1) {
            int factors[]={1,2,3,5,10,20,30,60, 2*60,3*60,5*60,10*60,20*60,30*60,60*60};
            int total_factors = sizeof(factors)/sizeof(int);
            int pos = 0;
            while(pos + 1 < total_factors && range / (1.0 / factors[pos+1]) < 1)
                pos++;
            tick = 1.0 / factors[pos];
			if(factors[pos] <= 1)
				r.digits = 0;
			else if(factors[pos] <= 60)
				r.digits = 1;
			else
				r.digits = 2;
        }
        r.tick = tick;
        r.tick2 = tick / 6;
        r.tick3 = tick / 12;
    }
    r.start = floor( vmin/tick) * tick;
    r.end   = ceil(  vmax/tick) * tick;
    return r;
}

std::vector<std::vector<bool> > load_grid_lat_lon(int tsize)
{
    std::vector<std::vector<bool> > grid(tsize,std::vector<bool>(tsize));
    static const int y_spacing = 12 * 10;
    tick_grid gy = get_optimal_grid_step(lat1,lat2,tsize,y_spacing);
    tick_grid gx;
    int digits_lon = 0;
    for(;;) {
        int spacing;
        for(;;){
            spacing = get_coord_spacing(gy.digits,digits_lon) + 15;
            if(spacing * 2 > y_spacing)
                break;
            digits_lon++;
        }
        tick_grid gx_tmp = get_optimal_grid_step(lon1,lon2,tsize,spacing);
        if(gx_tmp.digits <= digits_lon) {
            gx = gx_tmp;
            break;
        }
        digits_lon ++;
    }
    double row_factor = tsize / (lat2 - lat1);
    double col_factor = tsize / (lon2 - lon1);
    auto lat2row=[=](double lat) -> int{
        return round((lat-lat1) * row_factor);
    };
    auto lon2col=[=](double lon) -> int{
        return round((lon-lon1) * col_factor);
    };

    double ytick2 = gy.tick2;
    double xtick2 = gx.tick2;
    
    double ytick3 = gy.tick3;
    double xtick3 = gx.tick3;

    for(double lat=gy.start;lat<=gy.end;lat+=gy.tick) {
        int row = lat2row(lat);
        if(0<=row && row < tsize) {
            for(int c=0;c<tsize;c++)
                grid[row][c]=true;
        }
    }
    for(double lon=gx.start;lon<=gx.end;lon+=gx.tick) {
        int col = lon2col(lon);
        if(0<=col && col < tsize) {
            for(int r=0;r<tsize;r++)
                grid[r][col]=true;
        }
    }

    for(double lat=gy.start;lat<=gy.end;lat+=ytick2) {
        int row = lat2row(lat);
        if(row < 1 || row +1 >= tsize) 
            continue;
        for(double lon=gx.start;lon<=gx.end;lon+=xtick2) {
            int col = lon2col(lon);
            if(col < 1 || col + 1 >= tsize) 
                continue;
            grid[row][col]  =true;
            grid[row][col+1]=true;
            grid[row][col-1]=true;
            grid[row+1][col]=true;
            grid[row-1][col]=true;
        }
    }

    for(double lat=gy.start;lat<=gy.end;lat+=ytick3) {
        int row = lat2row(lat);
        for(double lon=gx.start;lon<=gx.end;lon+=xtick2) {
            int col = lon2col(lon);
            if(col < 1 || col + 1 >= tsize || row < 1 || row +1 >= tsize) 
                continue;
            grid[row][col]  =true;
        }
    }
    
    for(double lat=gy.start;lat<=gy.end;lat+=ytick2) {
        int row = lat2row(lat);
        for(double lon=gx.start;lon<=gx.end;lon+=xtick3) {
            int col = lon2col(lon);
            if(col < 1 || col + 1 >= tsize || row < 1 || row +1 >= tsize) 
                continue;
            grid[row][col]  =true;
        }
    }
    
    for(double lat=gy.start;lat<=gy.end;lat+=gy.tick) {
        int row = lat2row(lat);
        for(double lon=gx.start;lon<=gx.end;lon+=gx.tick) {
            int col = lon2col(lon);
            std::string coord = get_coord(lat,lon,gy.digits,gx.digits);
            print_str(coord.c_str(),row, col + 2,1,grid);
        }
    }
    return grid;
}

void make_clipboard_map(int max_elev,std::vector<std::vector<int16_t> > const &elev)
{
    int tsize;
    if(cbm_size == -1) {
        if(map_size == 512)
            tsize = 2048;
        else
            tsize = 1024;
    }
    else
        tsize = cbm_size;

    std::vector<std::vector<bool> > grid;
    try {
        if(cbm_grid_type == "game") 
            grid = std::move(load_grid(tsize));
        else 
            grid = std::move(load_grid_lat_lon(tsize));
    }
    catch(std::exception const &e) {
        std::cout << "\n    Failed to load grid texture: " << e.what() 
                  << "\n    falling back to 1024x1024 size\n";
        
        tsize=1024;
        std::vector<std::vector<bool> > tmp = load_grid(1024);
        grid.swap(tmp);
    }
	
	std::vector<std::vector<bool> > labels;
	if(!cbm_labels.empty()) {
		labels.resize(tsize,std::vector<bool>(tsize));
		render_text rnd(lat1,lat2,lon1,lon2,tsize);
		rnd.load(cbm_labels);
		rnd.optimize();
		rnd.render(labels);
	}

    bmp::header h(tsize,tsize);
    make_map_color_index(h);
    std::string map_file = output_dir +"/" + map_name + ".bmp";
    outfile f(map_file);
    f.write(&h,sizeof(h));
    int factor_type = map_size * 8;
    int water_step = std::max(factor_type / tsize,1);
    for(int r=tsize-1;r>=0;r--) {
        std::vector<unsigned char> colors(tsize);
        for(int c=0;c<tsize;c++) {
				if(!labels.empty() && labels[r][c]) {
					colors[c] = 252;
					continue;
				}
                if(grid[r][c]) {
                    colors[c]=253; // grid;
                    continue;
                }
                int t_r = factor_type * r / tsize;
                int t_c = factor_type * c / tsize;
                int type = types[t_r][t_c];

                bool has_water=false;
                bool has_ground=false;

                int t_r_b = std::max(t_r,0);
                int t_c_b = std::max(t_c,0);
                int t_r_e = std::min(t_r + water_step + 1,factor_type);
                int t_c_e = std::min(t_c + water_step + 1,factor_type);
                
                for(int wr=t_r_b;wr<t_r_e;wr++) {
                    for(int wc=t_c_b;wc<t_c_e;wc++) {
						int alt = elev[wr][wc];
						if(alt <= 0)
							has_water = true;
						else if(alt > 0)
							has_ground = true;
                    }
                }
                int color;
                if(has_water && has_ground) {
                    color=254; // beach
                }
                else if(has_water) {
                    color=255; // water
                }
                else {
                    
                    double brightness = 0.5;
                    if(altitude_handling == shade_hills) {
                        double e_r_1 = double(r-1) / tsize;
                        double e_c_1 = double(c-1) / tsize;
                        double e_r_2 = double(r+1) / tsize;
                        double e_c_2 = double(c+1) / tsize;
                        double e_1=get_scaled_elevation(e_r_1,e_c_1);
                        double e_2=get_scaled_elevation(e_r_2,e_c_2);
                        double emboss = (e_2-e_1) * 2.0 * shading_factor / map_scale / max_elev;
                        brightness=std::max(0.0,std::min(1.0,(0.5+emboss)));
                    }
                    else if(altitude_handling == dark_hills) {
                        double relative_altitude = get_scaled_elevation(double(r)/tsize,double(c)/tsize) / max_elev;
                        brightness = 0.5 * (1.0 - 0.6 * shading_factor * relative_altitude);
                        if(brightness < 0.0)
                            brightness = 0.0;
                    }
                    
                    brightness=std::max(0.0,std::min(1.0,brightness));
                    color = get_color_from_type(type,brightness);
                }
                colors[c]=color;
        }
        f.write(&colors[0],tsize);
    }
    f.close();
}

int get_elevation_from_table(double lat,double lon,std::vector<std::vector<int16_t> > const &el)
{
    int el_size = map_size * 2;
    int r = int(round((lat - lat2)/(lat1 - lat2) * el_size));
    int c = int(round((lon - lon1)/(lon2 - lon1) * el_size));
    if(r < 0 || r>= el_size || c < 0 || c>= el_size)
        return 0;
    return el[r][c];
}


int get_elevation(double lat,double lon)
{
    return get_elevation_from_table(lat,lon,elevations);
}

struct river_mark_callback {
    std::vector<std::vector<int16_t> > const *data;
    std::map<int,int> *ignore_set;
    int limit;
    
    void operator()(int id,int segment,int r,int c) const
    {
        int el_size = map_size * 2;
        r = water_to_elev_row(r);
        c = water_to_elev_col(c);
        for(int vr=r;vr<=r+1;vr++) {
            for(int vc=c;vc<=c+1;vc++) {
                if(vr < 0 || vr>= el_size || vc < 0 || vc>= el_size)
                    continue;
                int diff = (*data)[vr][vc] - elevations[vr][vc];
                if(diff > limit) {
                    if(ignore_set->find(id) == ignore_set->end()) {
                        if(remove_entire_river)
                            ignore_set->insert(std::make_pair(id,std::numeric_limits<int>::max()));
                        else
                            ignore_set->insert(std::make_pair(id,segment+1));
                    }
                }
            }
        }
    }
};

struct removed_river_callback {
    std::vector<std::vector<int16_t> > *data;
    std::map<int,int> *ignore_set;
    
    void operator()(int id,int segment,int r,int c) const
    {
        auto p = ignore_set->find(id);
        if(p == ignore_set->end())
            return;
        if(segment>=p->second)
            return;
        int el_size = map_size * 2;
        r = water_to_elev_row(r);
        c = water_to_elev_col(c);
        if(r < 0 || r>= el_size || c < 0 || c>= el_size)
            return;
        (*data)[r][c]=1;
    }
};
/*
void pass_one()
{
    std::vector<std::vector<int16_t> > fully_saved(elevations);
    std::map<int,int> ignore_set;
    std::cout << "- Pass 1: Collecting required water slope corrections " << std::endl;
    water_generator gen(lat1,lat2,lon1,lon2,map_size * 32);

    std::cout << "  --  Loading & processing shores & river data... " << std::endl;
    gen.load_land(shores);
    gen.load_rivers(rivers,rv_prop);
    gen.make_border();


    std::cout << "  --  Fixing sea altitudes" << std::endl;
    fix_sea_elevations(gen);
    std::vector<std::vector<int16_t> > river_elevations(elevations);
    surface_set surfaces = get_water_set(gen);
    int iterations = 0;
    while(!run_correction_iteration(surfaces))
    {
        iterations ++;
    }
    
    std::cout << "  --  Collecting rivers that created big altitude differences " << std::endl;\
    {
        river_mark_callback cb = { &river_elevations, &ignore_set, river_correction_limit };
        {
            water_generator::callback_guard guard(cb,gen);
            water_properties rv = rv_prop;
            rv.start_points = ignore_set;
            gen.load_rivers(rivers,rv);
        }
    }
    std::cout << "  -- Found " << ignore_set.size() << " rivers causing elevation drops above " <<  river_correction_limit << " feet"<< std::endl;
    std::cout << "  -- Saving removed_rivers.bmp... " << std::flush;
    {
        std::vector<std::vector<int16_t> > removed_rivers(map_size*2,std::vector<int16_t>(map_size*2,0));
        removed_river_callback cb = { &removed_rivers, &ignore_set};
        {
            water_generator::callback_guard guard(cb,gen);
            gen.load_rivers(rivers,rv_prop);
        }
        altitude_to_bmp(output_dir + "/removed_rivers.bmp", removed_rivers);
    }
    std::cout << "done" << std::endl;
    
    elevations.swap(fully_saved);
    rv_prop.start_points = ignore_set;
}
*/

std::string test_up_to_date(int check_updates)
{
    std::string res;
	try {
		makeahmap_version current = get_current_version();
        res = current.str();
        if(check_updates < 0)
            return res;
        time_t last_test = 0;
        {
            std::ifstream f("./data/last_update.txt");
            if(f) {
                f >> last_test;
                if(f && (time(0) - last_test < check_updates * 3600 * 24))
                    return res;
            }
        }
    	std::cout << "- Checking for makeahmap updates... " << std::flush;
		makeahmap_version latest = get_latest_version();
		if(current < latest) {
			std::cout << "  Newer Version Avalible " << latest << " !!!"<< std::endl;
            std::cout << "  ======================================================================" << std::endl;
            std::cout << "  ==                                                                  ==" << std::endl;
			std::cout << "  >>   Please use update_makeahmap.exe  to update to latest version   <<" << std::endl;
            std::cout << "  ==                                                                  ==" << std::endl;
            std::cout << "  ======================================================================" << std::endl;
		}
        else {
			std::cout << current << " is up to date " << std::endl;
        }
	}
	catch(std::exception const &e) {
		std::cout << " failed to get latest version number" << std::endl;
	}
	std::ofstream f("./data/last_update.txt");
	f << time(0);
    return res;
		
}

struct cluster_data {
    int size;
    int row;
    int col;
    bool operator<(cluster_data const &other) const
    {
        return size < other.size;
    }
};

cluster_data find_cluster(image<int> &img,int r,int c)
{
    int rows = img.height();
    int cols = img.width();
    double rs=0;
    double cs=0;
    cluster_data d=cluster_data();
    d.size = 1;
    rs+=r;
    cs+=c;

    std::queue<std::pair<int,int> > q;
    q.push(std::make_pair(r,c));
    static const int marked = 128;
    img[r][c]=marked;
    while(!q.empty()) {
        auto p=q.front();
        q.pop();
        for(int dr=-1;dr<=1;dr++) {
            for(int dc=-1;dc<=1;dc++) {
                int new_r=p.first+dr;
                int new_c=p.second+dc;
                if(new_r < 0 || new_c < 0 || new_r>=rows || new_c >= cols || img[new_r][new_c]==marked || img[new_r][new_c]==0)
                    continue;
                q.push(std::make_pair(new_r,new_c));
                rs+=new_r;
                cs+=new_c;
                d.size++;
                img[new_r][new_c]=marked;
            }
        }
    }
    d.row = rs / d.size;
    d.col = cs / d.size;
    return d;
}

void mark_clusters()
{
    std::priority_queue<cluster_data> clusters;
    int rows = types.height();
    int cols = types.width();
    for(int r=0;r<rows;r++) {
        for(int c=0;c<cols;c++) {
            if(types[r][c]==1) {
                cluster_data d=find_cluster(types,r,c);
                clusters.push(d);
            }
        }
    }
    int total = 10;
    while(!clusters.empty() && total > 0) {
        cluster_data d= clusters.top();
        clusters.pop();
        double lat = (lat2 - lat1) / rows * (rows - d.row) + lat1;
        double lon = (lon2 - lon1) / cols * d.col + lon1;
        char buf[256];
        snprintf(buf,sizeof(buf),"%3.3f,%3.3f",lat,lon);
        if(d.row > 0 && d.row < rows-1 && d.col > 0 && d.col < cols - 1) {
            types[d.row][d.col-1]=255;
            types[d.row][d.col]=255;
            types[d.row][d.col+1]=255;
            types[d.row-1][d.col]=255;
            types[d.row+1][d.col]=255;
        }

        print_str(buf,d.row,d.col,255,types);
        snprintf(buf,sizeof(buf),"%10.6f%% %9.3f,%9.3f",d.size*100.0 / rows / cols,lat,lon);
        std::cout << buf << std::endl;
        total --;
    }
}

void mark_type_and_save(int type)
{
    int rows = types.height();
    int cols = types.width();
    for(int r=0;r<rows;r++)
        for(int c=0;c<cols;c++)
            if(types[r][c]==type)
                types[r][c]=1;
            else
                types[r][c]=0;
    mark_clusters();
    write_reference_bmp();
}

void fix_periferal_elev(std::vector<std::vector<int16_t> > &elev)
{
    int rows = elev.size();
    int cols = elev[0].size();
    
    int size = 2*(rows-2+cols);

    auto get_r = [&](int pos) {
        pos %= size;
        if(pos < cols)
            return 0;
        pos -= cols;
        if(pos < rows - 2)
            return pos + 1;
        pos -= rows - 2;
        if(pos < cols)
            return rows-1;
        pos -= cols;
        return rows - 2 - pos;
    };

    auto get_c = [&](int pos) {
        pos %= size;
        if(pos < cols)
            return pos;
        pos -= cols;
        if(pos < rows - 2)
            return cols - 1;
        pos -= rows - 2;
        if(pos < cols)
            return cols-1 - pos;
        return 0;
    };
    int start = -1;
    for(int pos=0;pos<size;pos++) {
        if(elev[get_r(pos)][get_c(pos)] != dem::void_data) {
            start = pos;
            break;
        }
    }
    if(start == -1)
        throw std::runtime_error("Too much invalid DEM data - aborting");

    int16_t belev = 0;
    for(int n=0;n<2*(rows+cols);n++) {
        int begin = n + start;
        int r=get_r(begin);
        int c=get_c(begin);
        if(elev[r][c] != dem::void_data) {
            belev = elev[r][c];
            continue;
        }
        int end = begin;
        int16_t eelev;
        while((eelev=elev[get_r(end+1)][get_c(end+1)])==dem::void_data)
            end++;
        int d=end - begin + 2;
        int d2 = d/2;
        for(int pos = begin;pos <= end;pos++) {
            int w1 = d - (pos - begin + 1);
            int w2 = pos - begin + 1;
            int value = (belev * w1 +  eelev * w2 + d2) / d;
            int r = get_r(pos),c=get_c(pos);
            elev[r][c] = value;
        }
        n += d - 2;
    }
}

void fix_elevations(std::vector<std::vector<int16_t> > &elev,int vertices)
{
    std::cout << "-- Fixing DEM missing data for " << vertices << " points... " << std::flush;
    auto local_opt = solver_options;
    local_opt.force_cpu = true;
    std::unique_ptr<surface_solver_base> solver = get_solver(local_opt);
    int rows = elev.size();
    int cols = elev[0].size();
    
    fix_periferal_elev(elev);

    std::vector<std::vector<char > >  bmask(rows,std::vector<char >(cols));
    std::vector<std::vector<float> >  felev(rows,std::vector<float>(cols));

    int16_t max_val = 0;
    for(int r=0;r<rows;r++) {
        for(int c=0;c<cols;c++) {
            max_val = std::max(max_val,elev[r][c]);
            bmask[r][c] = elev[r][c] != dem::void_data;
            if(!bmask[r][c])
                felev[r][c]=0;
            else
                felev[r][c]=elev[r][c];
        }
    }
    surface_solver_base::stats st = solver->run(bmask,felev,0.0001);
    for(int r=0;r<rows;r++) {
        for(int c=0;c<cols;c++) {
            if(!bmask[r][c])
                elev[r][c] = int(felev[r][c] + 0.5f);
        }
    }
    std::cout << "Done in " << st.time << std::endl;
}

int main(int argc,char **argv)
{
    std::ofstream log("./output/log.txt");
    logger lcout(std::cout,log);
    logger lcerr(std::cerr,log);
    try {
        int save_globcover_type = -1;
        std::string file_name = "config.ini";

        if(argc > 2)
            throw std::runtime_error("Usage makeahmap [ /path/to/config.ini ]");
        else if(argc == 2) {
            std::string param=argv[1];
            if(param.find("--type=")==0) {
                save_globcover_type = atoi(param.c_str() + 7);
            }
            else {
                file_name = param;
            }
        }

        load_profile(file_name,log);
		
		log << "- Version: " << test_up_to_date(check_updates) << std::endl;

        downloader::manager::instance().init(download_sources,temp_dir,auto_download_enabled,disable_ssl_check);
       
        std::cout << "- Latitude, longitude range and scale" << std::endl;
        std::cout << std::setprecision(3) << std::fixed;
        std::cout << "    Lat: " << std::setw(10) << lat1 << ' ' << std::setw(10) << lat2 << std::endl;
        std::cout << "    Lon: " << std::setw(10) << lon1 << ' ' << std::setw(10) << lon2 << std::endl;
        std::cout << "  Scale: " << std::setw(10) << map_scale << std::endl;
        
        load_custom_type_mapping();

        std::cout << "- Loading GlobCover Data. " << std::endl;
        load_globcover_data();
        resample_type();

        if(save_globcover_type != -1) {
            mark_type_and_save(save_globcover_type);
            return 0;
        }
 
        std::cout << "- Loading Digital Elevations Model data. " << std::endl;
        int void_data_count=0;
        elevations = std::move(dem::read(db_type,map_size*8,lat1,lat2,lon1,lon2,void_data_count));
        
        if(void_data_count > 0)
            fix_elevations(elevations,void_data_count);
        
        water_generator gen(lat1,lat2,lon1,lon2,map_size * 8);
        std::cout << "- Loading & processing shores data. " << std::endl;
        gen.load_land(shores,elevations,lake_or_island_min_size,water_alt_limit);
        
        if(rv_prop.max_level != 0) {
            std::cout << "- Loading & processing rivers data... " << std::endl;
            gen.load_rivers(rivers,elevations,rv_prop,water_alt_limit);
            std::cout << "  Complete" << std::endl;
        }
        
        std::cout << "- Updateing altitudes... " << std::endl;
        gen.update_elevations(elevations,water_to_land_slope,water_to_land_range,solver_options);
        

        std::cout << "- Fixing ground types... " << std::flush;
        write_reference_bmp();
        recolor();
		double time = low_pass_filter_splattype(splattype_filter_radius);
        std::cout << "Done in " << time << " s" << std::endl;
        // */
        // make_beaches();
        
        
        std::cout << "- Saving ground types: splattype.bmp... " << std::flush;
        write_gndtype();
        std::cout << "Done" << std::endl;
        
        /*
        std::cout << "- Fixing elevations near sea water... " << std::flush;
        fix_sea_elevations(gen);
        std::cout << "Done" << std::endl;
        
        if(fix_river_slopes) {
            fix_river_elevations(gen);
        }
        */
        
        std::cout << "- Saving elevation data .raw and .bmp... " << std::flush;
        save_height_file();
        unsigned max_alt = altitude_to_bmp(output_dir + "/" + map_name + "_elevations.bmp",elevations);
        std::cout << "Done" << std::endl;
        std::cout << "-- Maximal altitude (for use with TE bmp import) is " << max_alt << " feet" << std::endl;
        
        std::cout << "- Generating clibboard map... " << std::flush;
        make_clipboard_map(max_alt,elevations);
        std::cout << "Done" << std::endl;
        
        std::cout << "\n\nCompleted\n";

        lcerr.close();
        lcout.close();
#if defined(_WIN32) || defined(WIN32)
        std::cout << "Press Enter to exit..." << std::endl;
        std::cin.get();
#endif
        return 0;
    }
    catch(std::exception const &e) {
        std::cerr << "\nError:" << e.what() << std::endl;
        lcerr.close();
        lcout.close();
#if defined(_WIN32) || defined(WIN32)
        std::cerr << "Press Enter to exit..." << std::endl;
        std::cin.get();
#endif        
        return 1;
    }
    return 0;
}

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
