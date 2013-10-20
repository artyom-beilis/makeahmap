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
#include "bmp.h"
#include "gshhs.h"


std::vector<std::vector<uint16_t> > elevations;

int map_size = 512;
int river_width = 1;

std::vector<std::vector<unsigned char> > bare_types;
std::vector<std::vector<unsigned char> > types;

static const int lower_lat = -65;

int width;
int height;
int r0,r1;
int c0,c1;
double lat1 = -1000,lat2 = -1000;
double lon1 = -1000,lon2 = -1000;

std::string map_name = "map";
std::string tiff_file="./data/globcover/GLOBCOVER_L4_200901_200912_V2.3.tif";
std::string dem_prefix = "./data/srtm";
std::string shores = "./data/gshhs/gshhs_f.b";
std::string rivers = "./data/gshhs/wdb_rivers_f.b";
std::string output_dir = "./output";
char const *real_file=0;
char const *color_file=0;

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

static const int map_type_in[][3] = {
 { 11,  0x44 }, //Post-flooding or irrigated croplands (or aquatic)
 { 14,  0x55 }, //Rainfed croplands
 { 20,  0x24 }, //Mosaic cropland (50-70%) / vegetation (grassland/shrubland/forest) (20-50%)
 { 30,  0x42 }, //Mosaic vegetation (grassland/shrubland/forest) (50-70%) / cropland (20-50%) 
 { 40,  0x21 }, //Closed to open (>15%) broadleaved evergreen or semi-deciduous forest (>5m)
 { 50,  0x22 }, //Closed (>40%) broadleaved deciduous forest (>5m)
 { 60,  0x21 }, //Open (15-40%) broadleaved deciduous forest/woodland (>5m)
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

bool is_water(int c) { return c==210 || c==230; }

static unsigned char map_type[256];

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

void prepare_map()
{
    for(int i=0;map_in[i][0]!=-1;i++) {
        for(int j=0;j<3;j++) 
            tmap[map_in[i][0]][j]=map_in[i][j+1];
    }
    for(int i=0;map_type_in[i][0]!=-1;i++) {
        map_type[map_type_in[i][0]] = map_type_in[i][1];
    }
}


void load_profile(std::istream &in)
{
    bool via_scale = false;
    bool via_coord = false;
    double scale = -1;
    double lat_c=-1000,lon_c=-1000;
    int line = 0;
    while(!in.eof()) {
        std::string s;
        std::getline(in,s);
        line++;
        if(s[0]=='#')
            continue;
        bool found =false;
        for(size_t i=0;i<s.size() && !found;i++) {
            if(s[i]!=' ' && s[i] != '\t' && s[i] != '\r')
                found = true;
        }
        if(!found)
            continue;
        std::istringstream ss(s);
        std::string key,value;
        ss >> key >> value;
        if(!ss) {
            std::ostringstream ss;
            ss << "Invalid entry in confiuguration file line " << line;
            throw std::runtime_error(ss.str());
        }
        if(key == "map_size") {
            map_size = atoi(value.c_str());
            switch(map_size) {
            case 64:
            case 128:
            case 256:
            case 512:
                break;
            default:
                throw std::runtime_error("Invalid map size " + value);
            }
        }
        else if(key == "scale") {
            via_scale = true;
            scale = atof(value.c_str());
        }
        else if(key == "map_name") {
            map_name = value;
        }
        else if(key == "shores") {
            shores = value;
        }
        else if(key == "rivers") {
            rivers = value;
        }
        else if(key == "river_width") {
            river_width = atoi(value.c_str());
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
        else if(key == "lat1" || key == "lat2" || key=="lat") {
            double v = atof(value.c_str());
            if(!(lower_lat < v && v < 90)) {
                std::ostringstream ss;
                ss << "Invalid latitude " << value << " should be in rage [" << lower_lat << ",90]";
                throw std::runtime_error(ss.str());
            }
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
                throw std::runtime_error("Invalid longitude " + value + " should be in rage [-180,180]");
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
    }
    if(via_scale == via_coord) {
        throw std::runtime_error("Either lat/lon/scale or lat1/lat2/lon1/lon2 should be specified");
    }
    if(via_coord) {
        if(lat1 == -1000 || lat2 == -1000 || lon1 == -1000 || lon2 == -1000) {
            throw std::runtime_error("The latitude or longitude range and not fully defined");
        }
    }
    else { // via scale
        if(lat_c == -1000 || lon_c == -1000 || scale == -1) {
            throw std::runtime_error("The latitude, longitude or scane are not fully defined");
        }
        double nmiles = map_size * 1.60934 / 1.852;
        lat1=lat_c - (nmiles / 2 / 60 );
        lat2=lat_c + (nmiles / 2 / 60 );
        double diff = (nmiles / 2 / 60 ) / cos(lat_c / 180 * 3.14159);
        lon1=lon_c + diff;
        lon2=lon_c - diff;
        std::cout << lat1 << ' ' << lat2 << std::endl;
        std::cout << lon1 << ' ' << lon2 << std::endl;
    }
}



TIFF *init_image()
{
    TIFFSetWarningHandler(0);
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
    width = imw;
    height = imh;
    return in;
}

FILE *init_file(std::string const &name,char const *type,int w=0,int h=0)
{
    FILE *f=fopen(name.c_str(),"wb");
    if(!f) {
        throw std::runtime_error("Failed to open "+name);
    }
    if(w==0)
        w=r1-r0+1;
    if(h==0)
        h=c1-c0+1;
    fprintf(f,"%s\n%d %d\n255\n",type,h,w);
    return f;
}

int row_from_lat(double lat)
{
    return int((90 - lat)/(90 - lower_lat) * height);
}

int col_from_lon(double lon)
{
    return int((lon + 180) / 360 * width);
}
void calc_rows()
{
    r0=row_from_lat(lat1);
    r1=row_from_lat(lat2);
    if(r1 < r0)
        std::swap(r1,r0);
    c0=col_from_lon(lon1);
    c1=col_from_lon(lon2);
    if(c1 < c0)
        std::swap(c1,c0);
}

void resample_type()
{
    // resample
    int bmp_size = map_size * 8;
    types.resize(bmp_size,std::vector<unsigned char>(bmp_size,0));
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
    FILE *f=fopen(fname.c_str(),"wb");
    if(!f)  {
        throw std::runtime_error("Failed to open file " + fname);
    }
    int rows = types.size();
    int cols = types[0].size();
    bmp::header hdr(rows,cols);

    for(int i=0;map_in[i][0]!=-1;i++) {
        bmp::rgbq *c = &hdr.ih.colors[map_in[i][0]];
        c->r = map_in[i][1];
        c->g = map_in[i][2];
        c->b = map_in[i][3];
    }

    fwrite(&hdr,sizeof(hdr),1,f);
    
    for(int i=types.size()-1;i>=0;i--) {
        if(fwrite(&types[i][0],1,types[i].size(),f)!=types[i].size()) {
            throw std::runtime_error("Failed to write to file " +fname);
        }
    }
    if(fclose(f)!=0) {
        throw std::runtime_error("Failed to save file " +fname);
    }
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

void make_beaches()
{
    int bmp_size = map_size * 8;
    // make beaches... where possible
    for(int r=1;r<bmp_size-1;r++) {
        for(int c=r;c<bmp_size-1;c++) {
            if(types[r][c]==0) { // water
                for(int dr=-1;dr<=1;dr++) {
                    for(int dc=-1;dc<=1;dc++) {
                        if(types[r+dr][c+dc]!=0) // is not water
                            types[r+dr][c+dc]=beach_type;
                    }
                }
            }
        }
    }
}

void write_gndtype()
{
    std::string fname = output_dir + "/gndtype.bmp";
    FILE *f=fopen(fname.c_str(),"wb");
    if(!f)  {
        throw std::runtime_error("Failed to open file " + fname);
    }
    int rows = types.size();
    int cols = types[0].size();
    bmp::header hdr(rows,cols);
    fwrite(&hdr,sizeof(hdr),1,f);
    
    for(int i=types.size()-1;i>=0;i--) {
        if(fwrite(&types[i][0],1,types[i].size(),f)!=types[i].size()) {
            throw std::runtime_error("Failed to write to " + fname);
        }
    }
    if(fclose(f)!=0) 
        throw std::runtime_error("Failed to close " + fname);
}

void read_elevations()
{
    double lat_start  = lat1;
    double lat_end = lat2;
    double lon_start = lon1;
    double lon_end = lon2;
    if(lat_start < lat_end)
        std::swap(lat_start,lat_end);
    if(lon_start > lon_end)
        std::swap(lon_start,lon_end);
    
    int matrix_size = map_size * 2;
    
    elevations.resize(matrix_size,std::vector<uint16_t>(matrix_size,0));

    double dlat = (lat_end - lat_start) / matrix_size;
    double dlon = (lon_end - lon_start) / matrix_size;

    int f_start = int(floor(lon_start));
    int f_end =   int(floor(lon_end));
    int n_files = f_end - f_start + 1;
    static const int file_size = 1200;
    std::vector<std::vector<int16_t> > data(file_size + 1,std::vector<int16_t>(file_size * n_files + 1,0));

    int prev_lat_indx = 0;
    for(int r=0;r<matrix_size;r++) {
        double dr = r * dlat + lat_start;
        int req_lat_indx = int(floor(dr));
        double ddr = dr - req_lat_indx;
        if(r == 0 || prev_lat_indx != req_lat_indx) {
            prev_lat_indx = req_lat_indx;
            int pos = 0;
            for(int f=f_start;f<=f_end;f++,pos+=file_size) {
                char name[32];
                snprintf(name,sizeof(name),"%c%02ld%c%03ld",
                    ( req_lat_indx >= 0 ? 'N' : 'S' ),
                    labs(req_lat_indx),
                    ( f >= 0 ? 'E' : 'W' ),
                    labs(f) 
                );
                std::cout << "Loading " << name << "... " << std::flush;
                FILE *fin = fopen((dem_prefix+"/" + name + ".hgt").c_str(),"rb");
                if(!fin) {
                    std::cout << "missing, water?" << std::endl;
                    for(int j=0;j<file_size+1;j++) {
                        for(int k=0;k<file_size+1;k++) {
                            data[j][pos+k] = 0;
                        }
                    }
                    continue;
                }

                for(int j=0;j<file_size+1;j++) {
                    if(fread(&data[j][pos],2,file_size+1,fin)!=size_t(file_size+1)) {
                        fclose(fin);
                        throw std::runtime_error("Failed to read file - unexpected eof" + std::string(name));
                    }
                    for(int k=0;k<file_size+1;k++) {
                        // little endian
                        uint16_t a = data.at(j).at(pos+k);
                        uint16_t b = ((a >> 8) & 0xFF) | ((a << 8) & 0xFF00);
                        data.at(j).at(pos+k)=b;
                    }
                }
                fclose(fin);
                std::cout << "ok" << std::endl;
                
            }
        }
        for(int c=0;c<matrix_size;c++) {
            double dc = lon_start + c * dlon;
            int req_lon_indx = int(floor(dc));
            double ddc = dc - req_lon_indx;
            int cpos = (req_lon_indx - f_start) * file_size + int(round(ddc * file_size));
            int rpos = int(round((1.0-ddr) * file_size));
            
            int val = 0;
            val = data.at(rpos).at(cpos);
            if(val < 0)
                val = 0;
            double feet = val * 3.28084; //feet in meter
            uint16_t value = int(round(feet));
            elevations.at(r).at(c)=value;
        }
    }

}

void write_alt_map()
{
    unsigned max = 0;
    int matrix_size = elevations.size();
    for(int i=0;i<matrix_size;i++) {
        for(int j=0;j<matrix_size;j++) {
            if(elevations[i][j] > max)
                max = elevations[i][j];
        }
    }
    if(max == 0) {
        throw std::runtime_error("Maximal elevation is 0");
    }
    std::cout << "Maximal elevation is " << max << " feet" << std::endl;
    std::string elev_file = output_dir + "/" + map_name + "_elevations.bmp";
    FILE *f=fopen(elev_file.c_str(),"wb");
    if(!f) {
        throw std::runtime_error("Failed to open " + elev_file);
    }
    bmp::header hdr(matrix_size,matrix_size);
    fwrite(&hdr,sizeof(hdr),1,f);
    std::vector<uint8_t> row(matrix_size,0);
    for(int i=matrix_size-1;i>=0;i--) {
        for(int j=0;j<matrix_size;j++) {
            row[j] = elevations[i][j] * 255u / max;
        }
        if(fwrite(&row[0],1,matrix_size,f)!=size_t(matrix_size)) {
            fclose(f);
            throw std::runtime_error("Failed to write file to disk " + elev_file);
        }
    }
    if(fclose(f)!=0) {
        throw std::runtime_error("Failed to close file " + elev_file);
    }
    elev_file = output_dir +"/" + map_name + ".elv";
    f=fopen(elev_file.c_str(),"wb");
    if(!f) {
        throw std::runtime_error("Failed to open " + elev_file);
    }
    for(int i=matrix_size-1;i>=0;i--) {
        if(fwrite(&elevations[i][0],2,matrix_size,f)!=size_t(matrix_size)) {
            throw std::runtime_error("Failed to write file to disk " + elev_file);
        }
    }
    if(fclose(f)!=0) {
        throw std::runtime_error("Failed to close " + elev_file);
    }

}


void make_waterd()
{
    std::string fdname = output_dir + "/waterd.bmp";
    FILE *fd=fopen(fdname.c_str(),"wb");
    if(!fd) {
        throw std::runtime_error("Failed to open " + fdname);
    }
    std::string fcname = output_dir + "/waterc.bmp";
    FILE *fc=fopen(fcname.c_str(),"wb");
    if(!fc) {
        throw std::runtime_error("Failed to open " + fcname);
    }
    int water_size = map_size * 32;
    bmp::header hdr(water_size,water_size);
    fwrite(&hdr,sizeof(hdr),1,fd);
    fwrite(&hdr,sizeof(hdr),1,fc);

    std::vector<unsigned char> waterd(water_size,0);
    std::vector<unsigned char> waterc(water_size,0);
    
    int bh = bare_types.size();
    int bw = bare_types[0].size();
    int elev_size = map_size * 2;

    double r_factor = double(bh) / water_size;
    double c_factor = double(bw) / water_size;

    for(int r=water_size-1;r>=0;r--) {
        for(int c=0;c<water_size;c++) {
            double real_r = r*r_factor;
            double real_c = c*c_factor;
            int r0 = int(floor(real_r));
            int c0 = int(floor(real_c));
            double r0_w = 1.0 - (real_r - r0);
            double c0_w = 1.0 - (real_c - c0);
            double r1_w = 1.0 - r0_w;
            double c1_w = 1.0 - c0_w;
            int r1 = r0+1;
            int c1 = c0+1;
            if(r1 >=bh) r1=bh-1;
            if(c1 >=bw) c1=bw-1;
            double weight = 
                  c0_w * (is_water(bare_types[r0][c0]) * r0_w + is_water(bare_types[r1][c0]) * r1_w)
                + c1_w * (is_water(bare_types[r0][c1]) * r0_w + is_water(bare_types[r1][c1]) * r1_w);
            bool water_detected = weight >= 0.5;
            if(!water_detected) {
                waterd[c]=255;
                waterc[c]=0;
            }
            else {
                waterd[c]=0;
                waterc[c] = int(floor((1.0 - weight)*2 * 255));
                int elev_r = r / 16;
                int elev_c = c / 16;
                if(elevations[elev_r][elev_c] <= 300) {
                    // do now lower high altitude lakes
                    for(int ter=elev_r-1;ter<=elev_r+1;ter++) {
                        for(int tec=elev_c-1;tec<=elev_c+1;tec++) {
                            if(ter < 0 || ter >= elev_size || tec<0 || tec >= elev_size) 
                                continue;
                            elevations[ter][tec]=0;
                        }
                    }
                }
            }
        }
        if(fwrite(&waterd[0],1,water_size,fd)!=size_t(water_size)) {
            throw std::runtime_error("Failed to write to " + fdname);
        }
        if(fwrite(&waterc[0],1,water_size,fc)!=size_t(water_size)) {
            throw std::runtime_error("Failed to write to " + fcname);
        }
    }
    if(fclose(fd)!=0) {
        throw std::runtime_error("Failed to close " + fdname);
    }
    if(fclose(fc)!=0) {
        throw std::runtime_error("Failed to close " + fcname);
    }
}

void update_gndtype(water_generator &gen)
{
    int gnd_size = map_size * 8;
    std::cerr << gnd_size << " " << types.size() << std::endl;
    for(int r=0;r<gnd_size;r++) {
        for(int c=0;c<gnd_size;c++) {
            bool has_water = false;
            bool has_land = false;
            bool has_sea_or_lake = false;
            for(int dr = 0;dr<4;dr++) {
                for(int dc=0;dc<4;dc++) {
                    int v;
                    int water_r = r*4+dr;
                    int water_c = c*4+dc;
                    v = gen.pixel(water_r,water_c);
                    if(v == water_generator::land_mark)
                        has_land=true;
                    else {
                        has_water=true;
                        if(v==water_generator::sea_mark || v==water_generator::lake_mark)
                            has_sea_or_lake = true;
                    }
                }
            }
            int current_type = types[r][c];
            if(has_land && has_water) {
                    if(has_sea_or_lake || current_type == 0)
                        types[r][c] = beach_type;
                       
            }
            else if(has_land && current_type == 0) {
                types[r][c] = beach_type;
            }
            else if(has_water && current_type != 0) {
                types[r][c] = 0;
            }
        }
    }
}

int main(int argc,char **argv)
{
    try {
        std::string file_name = "config.ini";
        if(argc == 2)  {
            file_name = argv[1];
        }
        std::ifstream cfg(file_name.c_str());
        if(!cfg) {
            throw std::runtime_error("Failed to open file " + file_name);
        }
        load_profile(cfg);
        cfg.close();

        TIFF *in = init_image();
        if(!in)
            return 1;
        calc_rows();

        size_t len = TIFFScanlineSize(in);
        unsigned char *buf=new unsigned char[len];

        prepare_map();

        int h = r1-r0 + 1;
        int w = c1-c0 + 1;

        bare_types.resize(h,std::vector<unsigned char>(w,0));

        for(int r=r0,inr=0;r<=r1;r++,inr++) {
            TIFFReadScanline(in,buf,r);
            memcpy(&bare_types[inr][0],buf+c0,w);
        }
        TIFFClose(in);

        resample_type();
        write_reference_bmp();
        recolor();
        
        water_generator gen(lat1,lat2,lon1,lon2,map_size * 32);
        gen.load_land(shores);
        gen.load_rivers(rivers,river_width);
        gen.save_water_map(output_dir + "/waterd.bmp",true);
        gen.save_water_map(output_dir + "/waterc.bmp",false);
        update_gndtype(gen);
        // make_beaches();
        write_gndtype();
        read_elevations();
        // make_waterd();
        write_alt_map();
    }
    catch(std::exception const &e) {
        std::cerr << "Error:" << e.what() << std::endl;
        std::cerr << "Press Enter to exit..." << std::endl;
        std::cin.get();
        return 1;
    }
    return 0;
}

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
