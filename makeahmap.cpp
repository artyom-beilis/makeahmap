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
#include <iomanip>
#include "bmp.h"
#include "gshhs.h"

#include "dem.h"


std::vector<std::vector<uint16_t> > elevations;

int map_size = 512;
int river_width = 1;
int river_level = -1;
int lake_water_color = 0;
int river_water_color = 0;
int land_water_color = 0;

dem::db_properties db_type;

std::vector<std::vector<unsigned char> > bare_types;
std::vector<std::vector<unsigned char> > types;

static const int lower_lat = -65;

double lat1 = -1000,lat2 = -1000;
double lon1 = -1000,lon2 = -1000;

std::string map_name = "map";
std::string tiff_file="./data/globcover/GLOBCOVER_L4_200901_200912_V2.3.tif";
std::string shores = "./data/gshhs/gshhs_f.b";
std::string rivers = "./data/gshhs/wdb_rivers_f.b";
std::string custom_mapping;
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

void load_custom_type_mapping()
{
    prepare_map();
    if(custom_mapping.empty())
        return;
    std::cout <<"- Loading custom ground type mapping" << std::endl;
    std::ifstream csv(custom_mapping.c_str());
    if(!csv) {
        throw std::runtime_error("Failed to open file " + custom_mapping);
    }
    std::string line;
    int lineno=0;
    while(std::getline(csv,line)) {
        lineno++;
        std::istringstream ss(line);
        ss.unsetf(std::ios_base::basefield);
        unsigned gcover=0,ah=0;
        char c;

        ss >> gcover >> c >> ah;
        if(!ss || c!=',' || gcover >=256 || ah >= 256) {
            std::ostringstream tmp;
            tmp << "Invalid line " << lineno << " in file " << custom_mapping;
            throw std::runtime_error(tmp.str());
        }
        map_type[gcover]=ah;
    }
}

void load_profile(std::istream &in)
{
    std::string dem_prefix;
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
        else if(key=="type_mapping") {
            custom_mapping = value;
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
                throw std::runtime_error("Invalid database type name " + value + " valid are srtm3, srtm30 or gtopo30");
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
        else if(key == "river_level") {
            river_level = atoi(value.c_str());
        }
        else if(key == "river_water_color") {
            river_water_color = atoi(value.c_str());
        }
        else if(key == "lake_water_color") {
            lake_water_color = atoi(value.c_str());
        }
        else if(key == "land_water_color") {
            land_water_color = atoi(value.c_str());
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
        else {
            throw std::runtime_error("Unknown key " + key);
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
        double nmiles = scale * map_size * 1.60934 / 1.852;
        lat1=lat_c - (nmiles / 2 / 60 );
        lat2=lat_c + (nmiles / 2 / 60 );
        double diff = (nmiles / 2 / 60 ) / cos(lat_c / 180 * 3.14159);
        lon1=lon_c + diff;
        lon2=lon_c - diff;
    }
    if(db_type.rows == 0) {
        throw std::runtime_error("Undefined dem - should be one of srtm30, srtm3, gtopo30");
    }
    if(!dem_prefix.empty()) {
        db_type.directory = dem_prefix;
    }
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

unsigned write_alt_map()
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
    return max;
}


void update_gndtype(water_generator &gen)
{
    int gnd_size = map_size * 8;
    for(int r=0;r<gnd_size;r++) {
        for(int c=0;c<gnd_size;c++) {
            bool has_water = false;
            bool has_land = false;
            for(int dr = 0;dr<4;dr++) {
                for(int dc=0;dc<4;dc++) {
                    int v;
                    int water_r = r*4+dr;
                    int water_c = c*4+dc;
                    v = gen.pixel(water_r,water_c);
                    if(v == water_generator::land_mark)
                        has_land=true;
                    else 
                        has_water=true;
                }
            }
            int current_type = types[r][c];
            if(has_land && has_water) {
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
        memcpy(&bare_types[inr][0],&buf[c0],w);
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

void fix_sea_elevations(water_generator &gen,int pixels = 33)
{
    int max_water = map_size * 32;
    for(int r=0;r<map_size*2;r++) {
        for(int c=0;c<map_size*2;c++) {
            bool has_sea=false;
            for(int dr = -pixels;dr<=pixels;dr++) {
                for(int dc=-pixels;dc<pixels;dc++) {
                    int wr=r*16 + dr;
                    int wc=c*16 + dc;
                    if(wr < 0 || wr>= max_water || wc < c || wc >= max_water) {
                        continue;
                    }
                    if(gen.pixel(wr,wc)==water_generator::sea_mark) {
                        has_sea = true;
                        goto break_point;
                    }
                }
            }
        break_point:
            if(has_sea)
                elevations[r][c] = 0;
        }
    }
}

void make_map_color_index(bmp::header &hdr)
{
    // 200-249 - beach, river bad, rock
    // 150-199 - grass, grocky grass, sandy grass, swamp
    // 100-149 - snow
    //  50- 99 - forest
    //   0- 49 - farm 
    int colors[5][3] = {
        { 192, 192, 192 },
        { 192, 192, 192 },
        { 192, 192, 192 },
        { 192, 192, 192 },
        { 192, 192, 192 },
        //{ 90, 206, 59 },
        //{ 216, 186, 62 },
        //{ 192, 192, 192 },
        //{ 74, 255, 150 },
        //{ 255, 210, 76 }
    };
    for(int c=0;c<5;c++) {
        for(int i=0;i<50;i++) {
            double f=0.5 + i / 100.0;
            bmp::rgbq *p = &hdr.ih.colors[c*50 + i ];
            p->r = int(floor(colors[c][0] * f));
            p->g = int(floor(colors[c][1] * f));
            p->b = int(floor(colors[c][2] * f));
        }
    }
    bmp::rgbq *p = &hdr.ih.colors[255];
    p->r=40;
    p->g=44;
    p->b=191;
}

int get_color_from_type(int type,double emboss)
{
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
    
    // 255 - water
    // 254 - beach
    // 253 - grid line
    // 252 - grid text
    // 251 - TBD
    // 250 - TBD
    // 200-249 - beach, river bad, rock
    // 150-199 - grass, grocky grass, sandy grass, swamp
    // 100-149 - snow
    //  50- 99 - forest
    //   0- 49 - farm 
    
    int basic;
    switch(type & 0xF) {
    case 4:
    case 5:
        basic = 0;
        break;
    case 2:
    case 3:
        basic = 1;
        break;
    case 1:
    case 7:
    case 8:
    case 9:
        basic = 2;
        break;
    case 6:
    case 10:
    case 11:
        basic = 3;
        break;
    case 12:
        basic = 4;
        break;
    case 0:
        return 255;
    default:
        throw std::runtime_error("Internal error invalid ground type");
    }
    if(emboss < -1.0)
        emboss = -1.0;
    else if(emboss > 1.0)
        emboss = 1.0;
    int color = int(floor((emboss + 1.0) / 2 * 49)) + 50 * basic;
    return color;
}

void make_clipboard_map(int max_elev)
{
    int tsize = 1024;
    bmp::header h(tsize,tsize);
    make_map_color_index(h);
    std::string map_file = output_dir +"/" + map_name + ".bmp";
    FILE *f = fopen(map_file.c_str(),"wb");
    if(!f) {
        throw std::runtime_error("Failed to open " + map_file);
    }
    fwrite(&h,sizeof(h),1,f);
    int factor_type = map_size * 8;
    int factor_elev = map_size * 2;
    int elev_size = map_size * 2;
    for(int r=tsize-1;r>=0;r--) {
        std::vector<unsigned char> colors(tsize);
        for(int c=0;c<tsize;c++) {
                int t_r = factor_type * r / tsize;
                int t_c = factor_type * c / tsize;
                int type = types.at(t_r).at(t_c);
                
                int e_r = factor_elev * r / tsize;
                int e_c = factor_elev * c / tsize;
                double e[2] = {0,0};
                for(int d=-1,p=0;d<=1;d+=2,p++) {
                    if(e_r + d < 0 || e_r + d >= elev_size || e_c + d < 0 || e_c+d >= elev_size)
                        continue;
                    try {
                        e[p]=double(elevations.at(e_r+d).at(e_c+d))/max_elev;
                    }
                    catch(...) {
                        std::cerr << "Elev " << e_r <<" " << e_c << " " << d << std::endl;
                        throw;
                    }
                }
                double emboss = (e[1]-e[0]) * 5;
                colors[c] = get_color_from_type(type,emboss);
        }
        fwrite(&colors[0],1,tsize,f);
    }
    fclose(f);
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

       
        std::cout << "- Latitude and longitude range " << std::endl;
        std::cout << std::setprecision(3) << std::fixed;
        std::cout << "    Lat: " << std::setw(10) << lat1 << ' ' << std::setw(10) << lat2 << std::endl;
        std::cout << "    Lon: " << std::setw(10) << lon1 << ' ' << std::setw(10) << lon2 << std::endl;

        load_custom_type_mapping();

        std::cout << "- Loading GlobCover Data... " << std::flush;
        load_globcover_data();
        resample_type();
        std::cout << "Done" << std::endl;
        
        water_generator gen(lat1,lat2,lon1,lon2,map_size * 32);
        std::cout << "- Loading & processing shores data... " << std::flush;
        gen.load_land(shores);
        std::cout << "Done" << std::endl;
        
        if(river_level != 0) {
            std::cout << "- Loading & processing rivers data... " << std::flush;
            gen.load_rivers(rivers,river_width,river_level);
            gen.make_border();
            std::cout << "Done" << std::endl;
        }
        
        std::cout << "- Generating waterd.bmp... " << std::flush;
        gen.save_waterd_map(output_dir + "/waterd.bmp");
        std::cout << "Done" << std::endl;
        
        std::cout << "- Generating waterc.bmp... " << std::flush;
        gen.save_waterc_map(output_dir + "/waterc.bmp",lake_water_color,river_water_color,land_water_color);
        std::cout << "Done" << std::endl;

        std::cout << "- Fixing ground types according to shorelines shapes... " << std::flush;
        write_reference_bmp();
        recolor();
        update_gndtype(gen);
        std::cout << "Done" << std::endl;
        
        // make_beaches();
        
        std::cout << "- Saving ground types: gndtype.bmp... " << std::flush;
        write_gndtype();
        std::cout << "Done" << std::endl;

        std::cout << "- Loading Digital Elevations Model data... " << std::endl;
        std::vector<std::vector<uint16_t> > tmp=dem::read(db_type,map_size*2,lat1,lat2,lon1,lon2);
        elevations.swap(tmp);
        std::cout << "  DEM is ready" << std::endl;
        std::cout << "- Fixing elevations near sea water... " << std::flush;
        fix_sea_elevations(gen);
        std::cout << "Done" << std::endl;
        
        std::cout << "- Saving elevation data .elv and .bmp... " << std::flush;
        unsigned max_alt = write_alt_map();
        std::cout << "Done" << std::endl;
        std::cout << "-- Maximal altitude (for use with TE bmp import) is " << max_alt << " feet" << std::endl;
        
        std::cout << "- Generating clibboard map... " << std::flush;
        make_clipboard_map(max_alt);
        std::cout << "Done" << std::endl;
        
        std::cout << "\n\nCompleted\n";

#if defined(_WIN32) || defined(WIN32)
        std::cout << "Press Enter to exit..." << std::endl;
        std::cin.get();
        return 0;
#endif
    }
    catch(std::exception const &e) {
        std::cerr << "\nError:" << e.what() << std::endl;
        std::cerr << "Press Enter to exit..." << std::endl;
        std::cin.get();
        return 1;
    }
    return 0;
}

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
