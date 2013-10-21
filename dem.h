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

#include <memory>
#include <string>
#include <stdio.h>
#include "fileio.h"
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

namespace dem {
    
    struct db_properties {
        int longitude_step;
        int latitude_step;
        int rows;
        int cols;
        int overlapping;
        bool top_left;
        bool may_be_missing;
        bool latitude_first;
        std::string directory;
        std::string suffix;
    };

    class tile {
    public:
        tile(db_properties const &prop,int column_no) 
        {
            row_ids_[0]=-1;
            row_ids_[1]=-1;
            file_id_ = -1;
            current_pos_ = 0;
            d=&prop;
            column_no_ = column_no;
            real_cols_ = prop.cols;
            if(prop.overlapping)
                real_cols_++;
        }

        int16_t get(int row,int col)
        {
            if(row_ids_[0] == row)
                return rows_[0][col];
            if(row_ids_[1] == row)
                return rows_[1][col];
            int file_id = row / d->rows;
            assert(file_id_ == -1 || file_id_ <= file_id);
            if(file_id != file_id_) {
                open_file(file_id);
                file_id_ = file_id;
            }
            int required_pos = row % d->rows;
            assert(required_pos >= current_pos_);
            while(current_pos_ < required_pos) {
                if(f_)
                    f_.skip(real_cols_*2);
                current_pos_++;
            }
            row_ids_[0] = row_ids_[1];
            rows_[0].swap(rows_[1]);
            std::vector<int16_t> new_row;
            
            new_row.swap(rows_[1]);
            new_row.resize(real_cols_);
            
            if(f_) {
                if(!f_.read(&new_row[0],real_cols_ * 2)) {
                    throw std::runtime_error("Failed to read data - unexpected EOF" +get_name(file_id));
                }
                for(int i=0;i<real_cols_;i++) {
                    uint16_t a=new_row[i];
                    uint16_t b = ((a >> 8) & 0xFF) | ((a << 8) & 0xFF00);
                    new_row[i]=b;
                }
                // fix void data
                for(int i=0;i<real_cols_;i++) {
                    static const int16_t void_data = -32768;
                    if(new_row[i]==void_data) {
                        if(i+1 < real_cols_ && new_row[i+1]!=void_data)
                            new_row[i]=new_row[i+1];
                        else if(i>0)
                            new_row[i]=new_row[i-1];
                        else if(row_ids_[0]!=-1) // valid row before
                            new_row[i]=rows_[0][i];
                        else
                            new_row[i]=0;
                    }
                }
                current_pos_++;
            }
            else
                memset(&new_row[0],0,real_cols_ * 2);

            row_ids_[1] = row;
            rows_[1].swap(new_row);
            return rows_[1][col];
        }

        std::string get_name(int file_id)
        {
            int lat = 90  -  file_id * d->latitude_step;
            if(!d->top_left)
                lat -= d->latitude_step;

            char lat_d = 'N';
            if(lat < 0) {
                lat = -lat;
                lat_d = 'S';
            }
            int lon = -180 + d->longitude_step * column_no_;
            char lon_d = 'E';
            if(lon < 0) {
                lon = -lon;
                lon_d = 'W';
            }
            std::ostringstream s_lat,s_lon;
            s_lat << lat_d << std::setw(2) << std::setfill('0') << lat;
            s_lon << lon_d << std::setw(3) << std::setfill('0') << lon;

            std::string file_name = d->directory + "/";
            if(d->latitude_first) {
                file_name += s_lat.str();
                file_name += s_lon.str();
            }
            else {
                file_name += s_lon.str();
                file_name += s_lat.str();
            }
            file_name += d->suffix;
            return file_name;
        }
        
    private:
        
        void open_file(int file_id)
        {
            f_.close();
           
            std::string file_name = get_name(file_id); 
            f_.open(file_name,!d->may_be_missing);
            if(!f_) {
                std::cerr << "\nFile " <<  file_name << " ... water?" << std::endl;
            }
            current_pos_ = 0;
        }

        db_properties const *d;

        std::vector<int16_t> rows_[2];
        int row_ids_[2];
        int file_id_;
        int column_no_;
        int current_pos_;
        int real_cols_;
        fileio f_;
    };

    std::vector<std::vector<uint16_t> > read(db_properties const &p,int points,double lat1,double lat2,double lon1,double lon2)
    {
        if(lat1 > lat2)
            std::swap(lat1,lat2);
        if(lon1 > lon2)
            std::swap(lon1,lon2);
        std::vector<std::vector<uint16_t> > elevations;
        
        elevations.resize(points,std::vector<uint16_t>(points,0));

        int tiles_no = 360 / p.latitude_step;
        std::vector<std::shared_ptr<tile> > tiles(tiles_no);
        for(int i=0;i<tiles_no;i++)
            tiles[i].reset(new tile(p,i));

        double lat_2_row = double(p.rows) / p.latitude_step;
        double lon_2_col = double(p.cols) / p.longitude_step;
    
        for(int r=0;r<points;r++) {
            double real_lat = (lat1 - lat2) / points * r + lat2; // from top to down
            double real_r = lat_2_row * (90-real_lat);
            int tile_r = int(floor(real_r));
            double rw0 = 1-(real_r - tile_r);
            double rw1 = 1-rw0;
            for(int c=0;c<points;c++) {
                double real_lon = (lon2 - lon1) / points * c + lon1;
                double real_c = lon_2_col * (real_lon + 180);
                int tile_c = int(floor(real_c));
                
                double cw0 = 1-(real_c - tile_c);
                double cw1 = 1-cw0;
                
                int tid0  = tile_c / p.cols;
                int toff0 = tile_c % p.cols;
                double c0r0 = tiles[tid0]->get(real_r,  toff0);
                double c0r1 = tiles[tid0]->get(real_r+1,toff0);

                int tid1  = (tile_c+1) / p.cols;
                int toff1 = (tile_c+1) % p.cols;
                double c1r0 = tiles[tid1]->get(real_r,  toff1);
                double c1r1 = tiles[tid1]->get(real_r+1,toff1);

                double v = (c0r0 * cw0 + c1r0 * cw1) * rw0 + (c0r1 * cw0 + c1r1 * cw1) * rw1;
                int iv = int(round(v * 3.28084));
                if(iv < 0)
                    iv = 0;
                elevations[r][c]=iv;
            }
        }
        return elevations;
    }

    inline db_properties srtm3()
    {
        db_properties p;
        p.longitude_step = 1;
        p.latitude_step = 1;
        p.rows = 1200;
        p.cols = 1200;
        p.overlapping = true;
        p.may_be_missing = true;
        p.latitude_first = true;
        p.top_left = false;
        p.directory = "./data/srtm3";
        p.suffix = ".hgt";
        return p;
    }

    inline db_properties dem_30()
    {
        db_properties p;
        p.longitude_step = 40;
        p.latitude_step = 50;
        p.rows = 6000;
        p.cols = 4800;
        p.overlapping = false;
        p.may_be_missing = false;
        p.latitude_first = false;
        p.top_left = true;
        p.suffix = ".DEM";
        return p;
    }
    inline db_properties srtm30()
    {
        db_properties p=dem_30();
        p.directory="./data/srtm30";
        return p;
    }
    inline db_properties gtopo30()
    {
        db_properties p=dem_30();
        p.directory="./data/gtopo30";
        return p;
    }
}


// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

