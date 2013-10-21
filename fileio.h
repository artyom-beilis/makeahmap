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

#include <stdio.h>
#include <zlib.h>
#include <stdexcept>
#include <string>

class fileio {
    fileio(fileio const &);
    void operator=(fileio const &);
public:
    void open(std::string fin,bool throw_on_error = true)
    {
        close();
        std::string gzfile;
        std::string file;
        if(fin.size() > 3 && fin.substr(fin.size()-3) == ".gz") {
            file = fin.substr(0,fin.size()-3);
            gzfile = fin;
        }
        else {
            gzfile = fin + ".gz";
            file = fin;
        }
        gzf_ = gzopen(gzfile.c_str(),"rb");
        if(!gzf_) {
            f_ = fopen(file.c_str(),"rb");
        }
        if(throw_on_error && !gzf_ && !f_) {
            throw std::runtime_error("Failed to open neither " + gzfile + " nor " + file);
        }
    }
    fileio() : f_(0), gzf_(0)
    {
    }
    fileio(std::string file,bool throw_on_error = true) : f_(0), gzf_(0)
    {
        open(file,throw_on_error);
    }
    void close()
    {
        if(f_) {
            fclose(f_);
            f_ = 0;
        }
        if(gzf_) {
            gzclose(gzf_);
            gzf_ = 0;
        }
    }
    ~fileio()
    {
        close();
    }
    operator bool() const
    {
        return f_ || gzf_;
    }
    bool read(void *ptr,size_t n)
    {
        if(f_) {
            if(fread(ptr,1,n,f_)!=n)
                return false;
            return true;
        }
        if(gzf_) {
            int r = gzread(gzf_,ptr,n);
            if(r<=0)
                return false;
            if(size_t(r) != n)
                return false;
            return true;
        }
        return false;
    }
    void skip(size_t n)
    {
        if(f_)
            if(fseek(f_,n,SEEK_CUR) < 0) {
                throw std::runtime_error("Internal error fseek failed\n");
            }
        if(gzf_) {
            if(gzseek(gzf_,n,SEEK_CUR) < 0) {
                throw std::runtime_error("Internal error gzseek failed\n");
            }
        }
    }
private:
    FILE *f_;
    gzFile gzf_;
};

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4


