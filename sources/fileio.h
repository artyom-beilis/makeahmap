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
#include <errno.h>
#include <string.h>

class outfile {
    outfile(outfile const &);
    void operator=(outfile const &);
public:
    outfile() : f_(0) {}
    void open(std::string const &name) 
    {
        close();
        name_ = name;
        f_ = fopen(name.c_str(),"wb");
        if(!f_) {
            int err = errno;
            throw std::runtime_error("Failed to open file `" + name_ + "': " + strerror(err));
        }
    }
    outfile(std::string const &name) : f_(0)
    {
        open(name);
    }
    void write(void const *buf,size_t s,size_t n)
    {
        write(buf,s*n);
    }
    void write(void const *buf,size_t n)
    {
        if(!f_)
            throw std::runtime_error("Internal error: outfile - file is not open");
        if(fwrite(buf,1,n,f_)!=n) {
            int err = errno;
            throw std::runtime_error("Failed to write to file `" + name_ + "': " + strerror(err));
        }
    }
    void close()
    {
        if(!f_)
            return;
        if(fclose(f_)!=0) {
            int err = errno;
            f_=0;
            throw std::runtime_error("Failed to close to file `" + name_ + "': " + strerror(err));
        }
        f_ = 0;
    }
    ~outfile()
    {
        if(f_) {
            fclose(f_);
            f_=0;
        }
    }
private:
    FILE *f_;
    std::string name_;
};

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
            char buf[1024];
            while(n > 0) {
                size_t rem = sizeof(buf);
                if(rem > n)
                    rem = n;
                if(!read(buf,rem)) {
                    throw std::runtime_error("Internal error in gz stream skipping");
                }
                n-=rem;
            }
            /*if(gzseek(gzf_,n,SEEK_CUR) < 0) { fails on Windows WTF?
                throw std::runtime_error("Internal error gzseek failed\n");
            }*/
        }
    }
private:
    FILE *f_;
    gzFile gzf_;
};

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4


