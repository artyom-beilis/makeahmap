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

#include <stdint.h>

#pragma pack(push,1)
namespace bmp {
    struct rgbq {
        unsigned char r,g,b,res;
    };
    struct info_header {
        uint32_t size;
        uint32_t width;
        uint32_t height;
        uint16_t planes,bitcount;
        uint32_t compr;
        uint32_t size2;
        int32_t xppm,yppm;
        uint32_t clr_used;
        uint32_t clr_important;
        rgbq colors[256];

        info_header(int rows,int cols)
        {
            memset(this,0,sizeof(*this));
            size=40;
            width=cols;
            height=rows;
            planes=1;
            bitcount = 8;
            compr = 0;
            xppm=1000;
            yppm=1000;
            for(int i=0;i<256;i++) {
                colors[i].r=colors[i].g=colors[i].b=i;
            }
        }
    };
    struct header {
        char type[2];
        uint32_t size;
        uint16_t res1;
        uint16_t res2;
        uint32_t offset;
        info_header ih;
        header(int rows,int cols) : ih(rows,cols) {
            type[0]='B';
            type[1]='M';
            res1=res2=0;
            size = sizeof(*this) + rows * cols;
            offset = sizeof(*this); 
        }
    };
}
#pragma pack(pop)

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4


