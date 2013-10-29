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

#include <string>
#include <map>

namespace downloader {
    class manager {
    public:
        static manager &instance();
        void init(std::string profile,std::string temp_dir,bool enable_download);
        bool check(std::string real_file_name,std::string file_code,bool should_exist = true);
    private:

        static void download_file(std::string url,std::string to);
        static void unzip(std::string zip_file,std::string to_dir,std::string file1,std::string file2=std::string());
        static void gzip(std::string input,std::string output);
        static void untar(std::string tar_gz_file,std::string to_dir,std::string file_name);
        static std::string pretty(std::string url);
        static std::string dir_from_path(std::string full);
        static std::string name_from_path(std::string full);

        enum source_type { tiff, targz, zip };
        struct file_data {
            std::string url;
            source_type type;
            std::string file,file_extra;
        };
        
        std::map<std::string,file_data> profile_;
        std::string temp_dir_;
        bool enabled_;
    };

} // namespce

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

