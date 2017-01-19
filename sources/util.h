#pragma once
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

#include <dirent.h>
#include <string>
#include <set>
namespace util {

	inline std::set<std::string> dir(std::string const &path,std::string const &ext=std::string()) 
	{
		std::set<std::string> r;
		std::string temp_dir = path + "/";
		DIR *d=opendir(temp_dir.c_str());
		if(!d){
			throw std::runtime_error("Failed to open directory " + path);
		}
		struct dirent *de=0;
		while((de=readdir(d))!=0) {
			if(de->d_name[0]=='.') 
				continue;
			std::string name = de->d_name;
			if(!ext.empty()) {
				if(name.size() < ext.size() || name.substr(name.size() - ext.size()) != ext)
					continue;
			}
			r.insert(std::string(de->d_name));
		}
		closedir(d);
		return r;
	}

} // namespce

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

