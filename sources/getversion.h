#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <set>
#include <iostream>
#include "url.h"

inline std::string get_current_version()
{
	std::ifstream f("Changelog.txt");
	if(!f)
		throw std::runtime_error("Failed to open Changelog.txt");
	std::string vline;
	std::getline(f,vline);
	if(!f)
		throw std::runtime_error("Failed to read Changelog.txt");
	std::string prefix="Version ";
	std::string ver=vline.substr(vline.find(prefix) + prefix.size());
	while(!ver.empty() && (ver[ver.size()-1]==' ' || ver[ver.size()-1]=='\t' || ver[ver.size()-1]=='\r' || ver[ver.size()-1]=='\n'))
		ver.resize(ver.size()-1);
	if(ver.empty() || ver[0] < '0' || ver[0] > '9')
		throw std::runtime_error("Invalid version " + ver + " in Changelog.txt");
	return ver;
}

inline std::string get_latest_version()
{
	std::string base = "http://cppcms.com/files/makeahmap/";
	std::string temp = "./temp/index.html";
	remove(temp.c_str());
	downloader::download_file_generic(base,temp,false,false);
	std::ifstream f(temp.c_str());
	if(!f) 
		throw std::runtime_error("Failed to open downloaded " + temp);
	std::string line;
	std::set<std::string> versions;
	std::string prefix = "href=\"makeahmap-";
	std::string suffix = ".zip\">makeahmap-";
	while(std::getline(f,line)) {
		size_t start = line.find(prefix);
		if(start==std::string::npos)
			continue;
		line = line.substr(start + prefix.size());
		size_t size = line.find(suffix);
		if(size==std::string::npos)
			continue;
		std::string ver = line.substr(0,size);
		versions.insert(ver);
	}
	if(versions.empty())
		throw std::runtime_error("No makeahmap versions found at " + base);
	std::string res = *versions.rbegin();
	remove(temp.c_str());
	return res;
}
