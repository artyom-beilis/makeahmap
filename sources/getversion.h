#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <set>
#include <iostream>
#include <sstream>
#include "url.h"

#ifdef major
#undef major
#endif
#ifdef minor
#undef minor
#endif

#define BASE_URL std::string("http://cppcms.com/files/makeahmap/")


struct makeahmap_version;
inline std::ostream &operator<<(std::ostream &out,makeahmap_version const &v);

struct makeahmap_version {
	int major;
	int minor;
	int revision;
	std::string suffix;
	std::string str() const
	{
		std::ostringstream ss;
		ss << *this;
		return ss.str();
	}
	makeahmap_version(std::string const &v) : major(0), minor(0), revision(0) 
	{
		if(v.size() > 256)
			throw std::runtime_error("Invalid version " + v);
		char buf[256]={};
		if(sscanf(v.c_str(),"%d.%d.%d%s",&major,&minor,&revision,buf)==4
		   || sscanf(v.c_str(),"%d.%d.%d",&major,&minor,&revision)==3
		   || sscanf(v.c_str(),"%d.%d%s",&major,&minor,buf)==3
		   || sscanf(v.c_str(),"%d.%d",&major,&minor)==2)
		{
			suffix = buf;
			if(str()!=v) {
				throw std::runtime_error("Failed to parse version " + v);
			}
		}
		else
			throw std::runtime_error("Failed to parse version " + v);
	}
	bool operator==(makeahmap_version const &other) const
	{
		return !(*this<other) && !(other < *this);
	}
	bool operator<(makeahmap_version const &other) const
	{
		if(major != other.major)
			return major < other.major;
		if(minor != other.minor)
			return minor < other.minor;
		if(revision != other.revision)
			return revision < other.revision;
		return suffix < other.suffix;
	}
};

inline std::ostream &operator<<(std::ostream &out,makeahmap_version const &v)
{
	out << v.major << "." << v.minor;
	if(v.revision != 0)
		out << "." << v.revision;
	out<<v.suffix;
	return out;
}

inline makeahmap_version get_latest_version()
{
	std::string base = BASE_URL;
	std::string temp = "./temp/index.html";
	remove(temp.c_str());
	downloader::download_file_generic(base,temp,false,false);
	std::ifstream f(temp.c_str());
	if(!f) 
		throw std::runtime_error("Failed to open downloaded " + temp);
	std::string line;
	std::set<makeahmap_version> versions;
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
		versions.insert(makeahmap_version(ver));
	}
	f.close();
	remove(temp.c_str());
	if(versions.empty())
		throw std::runtime_error("No makeahmap versions found at " + base);
	makeahmap_version res = *versions.rbegin();
	remove(temp.c_str());
	return res;
}

inline makeahmap_version get_current_version()
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
	return makeahmap_version(ver);
}

