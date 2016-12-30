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

#include "downloader.h"

#if defined WIN32 && ! defined USE_CURL_ON_WINDOWS
#include <urlmon.h>
#else
#include <curl/curl.h>
#endif
#include <zlib.h>
#include <stdexcept>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <libgen.h>
#include <dirent.h>





namespace downloader {
    extern "C" {
        inline size_t downloader_write_data(void *ptr, size_t size, size_t nmemb, FILE *stream) 
        {
            return fwrite(ptr, size, nmemb, stream);
        }
        inline int downloader_progress_function(void *clientp, double dltotal, double dlnow, double /*ultotal*/, double /*ulnow*/)
        {
            int *last_value = static_cast<int*>(clientp);
            int percent;
            if(dltotal != 0)
                percent =  int(100 * dlnow / dltotal);
            else
                percent = -2;
            if(*last_value != percent) {
                if(*last_value!=-1)
                    std::cout << "\b\b\b\b\b\b\b\b";
                if(percent>=0)
                    std::cout << std::setw(6) << percent << "% " << std::flush;
                else {
                    char const *prefixes="BKMG";
                    int p=0;
                    while(dlnow > 1000 && prefixes[p+1]!=0 ) {
                        dlnow/=1000;
                        p++;
                    }
                    std::cout << std::setw(6) << std::fixed << std::setprecision(1) << dlnow << prefixes[p] << ' ' << std::flush;
                }
                *last_value = percent;
            }
            return 0;
        }
    }

	#if defined WIN32 && !defined USE_CURL_ON_WINDOWS
	
	class winprogress : public IBindStatusCallback  
	{
	public:

		winprogress() : percent_(-1) {}
		
		void done()
		{
			if(percent_!=100)
				downloader_progress_function(&percent_,1,1,0,0);
		}
		
		// IUnknown
		STDMETHOD_(ULONG,AddRef)()  { return 0; }
		STDMETHOD(QueryInterface)(REFIID,void **) { return E_NOTIMPL; }
		STDMETHOD_(ULONG,Release)() { return 0; }
		// IBindStatusCallback		

		STDMETHOD(GetBindInfo)(DWORD *,BINDINFO *) { return E_NOTIMPL; }
		STDMETHOD(GetPriority)(LONG *) { return E_NOTIMPL; }
		STDMETHOD(OnDataAvailable)(DWORD,DWORD,FORMATETC *,STGMEDIUM *) { return E_NOTIMPL; }
		STDMETHOD(OnLowResource)(DWORD) { return E_NOTIMPL; }
		STDMETHOD(OnObjectAvailable)(REFIID,IUnknown *)  { return E_NOTIMPL; }

		STDMETHOD(OnProgress)(ULONG ulProgress,ULONG ulProgressMax,ULONG,LPCWSTR)
		{
			downloader_progress_function(&percent_,ulProgressMax,ulProgress,0,0);
			return S_OK;
		}
		
		STDMETHOD(OnStartBinding)(DWORD , IBinding *) { return E_NOTIMPL; }
		STDMETHOD(OnStopBinding)(HRESULT,LPCWSTR){ return E_NOTIMPL; }
		
	private:
		int percent_;
	};

	void manager::download_file(std::string url,std::string to)
	{
		winprogress st;
		if(URLDownloadToFile(NULL,url.c_str(),to.c_str(),0,&st)!=S_OK) {
			throw std::runtime_error("Failed to download file " + url + " to " + to);
		}
		st.done();
	}
	
	#else
	
    void manager::download_file(std::string url,std::string to)
    {
        std::string part = to + ".part";
        FILE *fout = fopen(part.c_str(),"wb");
        if(!fout)
            throw std::runtime_error("Failed to open " + to);
        
        CURL *curl = curl_easy_init();
        int percent = -1;
        curl_easy_setopt(curl,CURLOPT_URL,url.c_str());
        
        curl_easy_setopt(curl,CURLOPT_WRITEDATA,fout);
        curl_easy_setopt(curl,CURLOPT_WRITEFUNCTION,downloader_write_data);
        curl_easy_setopt(curl,CURLOPT_NOPROGRESS,0l);
        curl_easy_setopt(curl,CURLOPT_PROGRESSDATA,&percent);
        curl_easy_setopt(curl,CURLOPT_PROGRESSFUNCTION,downloader_progress_function);   
        CURLcode res = curl_easy_perform(curl);
        bool write_fail = ferror(fout) || fclose(fout) != 0;
        bool download_fail = res != 0;
        curl_easy_cleanup(curl);
        curl = 0;
        bool rename_fail = rename(part.c_str(),to.c_str()) != 0;
        if(download_fail || write_fail || rename_fail) {
            remove(part.c_str());
        }
        if(download_fail) 
            throw std::runtime_error("Failed to download file " + url + " " + curl_easy_strerror(res));
        if(write_fail) 
            throw std::runtime_error("Failed to save to file " + to);
        if(rename_fail) 
            throw std::runtime_error("Failed rename file " + part + " to " + to);
        if(percent!=100)
            downloader_progress_function(&percent,1,1,0,0);
    }
	
	#endif

    void manager::untar(std::string tar_gz_file,std::string to_dir,std::string file_name)
    {
        std::string full_out = to_dir + "/" + file_name + ".gz";
        std::string part_out = full_out + ".part";
        gzFile in=0,out=0;
        try {
            in = gzopen(tar_gz_file.c_str(),"rb");
            if(!in)
                throw std::runtime_error("Failed to open " + tar_gz_file);
            static const int block_size = 512;
            char buf[block_size];
            char blank[block_size]={0};
            while(gzread(in,buf,sizeof(buf))==int(sizeof(buf))) {
                if(memcmp(buf,blank,sizeof(buf))==0)
                    break;
                char slen[13]={0};
                strncpy(slen,buf+124,12);
                long long len = 0,padded_len;
                #if defined(WIN32) || defined(_WIN32)
                sscanf(slen,"%I64o",&len);
                #else
                sscanf(slen,"%llo",&len);
                #endif
                padded_len = (len + block_size -1 )/block_size*block_size;
                if(strncmp(buf,file_name.c_str(),100)==0) {
                    out = gzopen(part_out.c_str(),"wb");
                    if(!out)
                        throw std::runtime_error("Failed to open " + part_out);
                    while(len > 0) {
                        if(gzread(in,buf,block_size)!=block_size) {
                            throw std::runtime_error("Failed to read from " + tar_gz_file + " unexpected EOF");
                        }
                        int chunk = len;
                        if(chunk > block_size)
                            chunk = block_size;
                        if(gzwrite(out,buf,chunk)<0)
                            throw std::runtime_error("Failed to write to " + part_out);
                        len -= chunk;
                    }
                    if(gzclose(out)!=0) {
                        out = 0;
                        throw std::runtime_error("Failed to write to " + part_out);
                    }
                    gzclose(in);
                    in = 0;
                    if(rename(part_out.c_str(),full_out.c_str())!=0) 
                        throw std::runtime_error("Failed to rename file " + part_out + " to " + full_out);
                    return;
                }
                while(padded_len > 0) {
                    if(gzread(in,buf,block_size)!=block_size) 
                        throw std::runtime_error("Failed to read from " + tar_gz_file + " unexpected EOF");
                    padded_len -= block_size;
                }
            }
            throw std::runtime_error("Couldn't find " + file_name + " in " + tar_gz_file);
        }
        catch(...) {
            if(in)
                gzclose(in);
            if(out)
                gzclose(out);
            remove(part_out.c_str());
            throw;
        }
    }
   
    std::string manager::file_target(std::string file)
    {
        std::string target;
        file_target(file,target);
        return target;
    }
    void manager::file_target(std::string &file,std::string &target)
    {
        size_t pos;
        if((pos=file.find("->"))!=std::string::npos) {
            target = file.substr(pos+2);
            file = file.substr(0,pos);
        }
        else {
            target = file;
        }
    }
    void manager::unzip(std::string zip_file,std::string to_dir,std::string file1,std::string file2/*=std::string()*/)
    {
        std::string command = "unzip -qo \""+zip_file+"\" \"" + file1 + "\" ";
        if(!file2.empty())
            command += " \"" + file2 + "\"";
        command += " -d \"" + to_dir + "\"";
        if(system(command.c_str())!=0) {
            remove((to_dir + "/" + file1).c_str());
            if(!file2.empty()) {
                remove((to_dir + "/" + file2).c_str());
                file2 = ", " + file2;
            }
            throw std::runtime_error("Failed to unzip file(s) " + file1 + file2 + " from archive " + zip_file);
        }
    }
    
    void manager::gzip(std::string input,std::string output)
    {
        FILE *in=0;
        gzFile out=0;
        std::string part = output+".part";
        output=output+".gz";
        try {
            in=fopen(input.c_str(),"rb");
            if(!in) 
                throw std::runtime_error("Failed to open " + input);
            out = gzopen(part.c_str(),"wb");
            if(!out) 
                throw std::runtime_error("Failed to open " + output);
            char buf[1024];
            size_t n;
            while((n = fread(buf,1,sizeof(buf),in)) > 0) {
                if(gzwrite(out,buf,n)!=int(n)) 
                    throw std::runtime_error("Failed to write to " + output);
            }
            fclose(in);
            in = 0;
            if(gzclose(out)!=0) 
                    throw std::runtime_error("Failed to write to " + output);
            out = 0;
            if(rename(part.c_str(),output.c_str())!=0)
                throw std::runtime_error("Failed to rename file " + part + " to " + output );
        }
        catch(...) {
            if(in)
                fclose(in);
            if(out)
                gzclose(out);
            remove(part.c_str());
            throw;
        }
    }

    std::string manager::dir_from_path(std::string full)
    {
        std::vector<char> name(full.begin(),full.end());
        name.push_back('\0');
        char *dir = dirname(&name[0]);
        if(!dir)
            return "";
        return std::string(dir);
    }
    std::string manager::name_from_path(std::string full)
    {
        std::vector<char> name(full.begin(),full.end());
        name.push_back('\0');
        char *s = basename(&name[0]);
        if(!s)
            return full;
        return std::string(s);
    }

    std::string manager::pretty(std::string url)
    {
        if(url.size() < 40)
            return url;
        std::string first = url.substr(0,17);
        std::string second = url.substr(url.size() - 18);
        return first + "[...]" + second;
    }

    bool manager::check(std::string real_file_name,std::string file_code,bool should_exist)
    {
        auto p = profile_.find(file_code);
        if(p==profile_.end()) {
            if(!should_exist)
                return false;
            throw std::runtime_error("Missing file " + file_code + " in downloader profile");
        }
        if(!enabled_)
            return true;

        if(real_file_name.size()>=3 && real_file_name.substr(real_file_name.size()-3)==".gz")
            real_file_name = real_file_name.substr(0,real_file_name.size()-3);

        std::string gzname = real_file_name + ".gz";

        file_data d = p->second;
        if(access(real_file_name.c_str(),R_OK)==0)
            return true;
        if(d.type != tiff && access(gzname.c_str(),R_OK)==0)
            return true;
        
        std::cout<<"  - Downloading " << d.url << std::endl;
        std::cout<<"      downloaded:" << std::flush;
        
        if(d.type == tiff) {
            download_file(d.url,real_file_name);
            std::cout<<std::endl;
            return true;
        }

        std::string dir = dir_from_path(real_file_name);
        std::string expected = name_from_path(real_file_name);
        if(expected != file_target(d.file) && expected!=file_target(d.file_extra))
            throw std::runtime_error("Can't download file, internal error: expected name " 
                    + expected + " but actual file to download is " 
                    + file_target(d.file) + (d.file_extra.empty() ? std::string() : " or " + file_target(d.file_extra)));
        
        std::string zip_file;

        switch(d.type) {
        case targz:
            zip_file = temp_dir_ + "/tmp.tar.gz";
            download_file(d.url,zip_file);
            std::cout<<std::endl;
            std::cout<<"  - extracting and compressing... " <<std::flush;
            untar(zip_file,dir,d.file);
            remove(zip_file.c_str());
            std::cout<<"complete" << std::endl;
            break;

        case zip:
            {
                zip_file = temp_dir_ + "/tmp.zip";
                download_file(d.url,zip_file);
                std::cout<<std::endl;
                std::cout<<"  - extracting... " <<std::flush; 
                std::string tgt1,tgt2;
                std::string src1=d.file;
                std::string src2=d.file_extra;
                file_target(src1,tgt1);
                file_target(src2,tgt2);
                unzip(zip_file,temp_dir_,src1,src2);
                std::cout << "complete\n" << std::flush;
                if(tgt1 == src1) {
                    std::cout<<"  - compressing " << d.file << "..." << std::flush;
                    gzip(temp_dir_ + "/" + d.file,dir + "/" + d.file);
                    remove((temp_dir_ +"/" + d.file).c_str());
                    std::cout<<"done" << std::endl;
                }
                else {
                    rename((temp_dir_ + "/" + src1).c_str(),(dir + "/" + tgt1).c_str());
                }
                if(!d.file_extra.empty()) {
                    if(tgt2 == src2) {
                        std::cout<<"  - compressing " << d.file_extra << "..." << std::flush;
                        gzip(temp_dir_ + "/" + d.file_extra,dir + "/" + d.file_extra);
                        remove((temp_dir_ +"/" + d.file_extra).c_str());
                        std::cout<<"done" << std::endl;
                    }
                    else {
                        rename((temp_dir_ + "/" + src2).c_str(),(dir + "/" + tgt2).c_str());
                    }
                }
                remove(zip_file.c_str());
                std::cout<<"complete" << std::endl;
            }
            break;
        default:
            ; // should not get there 
        }
        return true;
    }

    void manager::init(std::string profile,std::string temp_dir,bool enabled)
    {
        temp_dir_ = temp_dir;
        enabled_ = enabled;
        profile_.clear();

        std::ifstream in(profile.c_str());
        if(!in) 
            throw std::runtime_error("Failed to load download manager profile");
        std::string line;
        int count=0;
		
		DIR *d=opendir(temp_dir_.c_str());
		if(!d){
			throw std::runtime_error("Failed to open temporary directory :" + temp_dir_);
		}
		struct dirent *de=0;
		bool message_written = false;
		while((de=readdir(d))!=0) {
			if(de->d_name[0]=='.') 
				continue;
			if(!message_written) {
				std::cout << "Cleaning temporary directory:\n";
				message_written = true;
			}
			if(remove((temp_dir_ + "/" + de->d_name).c_str())!=0) {
				std::string msg = strerror(errno);
				std::cout << "- Warning, failed to remove:" << de->d_name << ", " << msg << '\n';
			}
			else {
				std::cout << "- Removed: " << de->d_name << '\n';
			}
		}
		closedir(d);
		
        while(std::getline(in,line)) {
            count++;
            if(line.empty() || line[0]=='#' || line.find_first_not_of(" \t\r\n")==std::string::npos)
                continue;
            
            std::istringstream ss(line);

            std::ostringstream errline;
            errline << "Invalid line " << count << " in download manager profile " << profile << ": ";

            std::string key,type;
            file_data d;
            ss >> key >>  type >> d.file >> d.file_extra >> d.url;
            if(!ss)
                throw std::runtime_error(errline.str());
            if( type == "tiff" && d.file_extra == "-")
                d.type=tiff;
            else if(type == "tgz" && d.file_extra == "-")
                d.type=targz;
            else if(type == "zip")
                d.type=zip;
            else 
                throw std::runtime_error(errline.str() + " invalid type " + type + " or too many files for type");
            if(profile_.find(key)!=profile_.end())
                throw std::runtime_error(errline.str() + " duplicate key " + key);
            if(d.file_extra=="-")
                d.file_extra="";
            profile_[key]=d;
        }
    }

    manager &manager::instance()
    {
        static manager m;
        return m;
    }





} // namespce

// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

