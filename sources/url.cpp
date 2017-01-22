#include "url.h"
#if defined WIN32 && ! defined USE_CURL_ON_WINDOWS
#include <urlmon.h>
#else
#include <curl/curl.h>
#endif
#include <stdexcept>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>

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

	void download_file_generic(std::string url,std::string to,bool,bool progress)
	{
		winprogress st;
		winprogress *st_p = progress ? &st : nullptr;
		if(URLDownloadToFile(NULL,url.c_str(),to.c_str(),0,st_p)!=S_OK) {
			throw std::runtime_error("Failed to download file " + url + " to " + to);
		}
		if(progress)
			st.done();
	}
	
	#else
	
    void download_file_generic(std::string url,std::string to,bool disable_ssl_check,bool progress)
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
        if(disable_ssl_check)
            curl_easy_setopt(curl,CURLOPT_SSL_VERIFYPEER,0l);
	if(progress) {
        	curl_easy_setopt(curl,CURLOPT_NOPROGRESS,0l);
        	curl_easy_setopt(curl,CURLOPT_PROGRESSDATA,&percent);
	        curl_easy_setopt(curl,CURLOPT_PROGRESSFUNCTION,downloader_progress_function);   
	}
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
        if(progress && percent!=100)
            downloader_progress_function(&percent,1,1,0,0);
    }

#endif
} // namespace

