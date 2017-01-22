#pragma once
#include <string>

namespace downloader {
    void download_file_generic(std::string url,std::string to,bool disable_ssl_check=true,bool progress_bar=true);
}

