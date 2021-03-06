#include "getversion.h"
#include <assert.h>
#include <stdio.h>
#if defined(WIN32)  || defined(_WIN32)
# include <windows.h>
# include <io.h>
#  ifndef access
#   define access(x,y) _access(x,y)
# endif
# ifndef F_OK
#  define F_OK 0
# endif
#else
# include <unistd.h>
#endif


void restore()
{
	rename("update_makeahmap.exe",".\\temp\\bad.exe");
	rename(".\\temp\\update_save.exe","update_makeahmap.exe");
	
}
bool checker()
{
	makeahmap_version current = get_current_version();
	std::cout << "Current Version : " << current << std::endl;
	makeahmap_version latest = get_latest_version();
	
	char const *latest_override = getenv("LATEST");
	if(latest_override)
		latest = makeahmap_version(latest_override);
	
	std::cout << "Latest Version  : " << latest << std::endl;
	if(current == latest || current > latest) {
		std::cout << "The version is up to date" << std::endl;
		return false;
	}
	std::string dir = ".\\temp\\makeahmap-" + latest.str();
	std::string url = BASE_URL + "makeahmap-" + latest.str() + ".zip";
	std::string zip = dir + ".zip";
	remove(zip.c_str());
	system("del /S /Q .\\temp\\makeahmap-*.zip .\\temp\\makeahmap-*");
	if(!latest_override) {
		std::cout << "Downloading: "  << url  << " to " << zip << std::endl;
		downloader::download_file_generic(url,zip);
	}
	else {
		CopyFile(("..\\makeahmap-" + latest.str() + ".zip").c_str(),zip.c_str(),FALSE);
	}
	if(system(("unzip -o " + zip + " -d .\\temp\\").c_str())!=0)
		throw std::runtime_error("Failed to unzip " + zip);	
	remove(".\\temp\\update_save.exe");
	if(rename("update_makeahmap.exe",".\\temp\\update_save.exe")!=0)
		throw std::runtime_error("Failed to move update_makeahmap.exe");
	std::cout << "Testing: " << (dir + "\\update_script.bat") << std::endl;
	if(access((dir + "\\update_script.bat").c_str(),F_OK)==0) {
		std::cout << "Running: " << (dir + "\\update_script.bat") << std::endl;
		if(system((dir + "\\update_script.bat").c_str())!=0) {
			restore();
			throw std::runtime_error(dir + "\\update_script.bat returned with error!");
		}
	}
	else {
		std::cout << "Using default upgrade procedure " << std::endl;
		if(system(("xcopy " + dir + "\\* .\\ /A /EXCLUDE:config.ini /E /Y").c_str())) {
			restore();
			throw std::runtime_error("Failed to copy files");
		}
		if(system(("copy /Y " + dir + ".\\config.ini .\\config_reference.ini").c_str())) {
			throw std::runtime_error("Faile to copy config.ini");
		}
		system("del /S /Q .\\temp\\makeahmap-*.zip .\\temp\\makeahmap-*");
	}
	current = get_current_version();
	if(current==latest) {
		std::cout << "Upgraded sucessefully to " << latest << std::endl;
		system("start notepad.exe Changelog.txt");
		return true;
	}
	else
		throw std::runtime_error("Failed to upgrade to " + latest.str() + " current version is " + current.str());
	return false;
}
int main()
{
	int r=0;
	try {
		if(checker())
			return 0;
	}
	catch(std::exception const &e) {
		r=1;
		std::cerr << "\n Error: "<< e.what() << std::endl;
	}
	std::cerr << "Press Enter to exit..." << std::endl;
	std::cin.get();
	return r;
}
