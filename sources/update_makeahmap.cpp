#include "getversion.h"

int main()
{
	try {
		std::cout << "Upgrading from [" << get_current_version() << "] to [" << get_latest_version() << "]" << std::endl;
	}
	catch(std::exception const &e) {
		std::cerr << e.what() << std::endl;
	}
}
