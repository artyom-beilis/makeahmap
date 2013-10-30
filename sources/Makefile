../makeahmap: bmp.h  dem.h  fileio.h gshhs.h makeahmap.cpp downloader.h downloader.cpp 
	g++ -std=c++0x -Wall -Wextra -O2 -g makeahmap.cpp downloader.cpp -o ../makeahmap -lcurl -ltiff -lz

clean:
	rm -f ../makeahmap	
