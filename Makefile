makeahmap: bmp.h  dem.h  fileio.h  groundcover.h  gshhs.h makeahmap.cpp
	g++ -std=c++0x -Wall -Wextra -O2 -g makeahmap.cpp -o makeahmap -ltiff -lz

clean:
	rm -f makeahmap	
