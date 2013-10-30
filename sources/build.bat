set CFLAGS=-std=c++0x -O3 -g -m32 -DDLL_EXPORT -DCURL_STATICLIB -Wall -Wextra -I E:\AH1-Terrain\code\packages\include
set LINK=-L E:\AH1-Terrain\code\packages\lib -lcurl -lwinmm -lwldap32 -lws2_32 -ltiff -lz

windres -F pe-i386 makeahmap.rc -O coff -o makeahmap.res
g++ -c %CFLAGS% makeahmap.cpp -o makeahmap.o
g++ -c %CFLAGS% downloader.cpp -o downloader.o
g++ -m32 makeahmap.o downloader.o makeahmap.res %LINK% -o ../makeahmap.exe

del makeahmap.o
del downloader.o
del makeahmap.res
