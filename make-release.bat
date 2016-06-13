set BP=C:\Users\artik\mingw\packages\bin
set DIST=C:\Users\artik\mingw\packages\dist
mkdir makeahmap-%1
mkdir makeahmap-%1\output
mkdir makeahmap-%1\data
mkdir makeahmap-%1\data\globcover
mkdir makeahmap-%1\data\gshhs
mkdir makeahmap-%1\data\srtm3
mkdir makeahmap-%1\data\srtm30
mkdir makeahmap-%1\data\gtopo30
mkdir makeahmap-%1\temp
mkdir makeahmap-%1\images
copy %BP%\unzip.exe makeahmap-%1\
copy %BP%\jpeg62.dll makeahmap-%1\
copy %BP%\libtiff3.dll makeahmap-%1\
copy %BP%\zlib1.dll makeahmap-%1\
copy %DIST%\*.dll makeahmap-%1\
copy copyright*.txt makeahmap-%1\
copy images\*.tif makeahmap-%1\images
copy images\*.txt makeahmap-%1\images
copy download_sources.txt makeahmap-%1\ 
copy build\makeahmap.exe makeahmap-%1\ 
copy build\makeahmap_cpu.exe makeahmap-%1\ 
copy config.ini makeahmap-%1\ 
copy groundmapping.csv makeahmap-%1\
copy KNOWN-ISSUES.txt makeahmap-%1\
copy README.txt makeahmap-%1\ 
copy Changelog.txt makeahmap-%1\