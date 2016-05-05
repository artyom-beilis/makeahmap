set BP=C:\Users\artik\mingw\packages\bin
mkdir release-%1
mkdir release-%1\output
mkdir release-%1\data
mkdir release-%1\data\globcover
mkdir release-%1\data\gshhs
mkdir release-%1\data\srtm3
mkdir release-%1\data\srtm30
mkdir release-%1\data\gtopo30
mkdir release-%1\temp
mkdir release-%1\images
copy %BP%\unzip.exe release-%1\
copy %BP%\jpeg62.dll release-%1\
copy %BP%\libtiff3.dll release-%1\
copy %BP%\zlib1.dll release-%1\
copy copyright*.txt release-%1\
copy images\*.tif release-%1\images
copy images\*.txt release-%1\images
copy download_sources.txt release-%1\ 
copy build\makeahmap.exe release-%1\ 
copy config.ini release-%1\ 
copy groundmapping.csv release-%1\
copy KNOWN-ISSUES.txt release-%1\
copy README.txt release-%1\ 
copy Changelog.txt release-%1\