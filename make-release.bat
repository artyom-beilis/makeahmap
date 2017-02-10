set BP=C:\Users\artik\mingw\packages\bin
set DIST=C:\Users\artik\mingw\packages\dist
del /Q makeahmap-%1.zip
del /Q /S makeahmap-%1
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
mkdir makeahmap-%1\resources\
mkdir makeahmap-%1\resources\groundmapping
copy %BP%\unzip.exe makeahmap-%1\
copy %BP%\jpeg62.dll makeahmap-%1\
copy %BP%\libtiff3.dll makeahmap-%1\
copy %BP%\zlib1.dll makeahmap-%1\
copy %DIST%\*.dll makeahmap-%1\
copy resources\groundmapping\*.csv makeahmap-%1\resources\groundmapping
copy copyright*.txt makeahmap-%1\
copy images\*.tif makeahmap-%1\images
copy images\*.txt makeahmap-%1\images
copy download_sources.txt makeahmap-%1\ 
copy build\makeahmap.exe makeahmap-%1\ 
copy build\libmakeahmap_gpu.dll makeahmap-%1\ 
copy build\update_makeahmap.exe makeahmap-%1\ 
copy config.ini makeahmap-%1\ 
copy KNOWN-ISSUES.txt makeahmap-%1\
copy README.txt makeahmap-%1\ 
copy Changelog.txt makeahmap-%1\
copy cleanup_script.bat makeahmap-%1\
echo set VER=%1>makeahmap-%1\update_script.bat
type update_script.bat >>makeahmap-%1\update_script.bat
C:\msys64\usr\bin\zip.exe makeahmap-%1.zip -r makeahmap-%1
