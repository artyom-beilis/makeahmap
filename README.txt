
                          M A K E A H M A P
        
            (C) 2013 Artyom Beilis <artyomtnk@yahoo.com>
        
                      For licensing details see:
                       copyright-makeahmap.txt

========================================================================
This program is designed to prepare data for AH2 Maps for Terrain Editor.
========================================================================

It Generates
==============

Files for TE:
-------------

mapname.elv  - the elevations file
gndtype.bmp  - ground covering
waterd.bmp - import water bmp file (depth)
waterc.bmp - import water bmp file (water color)

Additional Files
----------------

ground_coverage.bmp - the ground types distribution - colored as created from GLOBCOVER database.
mapname_elevations.bmp - visaul representation of the elevation maps
mapname.bmp - clipboard map (not fully useful yet as lacks greed and finer colors)

Required Data Sets
===================


1.  Digital Elevation Database. Any of the 3 below: GTOPO30, SRTM30, SRTM3, no need for all of them:

    (a) GTOPO30 - 30 arc second database (1/2 nautical mile ~900m), oldest dataset, but yet very good one.
    
        Its resolution is similar to the resolution AH uses (1/2 a mile ~800m), so it should fill
        all map makers needs.
   
        Download from: ftp://edcftp.cr.usgs.gov/data/gtopo30/global/ files named like eXXXnZZ.tar.gz
        
        You can extract the files using 7zip. You will only need DEM files from the archives.
        Once all the files are extracted, put them under data/gtopo30:
        
            ...
            data/gtopo30/E020N40.DEM
            data/gtopo30/E020N90.DEM
            ...

        
        Note: You don't need all the files, you can detect wich one you need using this image:
        
        ftp://edcftp.cr.usgs.gov/data/gtopo30/global/tiles.gif
   
    (b) SRTM30 - 30 arc second (1/2 nautical mile ~900m) recent database. It has similar format as
        GTOPO30 and the same size, but uses different data source.
        
        As GTOPO30, its resolution is similar to the resolution AH uses (1/2 a mile ~800m),
        so it should fill all map makers needs. In comparison to SRTM3 - it is much smaller,
        but has higher resolution. In fact SRTM30 is build from combination of SRMT3 and GTOPO30
   
        Download from: http://dds.cr.usgs.gov/srtm/version2_1/SRTM30/ files like eXXXnZZ.dem.zip
       
        Unzip them to data/srmt30 such that you have files like
       
            ...
            data/srtm30/E020N40.DEM
            data/srtm30/E020N90.DEM
            ...
    
        Note: You don't need all the files, you can detect wich one you need using this image:
        
        ftp://edcftp.cr.usgs.gov/data/gtopo30/global/tiles.gif
        
    (c) SRTM3 - 3 arc second database (~90 m resolution), it is actually the source for SRTM30,
        highly accurate but much bigger and most likly overkill for map makers. I consists of many
        smaller files describing 1x1 degree areas, you can detect wich one you need using this image:
        
        http://dds.cr.usgs.gov/srtm/version2_1/Documentation/Continent_def.gif
        
        The files themselves can be downloaded from: http://dds.cr.usgs.gov/srtm/version2_1/SRTM3/
        
        Once you unzip them you should place them under data/srtm3 directory like this:
        
            ...
            data/srtm3/N36E034.hgt
            data/srtm3/N36E035.hgt
            ...

    Note: the data/gtopo30, data/srtm30 and data/srtm3 are placeholders for covinience, you can put the
    files at any directory on the disk and set the dem_path option in the config.ini that points to them.
            
2.  GlobCover Database - ground type mapping - describes which kind of ground is in use - dessert, corps, 
    forest etc.

    You can obtain it from there:
   
        http://due.esrin.esa.int/globcover/
   
    Technically you only need the file GLOBCOVER_L4_200901_200912_V2.3.tif which can be downloaded from:
    
        http://due.esrin.esa.int/globcover/LandCover2009/GLOBCOVER_L4_200901_200912_V2.3.color.tif
    
    Put it under:
    
        data/globcover/GLOBCOVER_L4_200901_200912_V2.3.color.tif 
        
    Or any other location, pointed by globcover_tiff_path config.ini option.
   
3.  GSHHS Database - Shoreline and river database, the provides high resolution lines for the shores and rivers
    that needed for generation of waterd.bmp/waterc.bmp

    You can download them from there: http://www.soest.hawaii.edu/pwessel/gshhg/
 
    Note: you need to download BINARY database fromat: gshhg-bin-2.2.3.zip that can be downloaded
    from te following exact URL (the current latest version):
    
        ftp://ftp.soest.hawaii.edu/pwessel/gshhg/gshhg-bin-2.2.3.zip
   
    You need only two files gshhs_f.b and wdb_rivers_f.b from the zip file, put them under
    
        data/gshhs/gshhs_f.b
        data/gshhs/wdb_rivers_f.b

    As usual, you can provide altetnative pathes using 
   
Usege
=====

Edit config.ini file:

1. Configure correct elevations database gtopo30, srtm30 or srtm3
2. Configure the latitude, longitude and scale of the area you want to create
3. Make sure you setup the map size correctly and use correct map name

Save the config.ini and double click on the makeahmap.exe, if everything ok.
Press enter when the program run is completerd to close the windows.
Don't forget to reviewing the logs the program prints.

Once the files are ready (output_dir directory in config.ini):

Go to terrain editor create a map of required size, save it and exit.

Copy the files: mapname.elv gndtype.bmp  to the required directory.
Copy waterd.bmp and waterc.bmp to the imptexpt, for example

  ...\ahiiterr\mymap\imptexpt\waterc.bmp
  ...\ahiiterr\mymap\imptexpt\waterd.bmp
  ...\ahiiterr\mymap\gndtype.bmp
  ...\ahiiterr\mymap\mymap.elv
  
You can also put mymap.bmp to the textsrc directory, but it is not 
fully useful yet.
  
Open TE, and run "import water bitmap"

Once it is complete, save it, exit and enter TE again.

Enjoy!!!

Building 
========

If you want to change the code, and rebuild the program you need libtiff and zlib you build it:

    g++ -std=c++0x -Wall -DDLL_EXPORT -O2 -g -I path/to/libtiff/include makeahmap.cpp -L path/to/libtiff/lib -ltiff -lz

MSVC may work as well.


Saving Disk Space
=================

The elevation DEM files of SRTM30 and GTOPO30, hgt files of SRTM3 and ".b" files of GSHHS database
can be stored compressed using gzip compression - note not ZIP, but gzip. For example:

    data/srtm30/E020N40.DEM   -> data/srtm30/E020N40.DEM.gz
    data/srtm3/N36E035.hgt    -> data/srtm3/N36E035.hgt.gz
    data/gshhs/wdb_rivers_f.b -> data/gshhs/wdb_rivers_f.b.gz

This would allow significantly reduce the reqired disk space by
all database files.
