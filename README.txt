
                          M A K E A H M A P
        
            (C) 2016 Artyom Beilis <artyomtnk@yahoo.com>
        
                      For licensing details see:
                       copyright-makeahmap.txt


========================================================================
TEMPORARY NOTICE


Makeahmap version for AH3 - it is beta version

- Elevation with higher accuracy using SRTM3 - for newer model (vasly improves visual appearance)

Known - Issues

- It looks like globcover to splattype is suboptimal
- Terrains smaller than 512x512 had not been tested
- GTOPO30 does not work (source is gone and I can't find a new one)

========================================================================
This program is designed to prepare data for AH3 Maps for Terrain Editor.
========================================================================

It Generates
==============

Files for TE:
-------------

mapname.raw  - the SIGNED elevations file for import to TE
splattype.bmp  - ground types - for import in TE
mapname.bmp - clipboard map

Additional Files
----------------

ground_coverage.bmp - the ground types distribution - colored as created from GLOBCOVER database.
mapname_elevations.bmp - visaul representation of the elevation maps
sea_level.pgm - a map the represents visually sea level adjustments
water_bodies.pgm - location of water
flatten_terrain.pgn - location of the terrain areas that become below 0 after altitude adjustments

Required Data Sets
===================


+--------------------------------------------------------+
|NOTE: Starting from version 0.5 the required files are  |
|      downloaded automatically, no manual download      |
|      is required.                                      |
|                                                        |
|      You can still download them manually for offline  |
|      use                                               |
+--------------------------------------------------------+



1.  Digital Elevation Database. Any of the 3 below: GTOPO30, SRTM30, SRTM3, no need for all of them:

    (c) SRTM3 - 3 arc second database (~90 m resolution), it is actually the source for SRTM30,
        highly accurate but much bigger and it consists of many
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


    (b) SRTM30 - 30 arc second (1/2 nautical mile ~900m) recent database. It has similar format as
        GTOPO30 and the same size, but uses different data source.
		
		Advantage it allows to work for areas above 60N that isn't represented in SRTM3 
		for example murmansk area
	
    (c) NOT WORKING ANY MORE - online source gone
	
	    GTOPO30 - 30 arc second database (1/2 nautical mile ~900m), oldest dataset, but yet very good one.
    
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
   

	
2.  GlobCover Database - ground type mapping - describes which kind of ground is in use - dessert, corps, 
    forest etc.

    You can obtain it from there:
   
        http://due.esrin.esa.int/globcover/
   
    Technically you only need the file GLOBCOVER_L4_200901_200912_V2.3.tif which can be downloaded from:
    
        http://due.esrin.esa.int/globcover/LandCover2009/GLOBCOVER_L4_200901_200912_V2.3.color.tif
    
    Put it under:
    
        data/globcover/globcover.tif 
        
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

1. Configure correct elevations database srtm3 is best unless you are working on areas above N60
   than use srtm30
2. Configure the latitude, longitude and scale of the area you want to create
3. Make sure you setup the map size correctly and use correct map name

Save the config.ini and double click on the makeahmap.exe

  NOTE: if it fails to run telling that missing OpenCL.dll
  
  heck if your graphics card supports OpenCL (latest Intel, ATI and Nvidia do)
  make sure you have recent drivers
  
  Otherwise run makeahmap_cpu.exe - does the same but without GPU support
	
If everything ok.

Press enter when the program run is completerd to close the windows.
Don't forget to reviewing the logs the program prints.

Once the files are ready (output_dir directory in config.ini):

Go to terrain editor create a map of required size, save it and exit.

 
Open TE, and import signed raw altitude data and splattypes from the output 
directory, put the generated mymap.bmp - clipboard map to textures directory
for your TE.

Once it is complete, save it, exit and enter TE again.

Enjoy!!!
