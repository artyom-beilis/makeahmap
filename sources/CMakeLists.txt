cmake_minimum_required(VERSION 2.8)

# MinGW Specific
if(WINDOWS)
	set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -m32 -DDLL_EXPORT -DCURL_STATICLIB -Wno-deprecated-declarations")
	set(CMAKE_RC_FLAGS " -F pe-i386")
else()
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wno-deprecated-declarations")
endif()

# Libraries

find_path(CURL_INC curl/curl.h)
find_path(ZLIB_INC zlib.h)
find_path(TIFF_INC tiff.h)
find_path(OCL_INC CL/cl.h)

find_library(CURL curl)
find_library(TIFF tiff)
find_library(ZLIB z)
find_library(OCL OpenCL)

include_directories("${CURL_INC}")
include_directories("${ZLIB_INC}")
include_directories("${TIFF_INC}")

set(EXTRALIBS )
if(OCL_INC AND OCL)
	include_directories("${OCL_INC}")
	set(EXTRALIBS ${OCL})
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_CPU")
endif()

if(WINDOWS)
	set(EXTRASRC makeahmap.rc)
endif()

add_executable(makeahmap
	makeahmap.cpp
	downloader.cpp
	surface.cpp
	${EXTRASRC}
	)

if(WINDOWS)
	target_link_libraries(makeahmap  ${CURL} winmm wldap32 ws2_32 ${TIFF} ${ZLIB} ${EXTRALIBS})
else()
	target_link_libraries(makeahmap  ${CURL}  ${TIFF} ${ZLIB} ${EXTRALIBS})
endif()
