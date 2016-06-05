cmake_minimum_required(VERSION 2.8)

# MinGW Specific
if(WINDOWS)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O2 -g -Wall -m32 -DDLL_EXPORT -DCURL_STATICLIB -Wno-deprecated-declarations")
	set(CMAKE_RC_FLAGS " -F pe-i386")
	set(PLATFORM_LIBS winmm wldap32 ws2_32)
	set(EXTRASRC makeahmap.rc)
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -g -Wall -Wno-deprecated-declarations")
	set(PLATFORM_LIBS )
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

set(STDLIBS ${CURL} ${TIFF} ${ZLIB})

if(OCL_INC AND OCL)
	include_directories("${OCL_INC}")
	
	add_executable(makeahmap_ocl
		makeahmap.cpp
		downloader.cpp
		surface.cpp
		${EXTRASRC}
		)
	target_link_libraries(makeahmap_ocl ${STDLIBS} ${OCL} ${PLATFORM_LIBS})

endif()

add_executable(makeahmap
		makeahmap.cpp
		downloader.cpp
		surface.cpp
		${EXTRASRC})

set_target_properties(makeahmap PROPERTIES COMPILE_FLAGS "-DUSE_CPU")
target_link_libraries(makeahmap ${STDLIBS} ${PLATFORM_LIBS})


