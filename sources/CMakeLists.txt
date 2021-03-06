cmake_minimum_required(VERSION 2.8)

# MinGW Specific
if(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_SIMD -mmmx -msse -msse2 -std=c++0x -O3 -g -Wall -m32 -Wl,--enable-stdcall-fixup -DDLL_EXPORT -Wno-deprecated-declarations")
	set(CMAKE_RC_FLAGS " -F pe-i386")
	set(PLATFORM_LIBS winmm wldap32 ws2_32 urlmon wininet)
	set(EXTRASRC makeahmap.rc)
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_SIMD -std=c++11 -O3 -g -Wall -Wno-deprecated-declarations -pthread")
	set(PLATFORM_LIBS pthread dl)
endif()

# Libraries

if(NOT WIN32)
	find_path(CURL_INC curl/curl.h)
	find_library(CURL curl)
endif()

find_path(ZLIB_INC zlib.h)
find_path(TIFF_INC tiff.h)
find_path(OCL_INC CL/cl.h)
 
find_library(TIFF tiff)
find_library(ZLIB z)
find_library(OCL OpenCL)

include_directories("${CURL_INC}")
include_directories("${ZLIB_INC}")
include_directories("${TIFF_INC}")

set(STDLIBS ${TIFF} ${ZLIB})

if(NOT WIN32)
	set(STDLIBS ${STDLIBS} ${CURL})
endif()

add_executable(makeahmap
		makeahmap.cpp
		downloader.cpp
		url.cpp
		surface.cpp
		text_to_image.cpp
		${EXTRASRC})


if(OCL_INC AND OCL)
	include_directories("${OCL_INC}")
	
	add_library(makeahmap_gpu SHARED gpu_solver.cpp)

	target_link_libraries(makeahmap_gpu ${OCL})

	add_executable(gpumemtest gpumemtest.cpp)
	target_link_libraries(gpumemtest ${OCL})
	

endif()

target_link_libraries(makeahmap ${STDLIBS} ${PLATFORM_LIBS})

if(WIN32)
	add_executable(update_makeahmap update_makeahmap.cpp url.cpp)
	set_target_properties(update_makeahmap PROPERTIES COMPILE_FLAGS "-static -static-libgcc -static-libstdc++")
	set_target_properties(update_makeahmap PROPERTIES LINK_FLAGS "-static -static-libgcc -static-libstdc++")
	target_link_libraries(update_makeahmap urlmon wininet)
endif()
