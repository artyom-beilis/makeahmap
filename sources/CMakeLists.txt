cmake_minimum_required(VERSION 2.8)

# MinGW Specific
if(WINDOWS)
	set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -m32 -DDLL_EXPORT -DCURL_STATICLIB")
	set(CMAKE_RC_FLAGS " -F pe-i386")
else()
	set(CMAKE_CXX_FLAGS "-std=c++11 -Wall")
endif()

# Libraries

find_path(CURL_INC curl/curl.h)
find_path(ZLIB_INC zlib.h)
find_path(TIFF_INC tiff.h)

find_library(CURL curl)
find_library(TIFF tiff)
find_library(ZLIB z)

include_directories("${CURL_INC}")
include_directories("${ZLIB_INC}")
include_directories("${TIFF_INC}")

set(EXTRASRC )
if(WINDOWS)
	set(EXTRASRC makeahmap.rc)
endif()

add_executable(makeahmap
	makeahmap.cpp
	downloader.cpp
	${EXTRASRC}
	)

if(WINDOWS)
	target_link_libraries(makeahmap  ${CURL} winmm wldap32 ws2_32 ${TIFF} ${ZLIB})
else()
	target_link_libraries(makeahmap  ${CURL}  ${TIFF} ${ZLIB})
endif()
