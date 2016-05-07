cmake_minimum_required(VERSION 2.8)

# MinGW Specific
set(CMAKE_CXX_FLAGS "-std=c++0x -Wall -m32 -DDLL_EXPORT -DCURL_STATICLIB")
set(CMAKE_RC_FLAGS " -F pe-i386")

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

add_executable(makeahmap
	makeahmap.cpp
	downloader.cpp
	makeahmap.rc
	)

target_link_libraries(makeahmap  ${CURL} winmm wldap32 ws2_32 ${TIFF} ${ZLIB})
