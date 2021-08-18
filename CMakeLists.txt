cmake_minimum_required(VERSION 3.15)
project(stew C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fopenmp")

include_directories(include)
file(GLOB SOURCES src/*.c)
add_executable(stew ${SOURCES} include/log.h src/log.c include/ketopt.h include/kseq.h include/hll.h include/hll_private.h src/hll.c src/city.c include/city.h include/ascii.h)