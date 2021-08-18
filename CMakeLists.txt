cmake_minimum_required(VERSION 3.15)
project(stew C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fopenmp")

include_directories(include)
file(GLOB INCLUDES include/*.h)
file(GLOB SOURCES src/*.c)
add_executable(stew ${SOURCES} ${INCLUDES})