cmake_minimum_required(VERSION 3.28)
project(optimized_ntt CXX)

set(CMAKE_C_STANDARD 11)

include_directories(ref/dilithium2/ntt-opt)


add_executable(optimized_ntt ntt.hpp ntt.cpp)