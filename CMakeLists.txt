#
#	Created by Wei Wang on 2024-1-14.
#	Copyright © W. Wang. All rights reserved.
#

project(numerical_recipes_cpp)

# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.0)

# Cmake C/C++ compiler will use C++ 11, which equals to "set(CMAKE_CXX_STANDARD 11)".	
add_compile_options(-std=c++11)

# Add some new flags into original CXXFLAGS.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

# Add directory "src" into the project, CMake will also execute the CMakeLists.txt file in "src",
# so you can use the files and libraries in "src".
add_subdirectory(src)