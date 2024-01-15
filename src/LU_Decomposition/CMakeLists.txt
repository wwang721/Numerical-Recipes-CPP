#
#	Created by Wei Wang on 2024-1-14.
#	Copyright © W. Wang. All rights reserved.
#

project(lu_decomposition)


# Add all CPP files to the executable
file(GLOB SOURCES "*.cpp")

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)
add_executable(ludcmp.exe ${SOURCES})
