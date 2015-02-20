cmake_minimum_required (VERSION 2.6)
project (gcmethod_1d)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(${PROJECT_NAME}_headers gcmethod_1d.h)
set(${PROJECT_NAME}_sources main.cpp)


add_executable (${PROJECT_NAME} ${${PROJECT_NAME}_headers} ${${PROJECT_NAME}_sources})
target_link_libraries (${PROJECT_NAME)

