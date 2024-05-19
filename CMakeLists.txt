cmake_minimum_required(VERSION 3.10)

project(PhyCollision)

# 设置头文件路径
include_directories(${PROJECT_SOURCE_DIR}/include)

# 设置源文件路径
file(GLOB SOURCES "src/*.cpp")

find_package(Eigen3 3.4.0 REQUIRED NO_MODULE)

# win 电脑需要自己配置这一行
set(MY_EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/SDK/Eigen3/include/eigen3")

include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(${MY_EIGEN3_INCLUDE_DIR})

# 输出EIGEN3_INCLUDE_DIR
message(STATUS "EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "EIGEN3_INCLUDE_DIR: ${MY_EIGEN3_INCLUDE_DIR}")

# 添加一个库，将源文件编译为动态链接库
# add_library(PhyColiision SHARED ${SOURCES})

add_library(PhyColiision STATIC ${SOURCES})

# 设置库的版本
set_target_properties(PhyColiision PROPERTIES VERSION 1.0 SOVERSION 1)

# Create test executable
add_executable(testMain testMain.cpp)

# Link test executable with the shared library
target_link_libraries(testMain PhyColiision Eigen3::Eigen)