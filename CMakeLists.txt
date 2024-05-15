cmake_minimum_required(VERSION 3.10)

project(PhyCollision)

# 设置头文件路径
include_directories(${PROJECT_SOURCE_DIR}/include)

# 设置源文件路径
file(GLOB SOURCES "src/*.cpp")

find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)
include_directories("/opt/homebrew/include")

# 添加一个库，将源文件编译为动态链接库
add_library(PhyColiision SHARED ${SOURCES})

# 设置库的版本
set_target_properties(PhyColiision PROPERTIES VERSION 1.0 SOVERSION 1)

# Create test executable
add_executable(testMain testMain.cpp)

# Link test executable with the shared library
target_link_libraries(testMain PhyColiision)