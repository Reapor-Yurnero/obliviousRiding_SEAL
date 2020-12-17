# CMakeLists.txt
cmake_minimum_required(VERSION 3.0)
project(FHE_app)

find_package(SEAL 3.6 REQUIRED)

add_executable(main main.cpp)

target_link_libraries(main PRIVATE SEAL::seal)