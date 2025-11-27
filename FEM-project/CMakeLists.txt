cmake_minimum_required(VERSION 3.16)

project(FEM-project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(/bigobj)

add_executable(FEM-project main.cpp
    test_ttl.h
    test_eigen.h
    create_OBJ.h
    FEMobject.h FEMobject.cpp
    create_OBJ.cpp)

find_package(Eigen3 CONFIG REQUIRED)

include_directories("C:/Users/hogne/Programmering/DTE-3612_FEM_theory_and_implementation/TTL-master/include")

target_link_libraries(FEM-project PRIVATE
    "C:/Users/hogne/Programmering/DTE-3612_FEM_theory_and_implementation/BUILD/TTL_DEBUG/ttl.lib"
    Eigen3::Eigen
)
