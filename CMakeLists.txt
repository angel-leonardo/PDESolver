cmake_minimum_required(VERSION 4.0)
project(Week1)

set(CMAKE_CXX_STANDARD 20)

add_executable(Week1 main.cpp
        sparse_matrix.cpp
        sparse_matrix.h
        PoissonEQ.cpp
        PoissonEQ.h
        polygon_boundary.cpp
        polygon_boundary.h
        ParabolaPDE.cpp
        ParabolaPDE.h
        HyperbolicEQ.cpp
        HyperbolicEQ.h)
