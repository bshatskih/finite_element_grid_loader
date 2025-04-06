cmake_minimum_required(VERSION 3.26)
project(finite_element_grid_loader)

set(CMAKE_CXX_STANDARD 23)

add_executable(finite_element_grid_loader main.cpp
        MeshLoader.h
        MeshLoader.cpp
        AneuMeshLoader.h
        AneuMeshLoader.cpp
        Node.h
        FiniteElement.h
        PlainClass.h)
