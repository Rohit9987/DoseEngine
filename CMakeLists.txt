cmake_minimum_required(VERSION 3.20)
project(dose_engine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

FetchContent_Declare(
    head_model
    GIT_REPOSITORY https://github.com/Rohit9987/head_model.git
    GIT_TAG master
)

find_package(OpenMP REQUIRED)

FetchContent_MakeAvailable(head_model)

add_executable(dose_engine
	src/main.cc 
	src/TermaCalculator.cpp 
	src/KernelLoader.cpp 
	src/DoseConvolver.cpp
	src/KernelStencil.cpp
	src/StencilConvolver.cpp
	src/SpectrumKernelBuilder.cpp
)

target_include_directories(dose_engine PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(dose_engine
	PRIVATE head_model OpenMP::OpenMP_CXX
)
