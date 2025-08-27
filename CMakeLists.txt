
cmake_minimum_required(VERSION 3.20)
project(WeigthedDelaunay3D LANGUAGES CXX)

option(WD3D_ENABLE_CUDA "Build with CUDA support" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (WD3D_ENABLE_CUDA)
  enable_language(CUDA)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(wd3d_common
    src/io/xyzw_io.cpp
    src/cpu/regular_cpu.cpp
    src/cpu/regular_cpu_bw.cpp
)

target_include_directories(wd3d_common PUBLIC include)
target_compile_definitions(wd3d_common PUBLIC $<$<NOT:$<BOOL:${WD3D_ENABLE_CUDA}>>:WD3D_NO_CUDA>)

if (WD3D_ENABLE_CUDA)
  add_library(wd3d_gpu
    src/gpu/power_predicates.cuh
    src/gpu/triangulation_kernels.cu
  )
  set_target_properties(wd3d_gpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  target_include_directories(wd3d_gpu PUBLIC include)
  target_link_libraries(wd3d_gpu PUBLIC wd3d_common)
endif()

add_executable(WeigthedDelaunay3D src/apps/WeigthedDelaunay3D.cpp)
target_link_libraries(WeigthedDelaunay3D PRIVATE wd3d_common $<$<BOOL:${WD3D_ENABLE_CUDA}>:wd3d_gpu>)

add_executable(EdgesWeigthedDelaunay3D src/apps/EdgesWeigthedDelaunay3D.cpp)
target_link_libraries(EdgesWeigthedDelaunay3D PRIVATE wd3d_common $<$<BOOL:${WD3D_ENABLE_CUDA}>:wd3d_gpu>)

enable_testing()
add_executable(test_tiny tests/test_tiny.cpp)
target_link_libraries(test_tiny PRIVATE wd3d_common)
add_test(NAME tiny COMMAND test_tiny ${CMAKE_CURRENT_SOURCE_DIR}/sample_data/tetra5.xyzw)
