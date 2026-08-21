# =============================================================================
# 2D ideal MHD on AMReX — CMake build
#
# Поддерживаемые платформы: Ubuntu (gcc/clang), macOS (Apple clang + Homebrew).
# Опции:
#   -DMHD_MPI=ON|OFF      (default ON)  — распараллеливание MPI (AMReX)
#   -DMHD_OPENMP=ON|OFF   (default ON)  — OpenMP-тайлинг внутри ранга
#   -DMHD_HDF5=ON|OFF     (default OFF) — вывод HDF5 (нужен libhdf5; ParaView
#                                          читает и нативные plotfile'ы AMReX)
#   -DMHD_FETCH_AMREX=ON|OFF (default ON) — скачать и собрать AMReX автоматически;
#                                          при OFF используется find_package(AMReX)
#
# Примеры:
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
#   cmake --build build -j 8
#   ./build/mhd2d inputs/orszag_tang.json
#
# macOS: OpenMP у Apple clang требует libomp:  brew install libomp open-mpi hdf5
# =============================================================================
cmake_minimum_required(VERSION 3.20)
project(mhd2d CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(MHD_MPI         "Enable MPI"                    ON)
option(MHD_OPENMP      "Enable OpenMP"                 ON)
option(MHD_HDF5        "Enable HDF5 plotfiles"         OFF)
option(MHD_FETCH_AMREX "Fetch and build AMReX via FetchContent" ON)

include(FetchContent)

# --- nlohmann/json (single-header, только заголовки) -------------------------
FetchContent_Declare(nlohmann_json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
)
FetchContent_MakeAvailable(nlohmann_json)

# --- AMReX --------------------------------------------------------------------
# Документация по сборке: https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html
if(MHD_FETCH_AMREX)
  set(AMReX_SPACEDIM 2          CACHE STRING "" FORCE)   # 2D-решатель
  set(AMReX_MPI      ${MHD_MPI} CACHE BOOL   "" FORCE)
  set(AMReX_OMP      ${MHD_OPENMP} CACHE BOOL "" FORCE)
  set(AMReX_HDF5     ${MHD_HDF5} CACHE BOOL  "" FORCE)
  set(AMReX_FORTRAN  OFF        CACHE BOOL   "" FORCE)
  set(AMReX_PIC      ON         CACHE BOOL   "" FORCE)
  set(AMReX_BUILD_TUTORIALS OFF CACHE BOOL   "" FORCE)
  FetchContent_Declare(amrex
    GIT_REPOSITORY https://github.com/AMReX-Codes/amrex.git
    GIT_TAG        25.01           # фиксированная версия — проверенный API
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(amrex)
  set(AMREX_TARGET AMReX::amrex_2d)
  if(NOT TARGET AMReX::amrex_2d)
    set(AMREX_TARGET AMReX::amrex)   # старые версии экспортируют один таргет
  endif()
else()
  find_package(AMReX REQUIRED COMPONENTS 2D)
  set(AMREX_TARGET AMReX::amrex_2d)
endif()

# --- исполняемый файл ----------------------------------------------------------
add_executable(mhd2d
  src/main.cpp
  src/MhdAmr.cpp
  src/Config.cpp
)
target_include_directories(mhd2d PRIVATE src)
target_compile_definitions(mhd2d PRIVATE MHD_WITH_AMREX)
target_link_libraries(mhd2d PRIVATE ${AMREX_TARGET} nlohmann_json::nlohmann_json)

if(MHD_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(mhd2d PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

# --- автономный верификационный драйвер (без AMReX) -----------------------------
# Использует те же ядра (src/kernels); удобен для быстрых проверок схемы.
add_executable(mhd2d_verify tests/standalone_verify.cpp)
target_include_directories(mhd2d_verify PRIVATE src/kernels)

# --- предупреждения --------------------------------------------------------------
foreach(tgt mhd2d mhd2d_verify)
  target_compile_options(${tgt} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra>)
endforeach()
