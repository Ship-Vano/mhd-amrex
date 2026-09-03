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

include(CTest)

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

# Независимый эталон Брио–Ву (фаза T06): центральная схема Куртганова–Тадмора.
# Намеренно НЕ подключает src/kernels — эталон не должен разделять код со
# схемой, которую проверяет.
add_executable(mhd2d_briowu_reference tests/briowu_reference.cpp)

add_executable(mhd2d_kernel_unit tests/kernel_unit.cpp)
target_include_directories(mhd2d_kernel_unit PRIVATE src/kernels)

add_executable(mhd2d_kernel_numerics tests/kernel_numerics.cpp)
target_include_directories(mhd2d_kernel_numerics PRIVATE src/kernels)

# Разбор конфигурации изолирован от AMReX, поэтому schema errors проверяются
# быстро и до выделения mesh/field data.
add_executable(mhd2d_config_validate
  tests/config_validate.cpp
  src/Config.cpp
)
target_include_directories(mhd2d_config_validate PRIVATE src)
target_link_libraries(mhd2d_config_validate PRIVATE nlohmann_json::nlohmann_json)

# --- предупреждения --------------------------------------------------------------
foreach(tgt mhd2d mhd2d_verify mhd2d_kernel_unit mhd2d_kernel_numerics mhd2d_config_validate mhd2d_briowu_reference)
  target_compile_options(${tgt} PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra>)
endforeach()

if(BUILD_TESTING)
  find_package(Python3 REQUIRED COMPONENTS Interpreter)

  add_test(NAME kernel.unit COMMAND mhd2d_kernel_unit)
  add_test(NAME kernel.numerics COMMAND mhd2d_kernel_numerics)
  add_test(NAME standalone.briowu
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_verify.py
                   --executable $<TARGET_FILE:mhd2d_verify> --case briowu)
  add_test(NAME standalone.alfven32
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_verify.py
                   --executable $<TARGET_FILE:mhd2d_verify> --case alfven --resolution 32)
  add_test(NAME standalone.alfven_order
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_alfven_order.py
                   --executable $<TARGET_FILE:mhd2d_verify> --coarse 64 --fine 128 --min-order 1.8)
  add_test(NAME standalone.dai_woodward
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_verify.py
                   --executable $<TARGET_FILE:mhd2d_verify> --case dw1d --resolution 200)
  add_test(NAME standalone.loop
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_verify.py
                   --executable $<TARGET_FILE:mhd2d_verify> --case loop --resolution 64)
  add_test(NAME config.valid
           COMMAND mhd2d_config_validate ${CMAKE_CURRENT_SOURCE_DIR}/tests/config/valid_minimal.json)
  add_test(NAME config.rejects_unknown_key
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/assert_config_rejected.py
                   --executable $<TARGET_FILE:mhd2d_config_validate>
                   --config ${CMAKE_CURRENT_SOURCE_DIR}/tests/config/invalid_unknown_key.json)
  add_test(NAME solver.compute_only
           COMMAND mhd2d ${CMAKE_CURRENT_SOURCE_DIR}/tests/config/valid_minimal.json)
  # NEW-004: режимы аблации задаются конфигом и действительно влияют на результат.
  add_test(NAME config.ablation_modes
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_ablation_modes.py
                   --executable $<TARGET_FILE:mhd2d>)

  # Гейт T06: две независимые схемы (центральная KT без решателя Римана и
  # HLLD+MUSCL) должны сходиться к одному слабому решению, включая составную
  # волну. Тест падает и если они разошлись, и если эталон не сошёлся сам.
  add_test(NAME briowu.independent_reference
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_briowu_reference.py
                   --reference $<TARGET_FILE:mhd2d_briowu_reference>
                   --verify $<TARGET_FILE:mhd2d_verify>
                   --fronts ${CMAKE_CURRENT_SOURCE_DIR}/scripts/briowu_fronts.py)

  # ADR 0001 (гейт TDG0): слой ядер не должен структурно зависеть от AMReX.
  add_test(NAME arch.kernel_purity
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_kernel_purity.py
                   --kernel-dir ${CMAKE_CURRENT_SOURCE_DIR}/src/kernels)

  # Обязательная матрица канонических задач (PROGRAM_PLAN). Диапазоны rho и p --
  # литературные цветовые шкалы, поэтому тест ловит не только падение, но и
  # смещение фронта или потерю признака.
  add_test(NAME canonical.constant_state
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_solver_case.py
                   --executable $<TARGET_FILE:mhd2d>
                   --config ${CMAKE_CURRENT_SOURCE_DIR}/inputs/uniform_const.json
                   --exact-constant)
  add_test(NAME canonical.orszag_tang
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_solver_case.py
                   --executable $<TARGET_FILE:mhd2d>
                   --config ${CMAKE_CURRENT_SOURCE_DIR}/inputs/orszag_tang_uniform.json
                   --rho-min 0.0913 --rho-max 0.4915
                   --p-min 0.0281 --p-max 0.5045 --tol 0.08)
  add_test(NAME canonical.rotor
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_solver_case.py
                   --executable $<TARGET_FILE:mhd2d>
                   --config ${CMAKE_CURRENT_SOURCE_DIR}/inputs/rotor_uniform.json
                   --rho-min 0.5885 --rho-max 11.55
                   --p-min 0.0354 --p-max 1.855 --tol 0.08)
  # Взрыв в среде с beta = 0.2: здесь откат HLLD->HLL обязан срабатывать, но
  # состояние должно остаться положительным, а div B -- на уровне округления.
  add_test(NAME canonical.mhd_blast
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_solver_case.py
                   --executable $<TARGET_FILE:mhd2d>
                   --config ${CMAKE_CURRENT_SOURCE_DIR}/inputs/mhd_blast.json
                   --rho-min 0.1526 --rho-max 2.506
                   --p-min 0.0785 --p-max 0.9409 --tol 0.08
                   --max-fallbacks 200 --max-floors 500)

  # Гейт T07: согласование газовых потоков на стыке уровней. Тест сам проверяет,
  # что стык уровней вообще существует и что при reflux=false дефект появляется.
  # Перестроение сетки обязано не двигать интегралы: cell_cons_interp
  # консервативен по построению. Расхождение тепловой и магнитной частей при
  # этом ожидаемо (ошибка интерполяции нелинейных величин) и здесь не гейтится.
  add_test(NAME amr.regrid_conservation
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_amr_conservation.py
                   --executable $<TARGET_FILE:mhd2d>
                   --config ${CMAKE_CURRENT_SOURCE_DIR}/inputs/amr_regrid_conservation.json
                   --tol 1e-6 --max-regrid-jump 1e-12 --skip-defect-check)
  add_test(NAME amr.conservation
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_amr_conservation.py
                   --executable $<TARGET_FILE:mhd2d>
                   --config ${CMAKE_CURRENT_SOURCE_DIR}/inputs/amr_conservation.json)
  add_test(NAME manifest.repeatable
           COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_manifest_repeat.py
                   --runner ${CMAKE_CURRENT_SOURCE_DIR}/scripts/write_run_manifest.py
                   --config ${CMAKE_CURRENT_SOURCE_DIR}/tests/config/valid_minimal.json
                   --out-dir ${CMAKE_CURRENT_BINARY_DIR}/test-manifests)
endif()
