cmake_minimum_required(VERSION 3.18)
set(LLVM_MOS_PLATFORM sim)
find_package(llvm-mos-sdk REQUIRED)
project(atto-pi)

add_executable(atto-pi.prg atto-pi.c)

# Use recommended optimization settings from AI.md
set_target_properties(atto-pi.prg PROPERTIES
    CMAKE_BUILD_TYPE MinSizeRel
)
target_compile_options(atto-pi.prg PRIVATE -Os -flto -g)
target_link_options(atto-pi.prg PRIVATE -flto)
