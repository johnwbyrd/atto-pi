cmake_minimum_required(VERSION 3.18)
set(LLVM_MOS_PLATFORM sim)
find_package(llvm-mos-sdk REQUIRED)
project(pi)

add_executable(pi.prg pi.c)

# Use recommended optimization settings from AI.md
set_target_properties(pi.prg PROPERTIES
    CMAKE_BUILD_TYPE MinSizeRel
)
target_compile_options(pi.prg PRIVATE -Os -flto -g)
target_link_options(pi.prg PRIVATE -flto)
