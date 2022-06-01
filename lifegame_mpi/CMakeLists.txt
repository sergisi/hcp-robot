cmake_minimum_required(VERSION 3.23)
project(lifegame_mpi C)

set(CMAKE_C_STANDARD 23)

add_executable(lifegame_mpi main.c debug_tools.h)
