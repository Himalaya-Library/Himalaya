cmake_minimum_required(VERSION 2.6)
project(h3cc)

add_executable(h3cc main.cpp)

install(TARGETS h3cc RUNTIME DESTINATION bin)
