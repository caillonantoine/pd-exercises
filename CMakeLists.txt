cmake_minimum_required(VERSION 3.0)

project(pd-exercises)

set(SRC duck_tilde.c)

add_library(duck SHARED ${SRC})

set_target_properties(duck PROPERTIES PREFIX "" SUFFIX "~.pd_linux")


#target_link_libraries(duck pd)
