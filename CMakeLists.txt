cmake_minimum_required(VERSION 3.0)

project(pd_external)

add_library(duck SHARED duck_tilde.c)
set_target_properties(duck PROPERTIES PREFIX "" SUFFIX "~.pd_linux")

add_library(pan SHARED signal_in_out_base.c)
set_target_properties(pan PROPERTIES PREFIX "" SUFFIX "~.pd_linux")
