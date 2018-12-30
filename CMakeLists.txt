cmake_minimum_required(VERSION 3.0)

project(pd-exercises)

add_library(duck SHARED duck_tilde.c)
set_target_properties(duck PROPERTIES PREFIX "" SUFFIX "~.pd_linux")

set(SRC_cross_synth cross_synth.c kiss_fft/kiss_fft.c)

add_library(cross_synth SHARED ${SRC_cross_synth})
set_target_properties(cross_synth PROPERTIES PREFIX "" SUFFIX "~.pd_linux")
