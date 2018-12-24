project(duck)

add_definitions(-DMSW)

set(SRC duck_tilde.c)

include_directories("C:/Program Files (x86)/pd/include")
link_directories("C:/Program Files (x86)/pd/bin")

add_library(duck SHARED ${SRC} duck.def)

target_link_libraries(duck pd)
