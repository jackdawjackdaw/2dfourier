cmake_minimum_required(VERSION 2.6)

add_library(decompf2d SHARED fourier-decomp)
target_link_libraries(decompf2d ${libs})

install_targets( /lib decompf2d)
install_files( /include decompf2d.h)

add_executable(f2d-driver-example example-driver)
target_link_libraries(f2d-driver-example decompf2d ${libs})

install_targets( /bin f2d-driver-example)
