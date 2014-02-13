2dfourier
=========
C.Coleman-Smith (cec24@phy.duke.edu)

Routines for carrying out a 2d Fourier decomposition for Heavy Ion events  and computing associated measures of event roughness. 

### Publication details
"Classification of initial state granularity via 2d Fourier Expansion" [1]

Christopher E. Coleman-Smith (Duke U.), Hannah Petersen (Duke U. & Frankfurt U., FIAS), Robert L. Wolpert (Duke U.). Apr 2012. 9 pp.

Published in J.Phys. G40 (2013) 095103 

# Getting Started

## Requires:

* CMAKE [2]
* A C compiler
* libGSL

## Building & Installing

This project uses CMAKE to generate Makefiles, it is canonical to do out of place builds using cmake. An "out of place" builds puts all the temporary files and compiler junk into a directory that is outside the source tree.

From the project root do:

    mkdir ./build
    cd ./build
    cmake ..
    make && make install

Cmake defaults to installing things in /usr/local, if you don't want that you should set invoke cmake as

    cmake -DCMAKE_INSTALL_PREFIX:PATH=/your/install/path ..


    


    



[1]: http://arxiv.org/abs/arXiv:1204.5774
[2]: http://www.cmake.org/


