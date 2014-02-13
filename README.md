2dfourier {#mainpage}
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

Some automatic documentation can be generated by running

    doxygen Doxyfile 

A shared library `libdecompf2d` is built and installed along with headers `decompf2d.h`. Consult the file
`example-driver.c` to see how to use the library fns.

## Running the example

`example_driver.c` shows how to read a file and compute the coefficients Amn. A file in the input format
required is provide in `example/example-event.dat`, this event is based on a 200 x 200 pt grid. 

To run the example analysis and compute up to the 8th order moments do, supposing you have installed the driver somewhere in your path

    f2d-driver-example ./example/example-event.dat 200 8 8

You should see an output like this:

     # reading from: ./example/example-event.dat
     # input grid size: (200 x 200)
     # moments computed up to: 8 8
     # com: -0.199229 -0.002554 (79 99)
     (-8 0)[1.326949 -0.100377] (-8 1)[-1.714913 0.082701] ...
     ...
     L2: 44.168916
     M1: 376.742586
     H1: 480.405883
     Rsq: 117.299565

The code prints out a grid of values: `(M, N) [ Re[A_MN], Im[A_MN]]`, and then lists the values of the various
norms for the event. The example event is one element of the ensemble of N_a=1 URQMD events discussed in the
article [1], after consulting fig.4 and fig.5 we se that the printed out values are fairly representative


[1]: http://arxiv.org/abs/arXiv:1204.5774
[2]: http://www.cmake.org/


