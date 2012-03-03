README FILE


This source code is based on the pure-Fortan QTCM1 version 2.3
(August 2002) from http://www.atmos.ucla.edu/~csi/.

To create the shared object files, do the following:

(1) Copy the correct makefile for your platform from ./Makefiles
to this directory, and rename it makefile.

(2) At the command line, type:
make clean && make _qtcm_full_365.so && make _qtcm_parts_365.so

The _qtcm_full_365.so and _qtcm_parts_365.so files will be put in
../lib automatically.
