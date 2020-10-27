# QGPU
Quality guide phase unwrapping based on linked lists


This repository contains the source code and data files of the phase unwrapping algorithms published in https://doi.org/10.1364/AO.57.003126. This is a separate QPGU algorithms implementations used for comparison againts the Pruning strategy. Here, it is provided two versions of the QGPU algorithm based on partitions of the quality values interval using linked lists. A traditional linked list was implemented using dynamic memory and a variant of a linked list using a pivot to speed up the insertion operation.   

See the LICENSE.txt file for the licensing statement and Disclaimer of Warranty.

The source code is written in C on the gcc using Linux Ubuntu operating system. The organization of this repository is the following:

1. I2L2: Implementation of the quality guide phase unwrapping algorithm using a linked list for the adjoin list. 
2. data: Contains the wrapped phase maps, quality maps and mask files used in the submitted paper.


Compiling:

Each implementation of the QGPU algorithms in the repository contains three directories, say include, lib and src. The header files in each implementation are in the include directory and the source files are in the src directory. To generate the executable program, locate in the corresponding src directory and type make in order to compile the project. To execute the program, type ./main- in the current src directory.

The project can be compiled either in MAC OS x and in Linux operating systems with minor modifications.
