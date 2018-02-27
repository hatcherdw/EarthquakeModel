# EarthquakeModel
Final project for PHY 517. 


This directory conatains the source code, Makefile, and data files for my earthquake model.

constants.mod is a file created from a module within the source code. It should NOT be edited. Changes can be made to the module within the source code.

Quake.f90 is the Fortran source code.

data.xyz is a data file with time steps, block numbers, and block velocities. The .xyz extension is easy to move out of the web directory.
To write to this file the output of Quake should be redirected using '>', otherwise output is placed on the screen.


Quake can be compiled with 'make' or by running run.sh (this also runs Quake).
