#!/bin/bash
g++ -std=c++11 -o ovcalc orth_calc.cxx
./ovcalc nx=1 ny=0 nz=0 lattice=read tmax=100 tstart=1 tstep=0.005 tol=0.000001
