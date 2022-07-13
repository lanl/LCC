Encapsulating orthogonal vectors
================================

Code for finding encapsulating orthogonal vectors for a given lattice.
Enter a given direction (surface normal) to be given two orthogonal vectors.
Scalars will be generated that correspond to an error value.
The error value is related to the periodicity of the orthogonal vectors.


Components in directory:
-  `README.md` - details
-  `orth_calc.cxx`, `func_v.h` - source code
-  `run.sh` - example run script
-  `latticefile.dat` - read in lattice vector file from command lattice=read


Compile the code with:

   $ g++ -std=c++11 orth_calc.cxx


see `run.sh` for example on commandline inputs
Inputs:
- `nx`	- x value for direction of interest (numerical)
- `ny`	- y value for direction	of interest (numerical)
- `nz`	- z value for direction	of interest (numerical)
- `random`	- Selects random direction of interest (yes, no)
- `lattice` - lattice vector input (fcc, bcc, cubic, read) 
- `t` 	- lattice vector scalar (numerical)
- `tmax`	- Max scalar length for direction of interest (numerical > 0)
- `tstart`  - Starting scalar for direction of interest (numerical > 0)
- `tol`	- Tolerance for Error function (numerical)
- `tstep`	- Loops through tsart to tmax with this step (numberical > 0)

