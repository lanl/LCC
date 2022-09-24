Using these python scripts: Quick start
=======================================

Get orthogonal Miller indices:
------------------------------

  ./getOtrthoMillers.py -hkl 1 0 0

Will get orthogonal Miller indices to (100)

Get orthogonal pseudo congruent vectors:
----------------------------------------

  ./getOtrthoCongruent.py -n0 1 1 0 -lattice "latticeFile.dat"  -tol 0.0001

Will get two pseudo congruent orthogonal vectors to (1,1,0), as 
well as an elongation factor "l" so that l x (1,1,0) is a pseudo 
congruent vector within a tolerace `tol`.

