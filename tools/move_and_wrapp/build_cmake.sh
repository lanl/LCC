mkdir build
cd build
cmake  -DCMAKE_Fortran_COMPILER="gfortran" -DEXTRA_FCFAGS="-I$HOME/LCC/bml/install/include" -DCMAKE_PREFIX_PATH="$HOME/LCC/qmd-progress/install/;$HOME/LCC/bml/install;$HOME/LCC/tools/move_and_wrapp/" ../
make 
