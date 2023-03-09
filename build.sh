mkdir build
cd build
MY_PATH=`pwd`
cmake -DEXTRA_LINK_FLAGS="-llapack -lblas" -DCMAKE_Fortran_FLAGS="-g -fopenmp" -DCMAKE_INSTALL_BINDIR="$MY_PATH/bin" -DCMAKE_PREFIX_PATH="$MY_PATH/../qmd-progress/install/;$MY_PATH/../bml/install" ../src/
make install 
