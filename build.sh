mkdir build
cd build
MY_PATH=`pwd`
cmake  -DCMAKE_INSTALL_BINDIR="$MY_PATH/bin" -DCMAKE_PREFIX_PATH="$MY_PATH/../qmd-progress/install/;$MY_PATH/..//bml/install" ../src/
make install 
