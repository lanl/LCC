mkdir build
cd build
MY_PATH=`pwd`
cmake  -DCMAKE_INSTALL_BINDIR="$MY_PATH" -DCMAKE_PREFIX_PATH="$HOME/LCC/qmd-progress/install/;$HOME/LCC/bml/install" ../src/
make install 
