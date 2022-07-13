mkdir build
cd build
cmake  -DCMAKE_PREFIX_PATH="$HOME/LCC/qmd-progress/install/;$HOME/LCC/bml/install" ../src/
make 
