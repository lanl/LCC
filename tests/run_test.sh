#!/bin/bash

# Script to test the LCC program.

MY_PATH=`pwd`               

RUN="../build/lcc_main" 

echo -e "\nTesting LCC ...\n"

#set -e     

for name in sc_bulk fcc_bulk vectors_bulk fcc_bravais_growth \
  fcc_planes; do

  INFILE=$name".in"
  REF="coords_"$name".pdb"
  COORDS=coords.pdb

  cp  ./tests_data/$INFILE test.in
  cp  ./tests_data/$REF .
  
  $RUN test.in > out
  diff $REF $COORDS > tmp 
  FAILED=`wc -l tmp | awk 'NF>1{print $1}'`

  if [ $FAILED != "0" ]; then
        echo  "Test for" $name "... FAILED"
        exit 1;
  else
        echo "Test for" $name "... PASSED"
  fi

  rm $REF test.in tmp

done

