#!/bin/sh

IFS=$'\n'; 
for line in $(cat $1); do 
  line="${line#"${line%%[![:space:]]*}"}"
  line="${line%"${line##*[![:space:]]}"}"
  echo "";
  echo "--------------------------------------------------------------------------------------------------------";
  echo "directory $line";
  echo "--------------------------------------------------------------------------------------------------------";
  echo "";
  cd $line
  ./BSM3GNormalizer Normalizer.in
  cd ./../..
done

./BSM3GPlotMaker BSM3GPlotMaker.in
./BSM3GCutFlowEffTable BSM3GCutFlowEffTable.in
