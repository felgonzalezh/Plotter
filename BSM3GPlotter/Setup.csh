#!/bin/bash

printf "\n"
echo Before running the plotting routine, we need to create the analysis directory. What do you want to name your directory?
printf "\n"

read dirname

#for y in `ls -d */`
#do

#  new_file=`echo $y | sed 's!/!!'`
#  cd ${new_file}
#  mkdir $dirname
#  cp default/* $dirname
#  cd ./..
#  cp BSM3GNormalizer ${new_file}/$dirname

#done

rm moveRootFilesToDirectories.csh
cp moveRootFilesToDirectories_default.csh moveRootFilesToDirectories.csh
sed -i -e s/TEMPDIRECTORY/"$dirname"/g moveRootFilesToDirectories.csh
./moveRootFilesToDirectories.csh

printf "\n"
echo The analysis root files have been moved to directory $dirname

rm main.in
cp main.in.default main.in
sed -i -e s/TEMPDIRECTORY/"$dirname"/g main.in

rm BSM3GPlotMaker.in
cp BSM3GPlotMaker.in.default BSM3GPlotMaker.in
sed -i -e s/TEMPDIRECTORY/"$dirname"/g BSM3GPlotMaker.in

rm BSM3GCutFlowEffTable.in
cp BSM3GCutFlowEffTable.in.default BSM3GCutFlowEffTable.in
sed -i -e s/TEMPDIRECTORY/"$dirname"/g BSM3GCutFlowEffTable.in

printf "\n"
echo The Plotter configuration files have been modified accordingly
printf "\n"

