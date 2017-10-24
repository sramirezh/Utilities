#!/usr/bin/env bash
#############################################
# This code is intended to analyze everything from the Measurement run
#############################################
Lets seee
dir=$(dirname $0) #to get the directory where the script and other source files are.

printf "\n##########################################################################\n"
echo "Analyzing the log file"
echo "##########################################################################"
bash $dir/Log_Analysis/Thermo_dump.sh log.lammps
python $dir/Log_Analysis/Thermo_Analyser.py

printf "\nGenerated Parameters.dat \n"

printf "\n##########################################################################\n"
echo "Analyzing the trajectory File"
echo "##########################################################################"
bash $dir/Trajectory_Analysis/1Trajectory_Splitter.sh trajectory.xyz
python $dir/Trajectory_Analysis/1Trajectory_Analizer.py
printf "\nGenerated 1Trajectory.xyz and Zshift.dat \n"

printf "\n##########################################################################\n"
echo "Analyzing the Chunk properties"
echo "##########################################################################"


printf "\n**************************************************************************\n"
echo "Averaging the Solute properties"
echo "**************************************************************************"
bash $dir/Chunk_Analysis/Chunk_Splitter.sh Sproperties.all
python $dir/Chunk_Analysis/Chunk_Analyser.py
mv Averages.dat SAverages.dat

printf "\nGenerated SAverages.dat  \n"

printf "\n**************************************************************************\n"
echo "Averaging the Solvent properties"
echo "**************************************************************************"
bash $dir/Chunk_Analysis/Chunk_Splitter.sh Lproperties.all
python $dir/Chunk_Analysis/Chunk_Analyser.py
mv Averages.dat LAverages.dat

printf "\nGenerated LAverages.dat  \n"

printf "\n**************************************************************************\n"
echo "Averaging the Fluid properties"
echo "**************************************************************************"
bash $dir/Chunk_Analysis/Chunk_Splitter.sh properties.all
python $dir/Chunk_Analysis/Chunk_Analyser.py
mv Averages.dat AAverages.dat


printf "\nGenerated AAverages.dat \n"

rm *.chunk*

printf "\n##########################################################################\n"
echo "Analyzing the properties"
echo "##########################################################################"

python $dir/Property_Analysis/Property_Analysis.py
