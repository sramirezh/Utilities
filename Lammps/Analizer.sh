#############################################
# This code is intended to analyze everything from the Measurement run 
#############################################

dir=$(dirname $0) #to get the directory where the script and other source files are.

printf "\n##########################################################################\n"
echo "Analizing the log file"
echo "##########################################################################"
bash $dir/Thermo_dump.sh log.lammps
python $dir/Thermo_analyzer.py

printf "\nGenerated Parameters.dat \n"

printf "\n##########################################################################\n"
echo "Analizing the trajectory File"
echo "##########################################################################"
bash $dir/Trajectory_Analysis/1Trajectory_Splitter.sh trajectory.xyz
python $dir/Trajectory_Analysis/1Trajectory_Analizer.py
printf "\nGenerated 1Trajectory.xyz and Zshift.dat \n"

printf "\n##########################################################################\n"
echo "Analizing the Chunk properties"
echo "##########################################################################"


printf "\n**************************************************************************\n"
echo "Analizing the Solute properties"
echo "**************************************************************************"
bash $dir/Chunk_Splitter.sh Sproperties.all
python $dir/Chunk_Analyzer.py
mv Averages.dat SAverages.dat
for f in *.chunk;do mv "$f" "$f"s;done  #To rename as .chunks to analyze afterwards with Force_Factor.py

printf "\nGenerated SAverages.dat  \n"

printf "\n**************************************************************************\n"
echo "Analizing the Solvent properties"
echo "**************************************************************************"
bash $dir/Chunk_Splitter.sh Lproperties.all
python $dir/Chunk_Analyzer.py
mv Averages.dat LAverages.dat

printf "\nGenerated LAverages.dat  \n"

printf "\n**************************************************************************\n"
echo "Analizing the Fluid properties"
echo "**************************************************************************"
bash $dir/Chunk_Splitter.sh properties.all
python $dir/Chunk_Analyzer.py
mv Averages.dat AAverages.dat


printf "\nGenerated AAverages.dat \n"


printf "\n##########################################################################\n"
echo "Analizing the force factor for the chemical potential"
echo "##########################################################################"

python $dir/Force_Factor.py


printf "\n##########################################################################\n"
echo "Analizing the Concentration Distribution"
echo "##########################################################################"
python $dir/Concen_dist.py

rm *.chunk* 



