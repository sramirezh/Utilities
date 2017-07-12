#This scripts runs the chunk splitter and analyzer for all the three components of the system.
#It is required to change first the zshift on the chunk_analyzer.py

Dir=$(dirname $0) #to get the directory where the script and other source files are.

echo"##################################################################"
echo "Analizing the Solute properties"
echo"##################################################################"
bash $Dir/Chunk_Splitter.sh Sproperties.all
python $Dir/Chunk_Analyzer.py
mv Averages.dat SAverages.dat
echo"##################################################################"
echo "Analizing the Solvent properties"
echo"##################################################################"
bash $Dir/Chunk_Splitter.sh Lproperties.all
python $Dir/Chunk_Analyzer.py
mv Averages.dat LAverages.dat
echo"##################################################################"
echo "Analizing the Fluid properties"
echo"##################################################################"
bash $Dir/Chunk_Splitter.sh properties.all
python $Dir/Chunk_Analyzer.py
mv Averages.dat AAverages.dat

rm *.chunk
