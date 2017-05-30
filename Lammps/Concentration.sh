#This scripsts splits both the Solute and total property files, it has two arguments:
# INPUT1=Solute Properties file.
# INPUT2=all property file.

File1=$1
File2=$2
echo "The file with the Solute properties is: $File1"
echo "The file with the fluid properties is: $File2"

echo "Splitting the files"
bash Chunk_Splitter.sh $File1

for f in *.chunk;do mv "$f" "$f"s;done
bash Chunk_Splitter.sh $File2
echo "Finished!!!"
