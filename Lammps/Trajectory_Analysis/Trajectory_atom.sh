#!/usr/bin/env sh
# This script reads the trajectory data and creates individual files for every timestep
# The name of the file to be splitted is given as input in the command line.
#Notice that the number of particles inside the file for every step can vary.

###############################################################################
#Function definition
###############################################################################
TotalSplit(){
csplit --digits=4 -z --quiet --prefix=outfile $name "/^ITEM: TIMESTEP/" "{*}" #Splits based on lines that contain at the beggining a number and after that number, the line finishes
#Now change all the filenames and create the file Times.dat
for f in outfile*
do
var=`sed -n '2p' $f` #gets the time
echo $var>>Times.dat
tail -n +10 $f>tmp
cut -f1 -d" " --complement tmp > $f
mv $f $var.cxyz
done
}

PartialSplit(){
  
  NAtoms=`sed -n '4p' $name` #Number of atoms
  Nlines=$(($NAtoms+9))
  head -$Nlines $name > 0.cxyz
  echo "0">Times.dat


}
usage() { echo "Usage: $0 [-i <inputFile>] [-s <Initial or Total>]" 1>&2; exit 1; }

###############################################################################
while getopts ":i:s:" o; do
    case "${o}" in
        i)
            name=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

if [ -z "${name}" ]  ; then  #If there is no value in the arguments
    usage
fi

rm *.cxyz 2>/dev/null
rm Times.dat 2>/dev/null

if [ "$s" == "Initial" ]; then
echo "Extracting the first configuration from the trajectory file..."
PartialSplit
fi

if [ -z "${s}" ]  ; then  #If there is no value in the arguments
    echo "Splitting the entire file ..."
    TotalSplit
fi

rm tmp 2>/dev/null 
