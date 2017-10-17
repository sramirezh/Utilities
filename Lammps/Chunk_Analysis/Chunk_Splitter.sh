# This script reads the averaged-chunk data and creates individual files for every timestep
# The name of the file to be splitted is given as input in the command line. 

rm outfile* 2>/dev/null
rm *.chunk 2>/dev/null
rm Times.dat 2>/dev/null
name=$1 

pattern=`head -4 $name| tail -1| awk {'print $2'}` #This lines gets the pattern that has to be used to split the file
csplit --digits=4 -z --quiet --prefix=outfile $name "/[0-9] $pattern [0-9]/" "{*}" #Splits between things with the

mv outfile0000 header

#Now change all the filenames and create the file Times.dat
for f in outfile*
do
var=`head -1 $f|awk {'print $1'}` #gets the time
echo $var>>Times.dat
mv $f $var.chunk
done
