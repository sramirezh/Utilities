# This script reads the trajectory data and creates individual files for every timestep
# The name of the file to be splitted is given as input in the command line. 
#Notice that the number of particles inside the file for every step can vary.

rm *.cxyz
rm Times.dat
name=$1 


csplit --digits=4 -z --quiet --prefix=outfile $name "/^[0-9]*$/" "{*}" #Splits based on lines that contain at the beggining a number and after that number, the line finishes

#Now change all the filenames and create the file Times.dat
for f in outfile*
do
var=`head -2 $f|tail -1|awk {'print $3'}` #gets the time
echo $var>>Times.dat
mv $f $var.cxyz
done
