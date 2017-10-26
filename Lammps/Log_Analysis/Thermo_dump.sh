#This script takes the log.lammps and generates a file just containing the parameters dumped by the thermo_style
#Creates a file called parameters.dat
#Just to check
File=$1
rm Parameters.dat 2>/dev/null
a=$(grep -n "Per MPI" $File| awk -F":" '{print $1}') #the awk gets the first column when they are separated by :
tail -n +$((a+1)) $File > out
b=$(grep -n "Loop time" out| awk -F":" '{print $1}')
if [ -z "$b" ];
then
    mv out  Parameters.dat
else
    head -$((b-1)) out> Parameters.dat
fi
rm out 2>/dev/null
