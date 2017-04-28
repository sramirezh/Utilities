#This script takes the log.lammps and generates a file just containing the parameters dumped by the thermo_style
#Creates a file called parameters.dat
rm parameters.dat
a=$(grep -n "Per MPI" log.lammps| awk -F":" '{print $1}') #the awk gets the first column when they are separated by :
tail -n +$((a+1)) log.lammps > out
b=$(grep -n "Loop time" out| awk -F":" '{print $1}')
head -$((b-1)) out> parameters.dat
rm out
