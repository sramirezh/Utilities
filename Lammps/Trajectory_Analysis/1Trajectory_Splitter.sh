# This script reads the trajectory data and creates just the initial trajectory
# The name of the file to be splitted is given as input in the command line. 
name=$1 
NAtoms=`head -1 $name` #Number of atoms 
Nlines=$(($NAtoms+2))
head -$Nlines trajectory.xyz > 1Trajectory.xyz
