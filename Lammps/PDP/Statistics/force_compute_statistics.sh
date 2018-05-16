#Runs dp_poly and FastAvearger.py for every force
#It has to be run inside E*
Interactions=$@
CurrentPath=$(pwd)
for interaction in $Interactions;
	do
	
	echo $interaction
	cd $interaction
	for force in dD*/;
		do
		echo $force
		cd $force
		/nodescratch/frenkelscratch/sr802/DiffusioP/programs/dp_poly -s 100000 -n 10000000
		python ~/Utilities/Others/Statistics/FastAverager.py vdata.dat
		cd ..
	done
	cd $CurrentPath
done
