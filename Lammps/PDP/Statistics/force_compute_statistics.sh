#Runs dp_poly and FastAvearger.py for every force
#It has to be run inside E*
	for force in dD*/;
		do
		cd $force
		/nodescratch/frenkelscratch/sr802/DiffusioP/programs/dp_poly -s 100000 -n 10000000
		python ~/Utilities/Others/Statistics/FastAverager.py vdata.dat
		cd ..
	done
