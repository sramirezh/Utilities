#This script gets all the statistical properties by running dp_poly and FastAvearger.py for every force
#It has to be run inside every N_X

CurrentPath=$(pwd)
rm Statistic_summary.dat 
touch Statistic_summary.dat
for f in */; 
	do
	echo "############################################################################" >>$CurrentPath/Statistic_summary.dat
	echo $f >>$CurrentPath/Statistic_summary.dat
	cd $f
	for force in dD*/;
		do
		echo " " >>$CurrentPath/Statistic_summary.dat
		echo $force >>$CurrentPath/Statistic_summary.dat
		cd $force
		#/nodescratch/frenkelscratch/sr802/DiffusioP/programs/dp_poly -s 100000 -n 10000000
		cat average_info.dat >>$CurrentPath/Statistic_summary.dat
		cat statistics.dat >>$CurrentPath/Statistic_summary.dat
		#python ~/Utilities/Others/Statistics/FastAverager.py vdata.dat 
		cd ..		
	done
	cd $CurrentPath

done
	




