#This script gets all the statistical properties by running dp_poly and FastAvearger.py for every force
#It has to be run inside every N_X
s=$1 #Initial step
d=$2 #Interval
n=$3 #FinalTimestep

CurrentPath=$(pwd)
rm Statistic_summary.dat 2>/dev/null
touch Statistic_summary.dat
for f in E_*/;
	do
	echo "############################################################################" >>$CurrentPath/Statistic_summary.dat
	echo $f >>$CurrentPath/Statistic_summary.dat
	cd $f
	echo $f
	for force in dD*/;
		do
		echo " " >>$CurrentPath/Statistic_summary.dat
		echo $force >>$CurrentPath/Statistic_summary.dat
		cd $force
		/nodescratch/frenkelscratch/sr802/DiffusioP/programs/dp_poly -s $s -d $d -n $n
		cat average_info.dat >>$CurrentPath/Statistic_summary.dat
		python ~/Utilities/Others/Statistics/FastAverager.py vdata.dat
		cat statistics.dat >>$CurrentPath/Statistic_summary.dat
		cd ..
	done
	cd $CurrentPath

done
