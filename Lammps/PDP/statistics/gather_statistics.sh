#This script gets all the statistical properties by running dp_poly and FastAvearger.py for every force
#It has to be run inside every N_X

CurrentPath=$(pwd)
rm Statistic_summary.dat 2>/dev/null
touch Statistic_summary.dat
for f in E_*/;
	do
	echo "############################################################################" >>$CurrentPath/Statistic_summary.dat
	echo $f >>$CurrentPath/Statistic_summary.dat
	cd $f
	for force in dD*/;
		do
		echo " " >>$CurrentPath/Statistic_summary.dat
		echo $force >>$CurrentPath/Statistic_summary.dat
		cd $force
		cat average_info.dat >>$CurrentPath/Statistic_summary.dat
		tail -n +2 statistics.dat >>$CurrentPath/Statistic_summary.dat
		cd ..
	done
	cd $CurrentPath

done
