# This scripts runs on every folder with results from the simulations and copies the results in a folder called *.total which contains all the results with the number of the run appended.

mkdir A_dir_total
mkdir B_dir_total
mkdir stress_dir_total
rm times.dat
for f in `ls -d RUN.*|sort -V`; # -d flag is to list only directories
	do 
	echo $f
	Nrun=`echo $f |cut -d '.' -f2` #To get the number of the run
	echo $Nrun >> times.dat
	cd $f
	for sf in `ls -d *_dir`;
		do
		dir=`echo $sf`
		dir+="_total"
		cp $sf/full-mean ../$dir/full-mean.$Nrun 
		done
	cd ..
done
