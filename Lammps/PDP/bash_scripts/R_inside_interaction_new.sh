#Run inside every iaxx
#Copies the mt.data from every directory and

for folder in dDP*/ ; do
	cd $folder
	echo $folder
	qsub run.qsub
	cd ..
done
