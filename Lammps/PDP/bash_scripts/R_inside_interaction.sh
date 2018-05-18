#Run inside every iaxx
#Copies the mt.data from every directory and 

for folder in */ ; do
	cd $folder
	echo $folder
	cp ../initial_conf/mt.data .
	qsub run.qsub
	cd ..
done
