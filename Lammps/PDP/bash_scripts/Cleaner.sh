#Run in the main Monomer, Polymer folder
#Runs Hiroaki cleaner in all the folders
for folder in */ ; do
	cd $folder
	echo $folder
	bash clean_all.sh
	cd ..
done
