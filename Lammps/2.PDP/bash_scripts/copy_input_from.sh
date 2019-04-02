#Run inside iaxx
#copies all the input file from the E_0.5_S_0.8, which should be the 
#More updated
for folder in d*/ ; do
	cd $folder
	echo $folder
	cp ../../E_0.5_S_0.8/$folder/input.lmp .
	cd ..
done
