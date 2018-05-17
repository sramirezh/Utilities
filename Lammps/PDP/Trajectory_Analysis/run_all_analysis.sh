#Now it can be run anywher as it detects all the poly.atom files
List_files=`find . -name "poly.atom"`
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )" #To get the absolute path of this script
CurrentPath=$(pwd)
echo This path is set as the starting directory
echo $CurrentPath
for file in $List_files; do
	echo $file
	cd "$(dirname "$file")"
	python $SCRIPTPATH/poly_analysis.py poly.atom --split True --Binsize 0.375
	rm *.cxyz 2>/dev/null

	cd $CurrentPath
done
