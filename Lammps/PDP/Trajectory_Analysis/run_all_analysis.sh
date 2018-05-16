#Due to the depths analysed, this has to be run in the main project folder, i.e where
#in the parent directory of N=1, N=30, etc
#_LISPDIRS=`find . -maxdepth 4 -mindepth 4 -type d | sort`
List_files=`find . -name "poly.atom"`
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )" #To get the absolute path of this script
CurrentPath=$(pwd)
echo This path is set as the starting directory
echo $CurrentPath
for file in $List_files; do
	echo $file
	cd "$(dirname "$file")"
	python $SCRIPTPATH/poly_analysis.py poly.atom --split True
	rm *.cxyz 2>/dev/null	
	
	cd $CurrentPath
done
