#Due to the depths analysed, this has to be run in the main project folder, i.e where
#The files N_1  N_30  N_90 are.
_LISPDIRS=`find . -maxdepth 4 -mindepth 4 -type d | sort`
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )" #To get the absolute path of this script
CurrentPath=$(pwd)
echo This path is set as the starting directory
echo $CurrentPath
for dir in ${_LISPDIRS}; do
    echo -e "$dir"
    cd ${dir} 
    python $SCRIPTPATH/poly_analysis.py poly.atom
	rm *.cxyz 2>/dev/null
    cd $CurrentPath
done
