#Due to the depths analysed, this has to be run in the main project folder, i.e where
#The files N_1  N_30  N_90 are.
_LISPDIRS=`find . -maxdepth 3 -mindepth 3 -type d | sort`
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )" #To get the absolute path of this script
CurrentPath=$(pwd)
dirname=plot_force
mkdir $dirname
echo This path is set as the starting directory
echo $CurrentPath
for dir in ${_LISPDIRS}; do
    cd ${dir}
 	Number=$(echo $dir | awk -F  "/" '{print $2}')
	Interaction=$(echo $dir | awk -F  "/" '{print $3}')
	Force=$(echo $dir | awk -F  "/" '{print $4}')
	tag=$(echo "${Number}_${Interaction}")
	mkdir -p $CurrentPath/$dirname/$Force/tip
	mkdir -p $CurrentPath/$dirname/$Force/rd
	cp tip_behaviour.pdf $CurrentPath/$dirname/$Force/tip/${tag}.pdf 
	cp radial_distribution.pdf  $CurrentPath/$dirname/$Force/rd/${tag}.pdf
    cd $CurrentPath
done
