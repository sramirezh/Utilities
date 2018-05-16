#Due to the depths analysed, this has to be run in the main project folder, i.e where
#The files N_1  N_30  N_90 are.
#_LISPDIRS=`find . -maxdepth 3 -mindepth 3 -type d | sort`
Files_tip=`find .  -name "tip_behaviour.pdf"  -path "*/E_*"`
Files_radial=`find .  -name "radial_distribution.pdf"  -path "*/E_*"`
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )" #To get the absolute path of this script
CurrentPath=$(pwd)
dirname=plot_force
mkdir  $dirname 2>dev/null
echo This path is set as the starting directory
echo $CurrentPath
for file in ${Files_tip}; do
	echo $file
 	Number=$(echo $file | awk -F  "/" '{print $2}')
	Interaction=$(echo $file | awk -F  "/" '{print $3}')
	Force=$(echo $file | awk -F  "/" '{print $4}')
	tag=$(echo "${Number}_${Interaction}")
	echo $Number
	echo $Force
#	mkdir -p $CurrentPath/$dirname/$Force/tip 2>dev/null
	#cp $file $CurrentPath/$dirname/$Force/tip/${tag}.pdf 
done
