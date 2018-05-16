#This script moves all the radial and tip behavior plots and classifies them by interactions and
#forces

#Due to the depths analysed, this has to be run in the main project folder, i.e where
#The files N_1  N_30  N_90 are.

#_LISPDIRS=`find . -maxdepth 3 -mindepth 3 -type d | sort`
Files_tip=`find .  -name "tip_behaviour.pdf"  -path "*/E_*"`
Files_radial=`find .  -name "radial_distribution.pdf"  -path "*/E_*"`
#SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )" #To get the absolute path of this script
CurrentPath=$(pwd)

dir_force=plot_force
dir_interaction=plot_interaction
mkdir -p $dir_force
mkdir -p $dir_interaction
echo This path is set as the starting directory
echo $CurrentPath

#Organizing the tip files
echo "Organizing the tip files"
for file in ${Files_tip}; do
 	Number=$(echo $file | awk -F  "/" '{print $2}')
	Interaction=$(echo $file | awk -F  "/" '{print $3}')
	Force=$(echo $file | awk -F  "/" '{print $4}')
	force_tag=$(echo "${Number}_${Interaction}")
	interaction_tag=$(echo "${Number}_${Force}")
	mkdir -p $CurrentPath/$dir_force/$Force/tip
	mkdir -p $CurrentPath/$dir_interaction/$Interaction/tip
	cp $file $CurrentPath/$dir_force/$Force/tip/${force_tag}.pdf
	cp $file $CurrentPath/$dir_interaction/$Interaction/tip/${interaction_tag}.pdf

done
#Organizing the radial distribution files
echo "Organizing the radial files"
for file in ${Files_radial}; do
 	Number=$(echo $file | awk -F  "/" '{print $2}')
	Interaction=$(echo $file | awk -F  "/" '{print $3}')
	Force=$(echo $file | awk -F  "/" '{print $4}')
	force_tag=$(echo "${Number}_${Interaction}")
	interaction_tag=$(echo "${Number}_${Force}")
	mkdir -p $CurrentPath/$dir_force/$Force/rd
	mkdir -p $CurrentPath/$dir_interaction/$Interaction/rd
	cp $file $CurrentPath/$dir_force/$Force/rd/${force_tag}.pdf
	cp $file $CurrentPath/$dir_interaction/$Interaction/rd/${interaction_tag}.pdf

done
