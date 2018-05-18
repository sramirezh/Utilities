#This script moves all the radial and tip behavior plots and classifies them by interactions and
#forces

#Due to the depths analysed, this has to be run in the main project folder, i.e where
#The files N_1  N_30  N_90 are.

Files_tip=`find .  -name "tip_behaviour.pdf"  -path "*/plots*"`
Files_radial=`find .  -name "radial_distribution.pdf"  -path "*/plots*"`
Files_cm=`find .  -name "xcm.pdf"  -path "*/plots*"`
#SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )" #To get the absolute path of this script
CurrentPath=$(pwd)

dir_force=plot_force
dir_interaction=plot_interaction

rm -r $dir_force 2>/dev/null
rm -r $dir_interaction 2>/dev/null

mkdir -p $dir_force
mkdir -p $dir_interaction

echo This path is set as the starting directory
echo $CurrentPath

#Organizing the tip files
echo "Organizing the tip files"
for file in ${Files_tip}; do

  path_file="$(dirname "$file")"
	sim_type=$(echo $file | awk -F  "/" '{print $(NF-5)}') 
 	Number=$(echo $file | awk -F  "/" '{print $(NF-4)}')
	Interaction=$(echo $file | awk -F  "/" '{print $(NF-3)}')
	Force=$(echo $file | awk -F  "/" '{print $(NF-2)}')
	force_tag=$(echo "${sim_type}_${Number}_${Interaction}")
	interaction_tag=$(echo "${sim_type}_${Number}_${Force}")
	mkdir -p $CurrentPath/$dir_force/$Force/tip
	mkdir -p $CurrentPath/$dir_interaction/$Interaction/tip
	cp $file $CurrentPath/$dir_force/$Force/tip/${force_tag}.pdf
	cp $file $CurrentPath/$dir_interaction/$Interaction/tip/${interaction_tag}.pdf

done


#Organizing the radial distribution files
echo "Organizing the radial files"
for file in ${Files_radial}; do

  	path_file="$(dirname "$file")"
 	sim_type=$(echo $file | awk -F  "/" '{print $(NF-5)}') 
	Number=$(echo $file | awk -F  "/" '{print $(NF-4)}')
	Interaction=$(echo $file | awk -F  "/" '{print $(NF-3)}')
	Force=$(echo $file | awk -F  "/" '{print $(NF-2)}')
	force_tag=$(echo "${sim_type}_${Number}_${Interaction}")
	interaction_tag=$(echo "${sim_type}_${Number}_${Force}")
	mkdir -p $CurrentPath/$dir_force/$Force/rd
	mkdir -p $CurrentPath/$dir_interaction/$Interaction/rd
	cp $file $CurrentPath/$dir_force/$Force/rd/${force_tag}.pdf
	cp $file $CurrentPath/$dir_interaction/$Interaction/rd/${interaction_tag}.pdf

done



#Organizing the displacements
echo "Organizing the cm files"
for file in ${Files_cm}; do

  	path_file="$(dirname "$file")"
 	sim_type=$(echo $file | awk -F  "/" '{print $(NF-5)}') 
	Number=$(echo $file | awk -F  "/" '{print $(NF-4)}')
	Interaction=$(echo $file | awk -F  "/" '{print $(NF-3)}')
	Force=$(echo $file | awk -F  "/" '{print $(NF-2)}')
	force_tag=$(echo "${sim_type}_${Number}_${Interaction}")
	interaction_tag=$(echo "${sim_type}_${Number}_${Force}")
	mkdir -p $CurrentPath/$dir_force/$Force/cm
	mkdir -p $CurrentPath/$dir_interaction/$Interaction/cm
	cp $file $CurrentPath/$dir_force/$Force/cm/${force_tag}.pdf
	cp $file $CurrentPath/$dir_interaction/$Interaction/cm/${interaction_tag}.pdf

done
