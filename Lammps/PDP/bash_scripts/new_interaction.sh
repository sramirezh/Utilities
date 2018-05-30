#This script creates a new interaction in the copyint the contet on E_0.5_S_1.0i
#The inputs are:
#the new Epsilon new Sigma.

Epsilon=$1
Sigma=$2


name=$(echo "E_${Epsilon}_S_${Sigma}")
echo $name
rm -r $name
cp -r Template $name

cd $name

#Creating the interaction file
sed 's/$Epsilon/'$Epsilon'/g' in.interaction>output
sed 's/$Sigma/'$Sigma'/g' output>in.interaction
rm output

#Modifying all the qsub
files=`find . -name "*.qsub"`
for file in $files; 
	do 
	echo 'Modifying' $file
	absolute_path="$( cd "$(dirname "$file")" ; pwd -P )"
	
	Number=$(echo $absolute_path | awk -F  "/" '{print $(NF-3)}')
	Interaction=$(echo $absolute_path | awk -F  "/" '{print $(NF-1)}')
	Sim_type=$(echo $absolute_path | awk -F  "/" '{print $(NF)}')
	name=$(echo "#PBS -N ${Number}_${Epsilon}_${Sigma}_${Sim_type}")
	sed -i "s/#PBS -N.*/${name}/g" $file
done



