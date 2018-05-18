#This script creates a new interaction in the copyint the contet on E_0.5_S_1.0i
#The inputs are:
#the new Epsilon new Sigma.

Epsilon=$1
Sigma=$2


name=$(echo "E_${Epsilon}_S_${Sigma}")
echo $name
cp -r Template $name

cd $name

#Creating the interaction file
sed 's/$Epsilon/'$Epsilon'/g' in.interaction>output
sed 's/$Sigma/'$Sigma'/g' output>in.interaction
rm output
#Equilibration qsub
sed 's/$Flag/'$name'/g' dSetBox/dSetBox.qsub1 >dSetBox/dSetBox.qsub
rm dSetBox/dSetBox.qsub1
