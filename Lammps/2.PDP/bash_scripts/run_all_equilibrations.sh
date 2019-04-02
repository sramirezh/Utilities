#runs python statistic_parameters.py in read mode in all folders

folders=`find . -type d -name "Equilibration"`
CurrentPath=$(pwd)
for folder in $folders;
do
echo "Running inside  ${folder}"
cd $folder
qsub run.qsub
cd $CurrentPath
done
