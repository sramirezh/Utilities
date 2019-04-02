#runs python statistic_parameters.py in read mode in all folders

folders=`find . -type d -name "N_*"`
CurrentPath=$(pwd)
for folder in $folders;
do
echo "Running inside  ${folder}"
cd $folder
python ~/Utilities/Lammps/PDP/statistics/statistics_parameters.py --source READ
cd $CurrentPath
done
