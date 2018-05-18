#Copies the structure and useful files from one system
#The inputs are the origin_folder and final_folder
origin_folder=$1
final_folder=$2
rsync -avz --exclude 'plots' --exclude 'restart.*' --exclude '*.dat' --exclude '*.gz' --exclude '*.atom' --exclude '*.o*' --exclude '*.lammps' --exclude '*.pbs' --progress --partial $origin_folder $final_folder
