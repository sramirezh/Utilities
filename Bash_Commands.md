# Disk usage in my home directory
du -sh ~/  

# To merge two pdf files, file1.pdf and file2.pdf:

pdftk file1.pdf file2.pdf cat output mergedfile.pdf

# To see the processes with a given name

pgrep, pkill - look up or signal processes based on name and other

for instance pgrep tmux


# To rsync only one type of file, in this case only the file called poly.atom in all the directories

rsync -avz -e ssh --include="*.pdf" --include="*/" --exclude="*" --progress --partial . nimbus:Foldername


# To find all the files *.pdf in folder and then execute grep to find something else

find .  -name '*.pdf'  -path './E_*' -exec grep 'SENTENCE TO FIND' {} /dev/null \;

also

find .  -name '*.pdf'  ! -path './E_*'  , The use of ! is to avoid that path

# To mount files from a remote server on my computer

mkdir dexter_scratch

sshfs dexter:/frenkelscratch/sr802 dexter_scratch


# Change the walltime of all *.qsub
list=\`find . -name "run.qsub"\`
list=\`find . -name 'Equilibration' -prune -o  -name "run.qsub"\` To exclude inside Equilibration

for file in $list; do sed -i '3s/.*/#PBS -l walltime=24:00:00/' $file;done 


# Delete all the tasks for a user
qselect -u \<username> | xargs qdel

# Change the bridghtness

xrandr -q | grep " connected" #Shows the connected monitors

generates this output:
DVI-I-0 connected primary 1920x1080+0+0 (normal left inverted right x axis y axis) 480mm x 270mm
DVI-D-0 connected 1680x1050+1920+0 (normal left inverted right x axis y axis) 459mm x 296mm

xrandr --output DVI-I-0 --brightness 0.5 #Changes DVI-I-0 to 50%

# Find files with regex

find . -type f -regextype sed -regex ".*.[0-9]$"   This was to delete the output from old simulations

# Find the 5 largest files

find -type f -exec du -Sh {} + | sort -rh | head -n 5

