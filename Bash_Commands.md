# Disk usage in my home directory
du -sh ~/  

# To merge two pdf files, file1.pdf and file2.pdf:

pdftk file1.pdf file2.pdf cat output mergedfile.pdf

# To see the processes with a given name

pgrep, pkill - look up or signal processes based on name and other

for instance pgrep tmux


# To rsync only one type of file, in this case only the file called poly.atom in all the directories

rsync -avz  --include="*/" --include="poly.atom" --exclude="*" "$from" "$to" --progress --partial  . ~/Trajectories
