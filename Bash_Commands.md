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
