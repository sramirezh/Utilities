#This script is used to change the prefix of the files 
#the input variables are the prefix and then the names of the files.

Prefix=$1
for file in "$@"
	do
	if [ "$file" = "$Prefix" ]; then 
	echo "Prefixing..."
	else
	name="$Prefix"_"$file"
	mv $file $name
	fi
	done

