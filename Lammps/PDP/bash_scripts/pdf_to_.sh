#this script converts all the pdfs inside the directory including subdirectories into the desired format
extension=$1
list_pdf=`find . -name "*.pdf"`
for file in $list_pdf;
do
path_file="$(dirname "$file")"
name=$(echo $file | awk -F  "/" '{print $(NF)}') 
name_noextension=`echo ${name%.*}`
echo $name_noextension 
convert $file  $path_file/$name_noextension.jpg
done
