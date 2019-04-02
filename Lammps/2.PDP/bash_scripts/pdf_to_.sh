#this script converts all the pdfs inside the directory including subdirectories into the desired format
extension=$1
list_pdf=`find . -name "*.pdf"`
for file in $list_pdf;
do
path_file="$(dirname "$file")"
name=$(echo $file | awk -F  "/" '{print $(NF)}') 
name_noextension=`echo ${name%.*}`
 
echo $file  
final_file=`echo $path_file/$name_noextension.$extension`

convert $file $final_file
done
