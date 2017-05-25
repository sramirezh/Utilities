#Generates the splitted configurations, if one wants to focus only in the first cluster, uncomment the line 17 and comment the next to get only the first desired cluster, using the variables in 2 and 3.

#The input is the name of the bond list to split

name=$1

rm *.bl

rm Times.dat

####Generating the Splitted Bonded List (The problem here is that the number of lines is variable)
csplit --digits=4 -z --quiet --prefix=outfile $name "/#/" "{*}" #Splits between things "#"

for f in outfile*
do
var=`head -1 $f|awk {'print $3'}` #gets the time
echo $var>>Times.dat
tail -n +2 $f>$var.bl
done

rm outfile*
