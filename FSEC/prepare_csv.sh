#!/bin/sh
#0.01 HM 150330

if [ ! -f "$1" ]
then
    echo "File not found!"
    echo "Please specify a file"    
exit
else
echo "Processing $1 ...."
fi
pdf=${2:-TRUE} #checks if the second argument is present if not sets the variable to TRUE
echo "Printing pdf: "$pdf
a=1
time=0
intensity=0
#Count the number of columns
columns=$(awk -F',' '{print NF; exit}' "$1")
runnumber=$((columns/8))
echo "Number of SEC runs: "$runnumber
string=""
#loop over the individual runs taking the the 7th and 8th columns which contain the fluorescence signal and the time
while [ $a -le $runnumber ]
do
   time=$((8*${a}-1))
   intensity=$((8*${a}))
   #string+=$time","$intensity","
   string=$string$time","$intensity","
   a=`expr $a + 1`
done
#Removing all other colums
echo "Cutting columns ...."
cut -d ',' -f ${string%?} "$1" > trimmed.csv
#Removing rows with missing values
echo "Removing commas ...."
grep -v ",," trimmed.csv>cleaned.csv
#Remove the 2nd line to avoid troubles in R later on
echo "Removing the 2nd line ...."
# Line to output the result of the file 
#sed '2d' cleaned.csv>${1%????}Ready.csv
sed '2d' cleaned.csv>Ready.csv
rm cleaned.csv
rm trimmed.csv
#Calling the R scripts which expects a file called Ready.csv as input
echo "Calling the R script\n"
Rscript SEC_curve_area_working.R $pdf
rm Ready.csv
