#!/bin/bash
#Declare a new array for pdb names 
declare -a StringArray2=("model2" "model55"  "model84")

#Iterate the string array using for loop
for i in ${StringArray2[@]}; do
	/usr/local/lmm-tools/bin/naccess2.1.1/naccess $i.pdb
	echo "NACCESS FOR $i COMPLETED"
done

