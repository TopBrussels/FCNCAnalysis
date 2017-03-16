#!/bin/bash 

cd output

for f in ../submit*.sh
 do
	echo "qsub " $f >> test.txt 
 done

cd -

