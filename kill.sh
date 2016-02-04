#!/bin/bash
count1=1	#keep at one!
count2=5273872	#first job number (look with qstat)
 for count1 in {1..200} 	#number of jobs to delete
 do
qdel $count2.cream02
    count1=$((count1+1))
    count2=$((count2+1))
done
exit 0
