#!/bin/bash
if [ "$#" -eq 5 ]; then
./clean_directory.sh $1 $2
mol=$1;
bas=$2;
low=$3;
hi=$4;
steps=$5;
step_size=($hi-$low)/$steps;
for ((i=0;i<=$steps-1;i++)); do
	dist=$(echo "scale=4; $low + $i*$step_size" | bc | awk '{printf "%f",$0}')
	./CI.sh $1 $2 $dist
done
else
echo "Input must be of the form: MOLECULE BASIS  LOWER_BOUND HIGHER_BOUND NUM_STEPS"
fi
