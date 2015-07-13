#!/bin/bash
if [ "$#" -eq 3 ]; then
	sed -e s:mmmmm:$1:g \
	-e s:bbbbb:$2:g \
	-e s:rrrrr:$3:g \
	input/input.txt > molecule/$1/basis/$2/inputs/input_$1_$2_$3.txt
else
	echo "Not enough arguments!"
fi
