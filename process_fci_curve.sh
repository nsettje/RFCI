#!/bin/bash
if [ "$#" -eq 2 ]; then
mol=$1;
bas=$2;
energy_dir=molecule/$1/basis/$2/data/FCI/energy
curve_file=$energy_dir/FCIcontour_$1_$2.txt
touch $curve_file
rm $curve_file
touch $curve_file
for file in $energy_dir/*energy*
do
	rxncoord=$(echo $file | sed s:${energy_dir}/FCIenergy_${mol}_${bas}_::)
	fci_energy=$(sed -n 1p $file)
	echo $rxncoord $fci_energy >> $curve_file
done
cat $curve_file
fi
