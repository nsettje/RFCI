#!/bin/bash
if [ "$#" -eq 2 ]; then
mol=$1;
bas=$2;
energy_dir=molecule/$1/basis/$2/data/RFCI/energy
cont_file=$energy_dir/RFCIcontour_$1_$2.txt
echo $cont_file
touch $cont_file
rm $cont_file
touch $cont_file
for file in $energy_dir/RFCIenergy*
do
	rxncoord=$(echo $file | sed s:${energy_dir}/RFCIenergy_${mol}_${bas}_::)
	fci_energy=$(cat molecule/${1}/basis/${2}/data/FCI/energy/FCIenergy_${1}_${2}_${rxncoord})
	nterms=$(wc -l < ${file})
	for ((i=1;i<=nterms;i++)); do
	rfci_energy=$(sed -n ${i}p $file)
	#energy_diff=$(echo "scale=4; $rfci_energy-$fci_energy" | bc | awk '{printf "%2.6f",$0}')
	#energy_diff=$(echo $energy_diff | tr -d -)
	#echo $rxncoord $energy_diff
	echo ${rxncoord}, ${i}, ${rfci_energy} 
	echo ${rxncoord}, ${i}, ${rfci_energy} >> $cont_file
	done
done
cat $cont_file
else
echo "Input must be of the form: MOLECULE BASIS"
fi
