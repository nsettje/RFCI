#!/bin/bash
if [ "$#" -eq 3 ]; then
./psi_integral_scripts/get_integrals.sh $1 $2 $3
./psi_integral_scripts/generate_RFCI_input.sh $1 $2 $3
mkdir -p molecule/$1/basis/$2/inputs
mkdir -p molecule/$1/basis/$2/outputs
mkdir -p molecule/$1/basis/$2/data/FCI/wfn
mkdir -p molecule/$1/basis/$2/data/FCI/energy
mkdir -p molecule/$1/basis/$2/data/RFCI/energy
mkdir -p molecule/$1/basis/$2/data/RFCI/wfn
mkdir -p molecule/$1/basis/$2/data/RFCI/tables
./RFCI molecule/$1/basis/$2/inputs/input_$1_$2_$3.txt molecule/$1/basis/$2/outputs/output_$1_$2_$3.txt 
echo "Final RFCI energies"
cat molecule/$1/basis/$2/data/RFCI/energy/RFCIenergy_$1_$2_$3
else
echo "Need exactly three arguments!"
fi
