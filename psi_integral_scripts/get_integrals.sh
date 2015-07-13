#!/bin/bash
echo "Molecule: $1"
echo "Basis: $2"
echo "Rxn Coord: $3"
if [ -e "molecule/$1/basis/$2/integrals/OEI_$1_$2_$3" -a -e "molecule/$1/basis/$2/integrals/TEI_$1_$2_$3" -a -e "molecule/$1/basis/$2/elec_constants/econst_$1_$2_$3" ]
then
    echo "Integrals and constants exist!"
else
    mkdir -p molecule/$1/basis/$2/integrals
    mkdir -p molecule/$1/basis/$2/elec_constants
    echo "Integrals and constants DO NOT exist!"
    cd ~/psi_plugins/get_integrals
    ./GET_INTEGRALS.sh $1 $2 $3
    cd ~/RFCI
fi
