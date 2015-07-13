#!/bin/bash
if [ "$#" -eq 2 ]; then
mkdir -p molecule/$1/basis/$2/data/FCI/wfn/old
mv -f molecule/$1/basis/$2/data/FCI/wfn/*wfn* molecule/$1/basis/$2/data/FCI/wfn/old
mkdir -p molecule/$1/basis/$2/data/FCI/energy/old
mv -f molecule/$1/basis/$2/data/FCI/energy/*energy* molecule/$1/basis/$2/data/FCI/energy/old
mkdir -p molecule/$1/basis/$2/data/RFCI/energy/old
mv -f molecule/$1/basis/$2/data/RFCI/energy/*energy* molecule/$1/basis/$2/data/RFCI/energy/old
mkdir -p molecule/$1/basis/$2/data/RFCI/wfn/old
mv -f molecule/$1/basis/$2/data/RFCI/wfn/*wfn* molecule/$1/basis/$2/data/FCI/wfn/old
mkdir -p molecule/$1/basis/$2/data/RFCI/tables/old
mv -f molecule/$1/basis/$2/data/RFCI/tables/*$1* molecule/$1/basis/$2/data/RFCI/tables/old
else
echo "Need two arguments: MOLECULE BASIS"
fi
