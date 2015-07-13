#!/bin/bash
molecule=h2o
basis="STO-3G"
for file in systems/${molecule}/basis/${basis}/inputs/*.dat
do
		echo $file
		radius=$(grep "#R = " $file | sed -e s:"#R = "::g | tr -d ' ')
		jobname="${molecule}_${basis}_${radius}_mdfci"
		echo $jobname
		sed -e s:nnnnn:$jobname:g \
		-e s:iiiii:${file}:g \
		-e s:ooooo:"systems/${molecule}/basis/${basis}/outputs/output_${molecule}_${basis}_${radius}.dat":g sub_.sh > sub.sh
		qsub sub.sh 
		sleep 150s
done
