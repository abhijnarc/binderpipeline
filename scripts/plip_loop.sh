#!/bin/bash
# script to run PLIP over a loop
PLIP=/home/user/miniconda3/envs/plip/bin/plip
PLIPDIR=./plip/
# run a loop over clone list
for pref in `cat ../clone_list`
do
	echo "STATUS: Analyzing interface for $pref.."
	# file names
	pdb=${pref}_0002.pdb
	# run PLIP
	${PLIP} -f ${pdb} --out ${PLIPDIR} --txt --xml \
	    --name $pref  --chains '[['A'],['B']]'
	echo "STATUS: Completed interface analysis for $pref.."
done
