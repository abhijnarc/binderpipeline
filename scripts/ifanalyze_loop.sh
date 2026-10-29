#!/bin/bash
ROSETTA=/home/user/tools/rosetta.binary.ubuntu.release-408/main/source/bin
IFXML=/data/scripts/bench20/interface_AB.xml
#
# run a loop over clone list
for pref in `cat ../clone_list`
do
	echo "STATUS: Analyzing interface for $pref.."
	# file names
	pdb=${pref}_0002.pdb
	score=./ifanalyze/${pref}_ifmetric.sc
	# run rosetta interface analyze
	${ROSETTA}/rosetta_scripts.default.linuxgccrelease \
    		-parser:protocol ${IFXML} \
    		-s ${pdb} \
    		-out:file:score_only ${score}
	echo "STATUS: Completed interface analysis for $pref.."
done
