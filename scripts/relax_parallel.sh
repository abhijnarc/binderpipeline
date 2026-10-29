#!/bin/bash
ROSETTA=/home/user/tools/rosetta.binary.ubuntu.release-408/main/source/bin
# run Rosetta relax in a loop over the five structures
# Run over two complex pdbs
# (You can run this over a set of pdbs from a list file using 
# for pref in `cat list_file`; do
for pref in 1a22 3mca; do
	# run over all models
	for i in {0..4}; do
		# filenames
		pdb=../AF3/AF3_models/${pref}_${i}.pdb
		score=${pref}_${i}_relax.sc
		log=${pref}_${i}_relax.log
		echo "Relaxing ${pdb} ..."
		time ${ROSETTA}/relax.default.linuxgccrelease \
    			-s ${pdb} \
    			-relax:fast \
    			-use_input_sc \
			-constrain_relax_to_start_coords \
    			-relax:coord_constrain_sidechains \
    			-ex1 -ex2 \
    			-nstruct 1 \
			-out:file:scorefile ${score} > ${log} 2>&1 &
	done
done
