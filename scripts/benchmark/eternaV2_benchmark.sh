#!/bin/bash
for ss in `cat ../../data/Eterna100/V2/eterna2_left.csv| cut -f 2 -d ","`; 
	do
		`python aRNAque.py -t $ss  -n 100 -g 20000 -msf 1  --job 20 -bp="GC4" >>../data/eterna/EternaV2_OP3.out`;
	      	echo "solved for target: $ss";

	done
