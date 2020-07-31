#!/bin/bash

for i in `seq 1 10`;
do
	let num=i;
	argString="--args "$num
	R CMD BATCH "$argString" run_gail_sim.R;
	sleep 12
done