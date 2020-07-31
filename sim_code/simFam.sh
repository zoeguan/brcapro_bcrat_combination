#!/bin/bash

for i in `seq 1 10`;
do
	let num=i;
	argString="--args "$num
	R CMD BATCH "$argString" simFam.R;
	sleep 10
done