#!/bin/bash

#SBATCH -n 1 --mem-per-cpu=16G -t 05:00:00
#SBATCH --dependency=afterok:<JOBID>
module load r

for f in *.depth
do
	Rscript downsample.R $f
	gzip $f
done
