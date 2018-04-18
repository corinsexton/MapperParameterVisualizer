#!/bin/bash

#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH -t 24:00:00
#SBATCH -J "ecoli mapq"
#SBATCH --array=4-31

module load samtools jdk

for f in bamfiles/${SLURM_ARRAY_TASK_ID}*.bam
do
	python viewMappingQuality.py $f
done
