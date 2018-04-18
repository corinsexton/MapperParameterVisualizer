#!/bin/bash

#SBATCH --mem-per-cpu=32G
#SBATCH --nodes=1
#SBATCH -t 72:00:00
#SBATCH -J "accuracy bowtie simulated"

module load samtools

time python checkAccuracy.py ../ecoli_sim_200x1.aln ../ecoli_sim_200x2.aln bamfiles/*.bam
