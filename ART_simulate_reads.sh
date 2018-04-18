#!/bin/bash

#SBATCH --time=48:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32G   # memory per CPU core
#SBATCH -J "simulate ecoli"   # job name

./art_illumina -ss HS25 -sam -i ~/compute/model/miseq_ecoli/ecoli.fasta -p -l 150 -f 200 -m 200 -s 10 -o ~/compute/model/ecoli_sim_200x
