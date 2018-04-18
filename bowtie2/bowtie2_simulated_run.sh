#!/bin/bash

#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH -t 08:00:00
#SBATCH -J "bwa sim picard"
#SBATCH --array=4-31
#SBATCH -o bowtie_%a_sim.out

module load bowtie2 samtools jdk

#for L in {4..31}
#do
	for N in {0..1}
	do
		for D in {5..25..5}
		do
			for R in {1..3}
			do
				outfile_prefix=sim_${SLURM_ARRAY_TASK_ID}_${N}_${D}_${R}
				echo "evaluating $outfile_prefix"
				
				time bowtie2 -L ${SLURM_ARRAY_TASK_ID} -N $N -D $D -R $R -x ecoli -1 ../ecoli_sim_200x1.fq -2 ../ecoli_sim_200x2.fq --threads 16 | samtools view -@ 16 -b - | samtools sort -@ 16 - > bamfiles/$outfile_prefix.bam
				samtools index bamfiles/$outfile_prefix.bam
				
				samtools depth -a bamfiles/$outfile_prefix.bam > depth/$outfile_prefix.depth
				samtools depth -Q 42 -a bamfiles/$outfile_prefix.bam > depth/${outfile_prefix}_42.depth

				java -jar picard.jar CollectAlignmentSummaryMetrics I=bamfiles/$outfile_prefix.bam O=picard_stats/$outfile_prefix.picard_stats R=ecoli.fasta
				echo
				echo
			done
		done
	done
#done
