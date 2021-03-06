#!/bin/bash

#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH -t 18:00:00
#SBATCH -J "bwa simulated"
#SBATCH --array=4-31
#SBATCH -o bwamem_%a_sim.out

module load bwa samtools jdk

#for k in {4..31}
#do
	seq 0.5 0.5 3.0 | while read r
	do	
		outfile_prefix=${SLURM_ARRAY_TASK_ID}_${r}
		echo "evaluating $outfile_prefix"
		
		time bwa mem -t 16 -k $SLURM_ARRAY_TASK_ID -r $r ecoli.fasta ../ecoli_sim_200x1.fq ../ecoli_sim_200x2.fq | samtools view -@ 16 -b - | samtools sort -@ 16 - > bamfiles/$outfile_prefix.bam
		samtools index bamfiles/$outfile_prefix.bam
		
		samtools depth -a bamfiles/$outfile_prefix.bam > depth/$outfile_prefix.depth
		samtools depth -Q 60 -a bamfiles/$outfile_prefix.bam > depth/${outfile_prefix}_42.depth

		java -jar picard.jar CollectAlignmentSummaryMetrics I=bamfiles/$outfile_prefix.bam O=picard_stats/$outfile_prefix.picard_stats R=ecoli.fasta

		echo
		echo
	done
#done

