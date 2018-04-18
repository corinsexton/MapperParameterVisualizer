# MapperParameterVisualizer
Executing `rshinyApp.R` begins the Mapper Parameter Visualizer app. This a general description of the pipeline I followed:
 1. Create the alignments using the scripts in `bowtie2` and `bwamem` directories
 2. Execute scripts in `analysis_scripts`
 3. Run `rshinyApp.R` to visualize data

### Description of Scripts
##### `ART_simulate_reads.sh`
script run to simulate Illumina paired-end reads using ART

##### `generate_figures.R`
script to generate all figures within the paper

##### `rshinyApp.R`
code for RShiny application

##### `bwamem/bwa_miseq_run.sh` and `bwamem/bwa_simulated_run.sh`
scripts for generating the BWA-MEM alignments with all different parameters

##### `bowtie2/bowtie2_miseq_run.sh` and `bowtie2/bowtie2_simulated_run.sh`
scripts for generating the Bowtie 2 alignments with all different parameters

#### Analysis Directory
For all scripts in the `analysis_scripts` directory, a primary script written in R or python is included as well as a bash script which shows how to execute the primary script. They are correspondingly named.

##### `checkAccuracy`
script used to compare the simulated read alignments to the correct ART alignment (written in python)

##### `downsample`
script used to downsample the depth of the reads (written in R)
average is taken across 10,000 base pair windows

##### `viewMappingQuality`
script used to isolate the mapping quality distributions of all reads (written in python, corresponds to mapq.sh)
