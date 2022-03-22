#!/bin/bash
#PBS -l nodes=1:ppn=6,vmem=30g,mem=30g,walltime=5:00:00

module load samtools/1.10

cd $PBS_O_WORKDIR

samtools view -@ 6 -hb --no-PG --write-index -o ${file}.reduced.bam##idx##${file}.reduced.bai ${file}

