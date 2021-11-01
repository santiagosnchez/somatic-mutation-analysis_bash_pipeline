#!/bin/bash

module load parallel/20210322

if [[ -z $1 ]]; then
    echo "provide path to files as argument list. Can use '*' wildcards"
    exit 1
fi

# consortium
# file_list
parallel --plus '
name=$(echo {/} | cut -d_ -f1)
lib=$(echo {/} | cut -d_ -f2)
lane=$(echo {/} | cut -d_ -f3)
echo $name,$lane,$lib,{}' ::: /hpf/largeprojects/tabori/WES/EDW13901.20200525/MMR1431*_files/fastq/CE7E6ANXX/*.fastq.gz | paste -d, - - | cut -d, -f1-4,8 > file_list.csv
# tumor and normal
cat file_list.csv | cut -d, -f1 | uniq | paste - - | awk '{print $2","$1 }' > tumors_and_normals.csv


# edit the columns of the file name that you want to keep                    ####################
#for i in $*; do echo $i; done | rev | cut -d'/' -f1 | rev | parallel --plus --keep --colsep="_" echo {10},{16},{15},{17} | tee file_list.csv
#for i in $*; do echo $i; done | paste -d, file_list.csv - | paste -d, - - | cut -d, -f1-2,4,9 > tmp && mv tmp file_list.csv

# generate read group headers for bam
#cat file_list.csv | cut -d, -f1-3 | sort | parallel --keep --colsep="," echo -e "@RG\\\tID:{2}\\\tSM:{1}\\\tLB:{3}\\\tPL:Illumina" > read_groups.txt

