#!/bin/bash

module load parallel/20210322

# make it a function
generate_csv_from_filepath(){
     echo $1 | perl -ne '
             @a = split "/", $_;
             if ($a[scalar(@a)-1] =~ m/^(.+?)_([S|ATGC].*?)_(L.*)_(R[12])_[0-9]*\.fastq\.gz/){
                 print $1 . "," . $4 . "," . $3 . "," . $2 .  "," . $_;
             } 
             elsif ($a[scalar(@a)-1] =~ m/^(.+?)_([S|ATGC].*?)_(R[12])_([0-9]*)\.fastq\.gz/){
                 print $1 . "," . $3 . "," . $4 .  "," . $2 . "," . $_;
             }
             elsif ($a[scalar(@a)-1] =~ m/^(.+?)_(R[12])\.fastq\.gz/){ 
                 print $1 . "," . $2 . "," . "L001" . "," . "S1" . "," . $_;
             } 
             else { 
                 print ",,,," . $_ 
             }'
}
export -f generate_csv_from_filepath

predict_TN(){
    patients=(`cat $1 | cut -d, -f1 | grep "M[DM]R*[0-9]*[TB][0-9]*" | sed 's/[TB][0-9]*//' | sort -u`)
    for patient in ${patients[@]}; do
        for tumor in $(cat $1 | cut -d, -f1 | grep -o "${patient}T[0-9]*" | uniq); do
            normals=(`cat $1 | cut -d, -f1 | grep -o "${patient}B[0-9]*" | uniq`)
            if [[ "${#normals[@]}" == 1 ]]; then
                echo "${tumor},${normals[0]}"
            elif [[ "${#normals[@]}" -gt 1 ]]; then
                for normal in ${normals[@]}; do
                    echo "${tumor},${normal}"
                done
            else
                echo ${tumor},NA
            fi
        done
    done
}
export -f predict_TN

# check if file exists
if [[ -e file_list.csv ]]; then
    read -p "file_list.csv exists. Do you wish to overwrite? [y/n]:" input1
    if [[ "${input1}" != "y" ]]; then
        read -p "skip to T/N guessing? [y/n]:" input2
        if [[ "${input2}" == "n" ]]; then
            echo "exiting..."
            exit 0
        fi
    fi
else
    input1="y"
fi

# check arguments and execute
if [[ ${input1} == "y" ]]; then
    if [[ -z $1 ]]; then
        echo "provide path to files as argument list. Can use '*' wildcards"
        exit 1
    elif [[ "$#" -gt 1 ]]; then
        # sort csv naturally
        parallel --keep 'generate_csv_from_filepath {}' ::: $* | sort -V | paste -d, - - | cut -d, -f1,3-5,10 > file_list.csv 
    else
        if [[ -e $1 ]]; then
            file $1 | grep "gzip compressed data"
            if [[ "$?" == 0 ]]; then
                generate_csv_from_filepath $1
            fi
            file $1 | grep "ASCII text"
            if [[ "$?" == 0 ]]; then
               cat $1 | parallel --keep 'generate_csv_from_filepath {}' | sort -V | paste -d, - - | cut -d, -f1,3-5,10 > file_list.csv
            fi
        fi
    fi
fi

# report
echo "file_list.csv is ready"

# check file
grep "^,,," file_list.csv &> /dev/null
if [[ "$?" == 0 ]]; then
    echo "check file_list.csv; file not properly constructed"
    exit 1
fi

# attempt to predict the T/N csv file
predict_TN file_list.csv > tumors_and_normals.csv
# check file
file tumors_and_normals.csv
if [[ "$?" == 0 ]]; then
    echo "tumors_and_normals.csv is ready"
fi

# check if T/N csv file is ready to go
total_samples=$(cat file_list.csv | cut -d, -f1 | sort -u | wc -l)
samples_in_TN=$(cat tumors_and_normals.csv | tr ',' '\n' | grep -v "^NA$" | wc -l)

if [[ "${total_samples}" == "${samples_in_TN}" ]]; then
    echo "ready to start pipeline"
else
    echo "revise/edit tumors_and_normals.csv manually before starting"
fi
