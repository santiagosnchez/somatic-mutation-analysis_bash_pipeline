#!/bin/bash

#### merges multiple working directories into one ####

# functions
make_fs_skeleton(){
  mkdir -p $1/bam $1/all_logfiles $1/analyses $1/contamination $1/mutect2/f1r2 $1/vcf/snpEff
}

make_md5sum(){

}

# merges single files into one
merge_single_files(){
  args=($@)
  file_to_merge=$1
  unset 'args[0]'
  cat $(for i in ${args[@]}; do echo $i/$file_to_merge; done)
}

# check if runs in each completed
check_complete_dir(){
  for dir in $*; do
    tail -1 $dir/main.log | grep "pipeline took" &> /dev/null
    if [[ "$?" != 0 ]]; then
      echo "Run in $dir did not complete. Aborting..."
      return 1
    fi
  done
}

# store args in array
dirs=($@)

# this block deals with the initial directory
# where files will be placed
# first check if arguments are empty
if [[ $# == 0 ]]; then
    echo "Specify the directories that need to be merged as arguments."
    echo "Optionally, the first argument could be the destination directory"
else
    if [[ "$#" -gt 1 ]]; then
        # first argument is the name of the new output dir
        if [[ -e $1 ]]; then
            echo "${1} exists. Do you want to replace? [y/n]:"
            read -r response
            if [[ "${response}" == "y"* ]]; then
                # check again
                echo "Are you sure? Files in ${1} will be lost. [y/n]:"
                read -r response
                if [[ "${response}" == "y"* ]]; then
                     # delete and make skeleton fs
                    rm -rf $1
                    mkdir $1
                    make_fs_skeleton $1
                    # set out dir var
                    out_dir=$1
                    # remove first element from array
                    unset 'dirs[0]'
                elif [[ "${response}" == "n"* ]]; then
                    echo "Exiting..."
                    exit 0
                fi
            elif [[ "${response}" == "n"* ]]; then
                echo "Do you want to create a new directory? [y/n]:"
                read -r response
                if [[ "${response}" == "y"* ]]; then
                    echo "What name? [make up a name]:"
                    read -r response
                    if [[ "${#response}" -gt 0 ]]; then
                        # make new dir, fs, and keep dirs in arguments
                        mkdir ${response}
                        make_fs_skeleton ${response}
                        # set out dir var
                        out_dir=${response}
                    else
                        echo "No name provided..."
                        exit 1
                    fi
            else
                echo "Unrecognized command..."
                exit 1
            fi

        else
            # make skeleton
            mkdir $1
            make_fs_skeleton $1
            # set out dir var
            out_dir=$1
            # remove first element from array
            unset 'dirs[0]'
        fi
    else
        echo "Single argument not allowed..."
        exit 1
    fi
fi

# check if the pipeline completed in each dir
check_complete_dir ${dirs[@]}
if [[ "$?" == 0 ]]; then
    # proceed
    # move files to destination dir
    # merge single files
    for f in file_list.csv tumors_and_normals.csv; do
        merge_single_files $f ${dirs[@]} | sort -u > $out_dir/$f
    done
    # do TBM and coverage separately
    # get head from first dir
    head -1 ${dirs[0]}/analyses/coverage_and_tmb.csv > $out_dir/analyses/head
    merge_single_files analyses/coverage_and_tmb.csv ${dirs[@]} | grep -v "obs_coverage," | cat $out_dir/analyses/head - > $out_dir/analyses/coverage_and_tmb.csv
    rm $out_dir/analyses/head
    # keep each main logfile as unique
    for d in ${dirs[@]}; do
        cp $d/main.log "$out_dir/main_${d}.log"
    done
    # safely move all files in these dirs to out skeleton
    # all these should have unique file names
    # other wise the latest version of a file will be kept
    for sd in all_logfiles bam vcf vcf/snpEff contamination mutect2 mutect2/f1r2; do
        for d in ${dirs[@]}; do
            mv -u $( find ${d}/${sd} -type f ) ${out_dir}/${sd}
        done
