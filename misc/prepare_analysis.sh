#!/bin/bash

# This script will prepare an NCBI dataset so that it can be analysed using my
# fmhdist code.

OPTIND=1
N_THREADS=1
RANDOM_SEEDS=("42")
FMHDIST="fmhdist.jar"

function display_help() {
    echo "Assumes the current working directory contains folders each representing a genome."
    echo "Each genome folder should contain a fasta/fna file and a gtf file."
    echo "Usage"
    echo "prepare_analysis.sh -w W -k K -s S -f F [-r R] [-t T]"
    echo "W - comma separated list of window sizes"
    echo "K - comma separated list of k-mer sizes"
    echo "S - comma separated list of scaling parameters"
    echo "F - PATH to fmhdist executable"
    echo "R - comma separated list of random seeds, default: 42"
    echo "T - number of threads to use for the sketching"
}


while getopts "w:k:s:r:t:f:" opt; do
    case "$opt" in
    "w")
        IFS="," read -ra WINDOW_SIZES <<< "$OPTARG"
        ;;
    "k")
        IFS="," read -ra K_SIZES <<< "$OPTARG"
        ;;
    "s")
        IFS="," read -ra SCALING_PARAMS <<< "$OPTARG"
        ;;
    "r")
        IFS="," read -ra RANDOM_SEEDS <<< "$OPTARG"
        ;;
    "t")
        N_THREADS=$OPTARG
        ;;
    "f")
        FMHDIST=$OPTARG
        ;;
    *)
        display_help
        ;;
    esac
done

function check_tools() {
    if ! command -v splitfasta &> /dev/null; then
        echo "splitfasta not found in PATH, please install"
        return 1
    fi
    if ! command -v macle &> /dev/null; then
        echo "macle not found in PATH, please install"
        return 1
    fi
    if ! command -v java &> /dev/null; then
        echo "java not found in PATH, please install"
        return 1
    fi
    return 0
}

function prepare_input_files() {
    # Remove summary files that are typically included in an NCBI ZIP
    rm -f *.csv
    rm -f *.tsv
    rm -f *.json
    rm -f *.jsonl

    # Rename files and prepare folders
    for dir in *; do
        mv -f $dir/*.fna $dir/genome.fasta
        mv -f $dir/*.gtf $dir/genome.gtf
        mkdir $dir/macle
        mkdir $dir/coordinates
    done
}

function prepare_complexities() {
    for dir in *; do
        # splitfasta needs us to be in the correct directory..
        cd $dir
        splitfasta genome.fasta
        for file in genome_split_files/*; do
            name=$(head -n 1 $file | sed 's/>//g' | awk '{print $1}')
            mv $file genome_split_files/$name.fasta
            for w in ${WINDOW_SIZES[@]}; do
                macle -w $w genome_split_files/$name.fasta >> macle/$w.txt
            done;
        done
        cd ..
    done
}

function prepare_list() {
    for dir in *; do
        echo "$(pwd)/$dir/genome.fasta,$dir" >> list.csv
    done;
}

prepare_coordinates() {
    for s in ${SCALING_PARAMS[@]}; do
        for rs in ${RANDOM_SEEDS[@]}; do
            for k in ${K_SIZES[@]}; do
                echo "Sketching with k=$k, s=$s and random seed=$rs...";
                java -jar $FMHDIST sketch --input list.csv --output . -k $k -s $s -rs $rs -c -t $N_THREADS
                rm *.sketch
                for c in *.coordinates; do
                    mv $c $(basename -s .sketch.coordinates $c)/coordinates/k_${k}-s_${s}-rs_${rs}.coordinates
                done;
            done;
        done;
    done;
}

function main() {
    echo "assuming $(pwd) as base directory"
    echo "checking software requirements"
    if ! check_tools; then
        exit 1
    fi
    echo "using $FMHDIST"
    echo "using random seeds ${RANDOM_SEEDS[@]}"
    echo "using window sizes ${WINDOW_SIZES[@]}"
    echo "using k ${K_SIZES[@]}"
    echo "using s ${SCALING_PARAMS[@]}"
    echo "using $N_THREADS threads"

    prepare_input_files
    prepare_complexities
    prepare_list
    prepare_coordinates
}

main
