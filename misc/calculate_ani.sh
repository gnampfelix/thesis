#!/bin/bash

sequences=(
    "/home/gnampfelix/Documents/term5/thesis/data/phytophthora/ref_db/queries/sequences/misc/GCA_022750515.1/GCA_022750515.1.fna"
    "/home/gnampfelix/Documents/term5/thesis/data/phytophthora/ref_db/queries/sequences/misc/GCA_008974285.1/GCA_008974285.1.fna"
    "/home/gnampfelix/Documents/term5/thesis/data/phytophthora/ref_db/queries/sequences/misc/GCA_014858625.1/GCA_014858625.1.fna"
    "/home/gnampfelix/Documents/term5/thesis/data/phytophthora/ref_db/queries/sequences/misc/GCA_018394375.1/GCA_018394375.1.fna"
    "/home/gnampfelix/Documents/term5/thesis/data/phytophthora/ref_db/queries/sequences/phy/GCA_001314365.1/GCA_001314365.1.fna"
    "/home/gnampfelix/Documents/term5/thesis/data/phytophthora/ref_db/sequences/GCF_000142945.1/GCF_000142945.1.fna"
)

blast_dir=/home/gnampfelix/.conda/envs/sourmash/bin/
ortho_ani_path=/home/gnampfelix/Downloads/OAT_cmd.jar
output=/tmp/dist.txt

for ((i=0;i<${#sequences[@]};i++)); do
    for ((j=i+1;j<${#sequences[@]};j++)); do
        s1=${sequences[i]}
        s2=${sequences[j]}
        echo "calculating pair $s1 $s2..."
        score=$(java -jar $ortho_ani_path -blastplus_dir $blast_dir -fasta1 $s1 -fasta2 $s2 | grep "OrthoANI" | awk '{print $3}')
        echo "$(basename -s .fna $s1) $(basename -s .fna $s2) $score" >> $output
    done
done