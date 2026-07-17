#!/bin/bash

#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=220G
#SBATCH --time=10:00:00
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=base
#SBATCH --mail-user=n.rangegowda@phytomed.uni-kiel.de
#SBATCH --mail-type=ALL

module load gcc12-env/12.1.0
module load nano
module load miniconda3/4.12.0

output="~/IPS/genome.proteins/Athaliana/Correcting_the_gff_Ath/Extesnion_results"
rm -r ${output}/actual__extended_homlog_jusitified__genes.csv
rm -r ${output}/Final_justified_genes_and_coordinateds.csv

uniq_txt=$(cat ${output}/blastp_out_extended_justified_hits.txt)
gff_csv=$(cat Extracting_extended_compared_justified_gff_updates.csv)


for i in ${uniq_txt}
do
    grep "$i" "$gff_csv" >> ${output}/actual__extended_homlog_jusitified__genes.csv
done

cat ${gff_csv} | cut -f 1,3 -d ',' | head -n 1 > ${output}/Final_justified_genes_and_coordinateds.csv
cat ${output}/actual__extended_homlog_jusitified__genes.csv | cut -f 1,3 -d ',' | tail -n +1 >> ${output}/Final_justified_genes_and_coordinateds.csv 




