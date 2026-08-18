#!/bin/bash

#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=220G
#SBATCH --time=20:00:00
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=base
#SBATCH --mail-user=n.rangegowda@phytomed.uni-kiel.de
#SBATCH --mail-user=nandee.jr@gmail.com
#SBATCH --mail-type=ALL

inputgff="~/IPS/genome.proteins/Athaliana/Correcting_the_gff_Ath/Extesnion_results/input.gff"
genome="~/IPS/genome.proteins/Athaliana/Athaliana_447_TAIR10.fa"
gffread="~/gffread/gffread/gffread"
nr="~/ncbi_nr/brassica/brassicaceae_nr"

output="~/IPS/genome.proteins/Athaliana/Correcting_the_gff_Ath/Extesnion_results"



module load gcc12-env/12.1.0
module load nano
module load miniconda3/4.12.0



rm -r ${output}/blastp_out_extended.txt
rm -r ${output}/All_compared_justified_extended_proteins.fasta
rm -r ${output}/All_compared_justified_extended_proteins.txt
rm -r ${output}/Extracting_proteins_added_with_ext_name.fasta
rm -r ${output}/Final_justified_genes_and_coordinateds.csv
rm -r ${output}/First_extended_protein_chaging_stop_codon.fasta
rm -r ${output}/First_extended_proteins_added_with_ext_name.fasta
rm -r ${output}/First_extendension.gff
rm -r ${output}/blastp_out_extended_justified_hits.txt
rm -r ${output}/final__jusitifed_genes.gff
rm -r ${output}/original__cds.fasta
rm -r ${output}/original__unmodified_header_cds.fasta
rm -r ${output}/original__unmodified_header_protein.fasta
rm -r ${output}/original_primary_protein.fasta
rm -r ${output}/whole__original__protein.fasta





conda activate 5_extension_env
python3 ${output}/First_gff_extension.py ${inputgff} ${output}/First_extendension.gff



# after extending the starting coordinated of all the genes, we will extract their cds and protiens sequences from respective gff and genome.fasta file using  gffread tool
${gffread} ${output}/First_extendension.gff -g ${genome} -x ${output}/First_extended_cds.fasta -y ${output}/First_extended_protein.fasta


# after extracting the extended proteins, using awk command, we will replace the "." with "*" in their sequences (no change in their headers) and then add extra names onto the sequence headers using sed command
awk '/^>/ {print; next} {gsub(/\./, "*"); print}' ${output}/First_extended_protein.fasta > ${output}/First_extended_protein_chaging_stop_codon.fasta
sed '/^>/s/$/\t\text/' ${output}/First_extended_protein_chaging_stop_codon.fasta | sed 's/\.Araport11\.447//'  > ${output}/First_extended_proteins_added_with_ext_name.fasta

conda deactivate
conda activate HMM 
seqtk subseq ${output}/First_extended_proteins_added_with_ext_name.fasta ${output}/protein_ids.txt  > ${output}/Extracting_proteins_added_with_ext_name.fasta


conda deactivate
conda activate 5_extension_env
# Later, we will compare all the extended sequences against original sequence to identify the proteins with extended N-terminal. Here, "clustal o" will help in multiply aligning two sequences and then try to find the "*" in the extended region and then if it finds then it will search for the "M" and then extracted all the amino acids till the end of the proteins. If it doesnot find any "*" then it will retian the original protien. Script we use for this process is "copilot_updating_gff.py"
python3 ${output}/Comparing_original_and_extended_proteins.py


cat ${output}/Extracting_extended_compared_justified_gff_updates.csv | tail -n +2 | cut -f 1 -d ',' | sort | uniq  > ${output}/All_compared_justified_extended_proteins.txt


conda deactivate 
conda activate HMM
seqtk subseq ${output}/Extracting_the_compared_and_justified_sequences.fasta ${output}/All_compared_justified_extended_proteins.txt > ${output}/All_compared_justified_extended_proteins.fasta


conda deactivate
conda activate blast
blastp -query ${output}/All_compared_justified_extended_proteins.fasta  -db ${nr} -out ${output}/blastp_out_extended.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen' -max_target_seqs 3  -evalue 1e-5 -num_threads 32

awk '($3 > 85) && (($8 - $7 + 1) >= 0.85 * $13)' ${output}/blastp_out_extended.txt  |  awk '$7 == $9 { print $1 }' | sort | uniq > ${output}/blastp_out_extended_justified_hits.txt 


#extreacting protein ids and their respective coordinated for all the compared and blastp_justified hits
${output}/exteracting-extended_justified_protein_txt_with_numberofbp_extension.sh 


conda activate 5_extension_env
python3 ${output}/Final_justified_genes_uppdating_gff.py ${inputgff} ${output}/Final_justified_genes_and_coordinateds.csv ${output}/final__jusitifed_genes.gff




${gffread} ${output}/final__jusitifed_genes.gff -g ${genome} -x ${output}/original__unmodified_header_cds.fasta -y ${output}/original__unmodified_header_protein.fasta

cat ${output}/original__unmodified_header_cds.fasta | sed 's/\.Araport11\.447//'  >  ${output}/original__cds.fasta
cat ${output}/original__unmodified_header_protein.fasta | sed 's/\.Araport11\.447//'  > ${output}/whole__original__protein.fasta


conda deactivate
conda activate HMM
seqtk subseq ${output}/whole__original__protein.fasta ${output}/protein_ids.txt > ${output}/original_primary_protein.fasta

conda deactivate

