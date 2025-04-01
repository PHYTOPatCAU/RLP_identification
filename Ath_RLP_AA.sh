#!/bin/bash

# Path variables ditected towards different directories, datasets, and softwares.
output='/work_beegfs/suaph302/output/Athout'
genome='/work_beegfs/suaph302/IPS/genome.proteins/Athaliana/Athaliana_Araport.primaryTransc.fasta'
c3c2f='/work_beegfs/suaph302/Pfam/NC2C3.hmm'
pfam='/work_beegfs/suaph302/Pfam'
signalp='/work_beegfs/suaph302/scripts/signalp-4.1/signalp'
phobius='/work_beegfs/suaph302/phobius/phobius.pl'
ips='/work_beegfs/suaph302/IPS/interproscan-5.65-97.0/interproscan.sh'
gff='/work_beegfs/suaph302/IPS/genome.proteins/Athaliana/Athaliana_447_Araport11.gene.gff3'
nu_genome='/work_beegfs/suaph302/IPS/genome.proteins/Athaliana/Athaliana_447_TAIR10.fa'
gffread='/work_beegfs/suaph302/gffread/gffread/gffread'




# Path variables directed towards files containing Pfam-33.0 domain IDs for NLRs, RLK and RLP resistant genes.
NLR=$(cat /work_beegfs/suaph302/Pfam/NLR_DM_list.txt)
LRR=$(cat /work_beegfs/suaph302/Pfam/Lrr_DM_list.txt)
RLK=$(cat /work_beegfs/suaph302/Pfam/Pkinase_DM_list.txt)

# path to scripts: needed for extracting results from server based DEEPTMHMM analysis.
scripts='/work_beegfs/suaph302/scripts/codes_RLP_AA'






# Flow of pipeline:  hmmer(c2-c3f) >> hmmer(to remove nlrs and rlks --> LRRs-NLRs-RLKs) >> ips >> deeptmhmm+phobius >> signalp >> relaxed + strict RLPs_set >> extracting nucleotide and protein fasta for 2 types of RLP set 
# At the end of the pipeline, we will have two sets of RLPs, relaxed and strict (complete and incomplete RLPs based on RLPs domains that they carry). 
module load gcc12-env
module load nano
module load miniconda3









# Using miniconda we created a HMM environemnt to include 'hmmer-3.1b2' software and relared dependencies. 
# With the help of stream line editior (sed) '*' sign from the protein sequences are removed. 
# sequences with asterisk '*' sign will not be analyzed by interproscan software.
# performed hmmsearch using genome protein data against c2-c3f-HMM/RLP-HMM database (E value of 1e-5) to filter proteins with C2-C3F domain.
# Again, 'sed' command is used to extract short result from hmmsearch output (untill the pattern 'hmmsearch' found in thge original output).
# Filtered proteins were extracted both as list of protein IDs and fasta sequences by using grep command and seqtk software.
# seqtk software is also embeded inside HMM environment with related dependencies.
# Later, removing the extra files that are not useful for future steps.


conda activate HMM
sed -i 's/\*//g' ${genome}
hmmsearch --tblout ${output}/hmmc2-c3F.out -E 1e-5 --cpu 10 ${c3c2f} ${genome} > ${output}/hmmc2-c3F.out
sed -n '/hmmsearch/q; p' ${output}/hmmc2-c3F.out > ${output}/hmmc2-c3F.sed.out
grep -P -o 'AT\d+\D\d+\.\d+' ${output}/hmmc2-c3F.sed.out | sort | uniq > ${output}/hmmc2-c3F.gene_list.txt
seqtk subseq ${genome} ${output}/hmmc2-c3F.gene_list.txt  > ${output}/hmmc2-c3F.gene_list.fasta

rm -r ${output}/hmmc2-c3F.out
rm -r ${output}/hmmc2-c3F.sed.out
rm -r ${output}/hmmc2-c3F.gene_list.txt








# Similar to RLP-HMM, database for filtering LRR are also generated.
# lrr.hmm profile and contected database are constructed using repective Pfam ids and profiles extracted from Pfam databases.
# 'IDs' means each LRR Pfam domain entries. There are 14 LRR entries in Pfam database. These entries are similar to Nagou et al. 2022 used in his paper concerted expansion and contraction of immune receptors in plant species.   
# I have added a small loop for extractig info for each LRR domain entries and building the LRR profile.
# Later, hmmpress command from hmmer-3.1b2 is used to convert the profiles into databases.


for IDs in ${LRR}; do
grep -B 2 -A 1300 "ACC[[:space:]]*${IDs}$" ${pfam}/Pfam-A.hmm | awk 'BEGIN{doPrint=1}{if (doPrint==1) {print $0}; if ($0=="//") {doPrint=0}; if ($0=="--") {doPrint=0}}' >> ${pfam}/lrr.hmm
done
hmmpress ${pfam}/lrr.hmm


# Next step in the pipeline is to extract the proteins with LRR domains as we focused on extracting resistant gene family having LRR domains.
# In this step, constructed LRR database is used to filter the proteins those initially passed C2-C3F domain prediction.
# Again, hmmersearch command is used to collect all the hit sequences against LRR database. hmmersearch command with an E value of 1e-3. The e value is taken from Nagou et al paper.
# Later, sed command is used to limit the result session. Further, grep and seqtk command is used in the same fashion as before to extract the protein sequence list and their respective fasta sequences.
# At the end of this step, we removed the files with no future importance.

hmmsearch --tblout ${output}/hmm.lrr.out -E 1e-3 --cpu 10 ${pfam}/lrr.hmm ${output}/hmmc2-c3F.gene_list.fasta > ${output}/hmm.lrr.out
sed -n '/hmmsearch/q; p' ${output}/hmm.lrr.out > ${output}/hmm.lrr.sed.out
grep -P -o 'AT\d+\D\d+\.\d+' ${output}/hmm.lrr.sed.out | sort | uniq > ${output}/hmm.lrr.genelist.txt
seqtk subseq ${genome} ${output}/hmm.lrr.genelist.txt > ${output}/hmm.lrr.genelist.fasta
rm -r ${output}/hmm.lrr.out
rm -r ${output}/hmm.lrr.sed.out








# Similar to LRR.HMM database, nlr.hmm database is constructed to extract all NOD-like receptors (they share LRR domains). But, the filtered proteins from this step are excluded from the main set of proteins. 
# nlr.hmm profile is developed using the Pfam domains unique to NLR genes but excluding LRR domains.  

for IDs in ${NLR}; do
grep -B 2 -A 1300 "ACC[[:space:]]*${IDs}$" ${pfam}/Pfam-A.hmm | awk 'BEGIN{doPrint=1}{if (doPrint==1) {print $0}; if ($0=="//") {doPrint=0}; if ($0=="--") {doPrint=0}}' >> ${pfam}/nlr.hmm
done
hmmpress ${pfam}/nlr.hmm



# Similar to lrr.hmm search, we used hmmscan for filtering nlr domain carrying protiens by passing via nlr.hmm database. The input protein set are the entires extracted in the previous lrr.filtering step.
# Here, hmmerscan command is used to collect all the hit sequences against nlr.hmm database. hmmerscan command with an E value of 1e-5.
# Later, sed command is used to limit the result session. Further, grep command is used in the same fashion as before to extract the protein sequence list.
# At the end of this step, we removed the files with no future importance.

hmmscan --tblout ${output}/hmm.nlr.out -E 1e-5 --cpu 15 ${pfam}/nlr.hmm ${output}/hmm.lrr.genelist.fasta > ${output}/hmm.nlr.out
sed -n '/hmmscan/q; p' ${output}/hmm.nlr.out > ${output}/hmm.nlr.sed.out
grep -P -o 'AT\d+\D\d+\.\d+' ${output}/hmm.nlr.sed.out | sort | uniq > ${output}/hmm.nlr-genelist.txt

rm -r ${output}/hmm.nlr.out
rm -r ${output}/hmm.nlr.sed.out








# Similar to LRR.HMM database, Pkinase.hmm database is constructed to extract all RLK receptors (they share LRR domains). But, the filtered proteins from this step are excluded from the main set of proteins.
# pkinase.hmm profile is developed using the Pfam domains unique to RLKs genes but excluding LRR domains.

for IDs in ${RLK}; do
grep -B 2 -A 1300 "ACC[[:space:]]*${IDs}$" ${pfam}/Pfam-A.hmm | awk 'BEGIN{doPrint=1}{if (doPrint==1) {print $0}; if ($0=="//") {doPrint=0}; if ($0=="--") {doPrint=0}}' >> ${pfam}/pkinase.hmm
done
hmmpress ${pfam}/pkinase.hmm


# Similar to nlr.hmm search, we used hmmscan for filtering piknase domain carrying protiens by passing via pkinase.hmm database. The input protein set are the entires extracted from lrr.filtering step.
# Here, hmmerscan command is used to collect all the hit sequences against pkinase.hmm database. hmmerscan command with an E value of 1e-5.
# Later, sed command is used to limit the result session. Further, grep command is used in the same fashion as before to extract the protein sequence list.
# At the end of this step, we removed the files with no future importance.

hmmscan --tblout ${output}/hmm.pkinase.out -E 1e-5 --cpu 10 ${pfam}/pkinase.hmm ${output}/hmm.lrr.genelist.fasta > ${output}/hmm.pkinase.out
sed -n '/hmmscan/q; p' ${output}/hmm.pkinase.out > ${output}/hmm.pkinase.sed.out
grep -P -o 'AT\d+\D\d+\.\d+' ${output}/hmm.pkinase.sed.out | sort | uniq > ${output}/hmm.pkinase.genelist.txt

rm -r ${output}/hmm.pkinase.sed.out
rm -r ${output}/hmm.pkinase.out





# Concatenating protein list having nlr and pkinase entires. Later, these proteins are deleted from the main set of LRR proteins to filter all probable rlps like proteins.
# cat, comm and seqtk commands are used to extract the rlps after deleting nlrs and rlks from the main set of lrr prteins.
# At the end, we have rlps that needs to be filtered for specific gene families via interproscan analyzis.


cat ${output}/hmm.pkinase.genelist.txt > ${output}/non-rlp.gene.list.txt
cat ${output}/hmm.nlr-genelist.txt   >> ${output}/non-rlp.gene.list.txt
comm -13 <(sort ${output}/non-rlp.gene.list.txt) <(sort ${output}/hmm.lrr.genelist.txt) > ${output}/rlpslist.txt
seqtk subseq ${genome} ${output}/rlpslist.txt > ${output}/rlpslist.fasta

rm -r ${output}/non-rlp.gene.list.txt
rm -r ${output}/hmm.pkinase.genelist.txt
rm -r ${output}/hmm.nlr-genelist.txt
rm -r ${output}/hmm.lrr.genelist.fasta
rm -r ${output}/hmm.lrr.genelist.txt






# Performing Interproscan and removing probable protein like kinases and others.
# IPS software is decoded into a new conda environment .i.e, annotation with all suppoorting dependencies.
# IPS is mainly directed to extract gene families like, 'SSF52058' and 'SSF52047' and remove all other protein families.
# At the end of the step, we remove files with no future help.
conda deactivate
conda activate annotation
${ips} -f TSV,GFF3 --appl SUPERFAMILY -i ${output}/rlpslist.fasta -b ${output}/rlpsips.out
grep -v -e 'SSF52047' -e 'SSF52058' ${output}/rlpsips.out.tsv | grep -P -o 'AT\d+\D\d+\.\d+' > ${output}/RMips.list.txt
comm -23 <(sort ${output}/rlpslist.txt) <(sort ${output}/RMips.list.txt) > ${output}/newrlp_ips_list.txt

conda deactivate
conda activate HMM
seqtk subseq ${genome} ${output}/newrlp_ips_list.txt > ${output}/newrlp_ips_list.fasta



rm -r ${output}/RMips.list.txt
rm -r ${output}/rlpsips.out.tsv
rm -r ${output}/rlpsips.out.gff3
rm -r ${output}/rlpslist.fasta
rm -r ${output}/rlpslist.txt






# Till here, the pipeline helped in extracting rlps (complete and incomplete proteins). Further, we segreated protiens into complete (SP+LRR+TM+SCT) and incomplete protiens (missing one SP/TM) using sofwares to identify Transmembrane and signal peptide domains.  

# running transmembrane analysis using DEEPTMHMM and PHOBIUS
# For both softwares, different conda environment has created to keep them separated with respective dependencies.
# Later, grep and cat commands are used to extract rlps with either SP or TM domains.

conda deactivate
conda activate deep

biolib run DTU/DeepTMHMM --fasta ${output}/newrlp_ips_list.fasta
grep 'SP+TM' ${scripts}/biolib_results/predicted_topologies.3line > ${output}/deep.TM.and.SP_TM.list.out
grep '| TM' ${scripts}/biolib_results/predicted_topologies.3line >> ${output}/deep.TM.and.SP_TM.list.out
cat ${output}/deep.TM.and.SP_TM.list.out | grep -P -o 'AT\d+\D\d+\.\d+' | sort | uniq > ${output}/deeptmhmm_TM.and.SP_TM.genelist.txt

conda deactivate
conda activate phobius

${phobius} -short ${output}/newrlp_ips_list.fasta > ${output}/phobius.out
cat ${output}/phobius.out | grep -v '0  Y' | grep -v '0  0' | grep -P -o 'AT\d+\D\d+\.\d+' | sort | uniq > ${output}/phobius.TM.and.SP_TM.gene.list.txt


# listing transmembrane having protein sequences
cat ${output}/deeptmhmm_TM.and.SP_TM.genelist.txt > ${output}/TM.and.SP_TM.gene.list.txt
cat ${output}/phobius.TM.and.SP_TM.gene.list.txt >> ${output}/TM.and.SP_TM.gene.list.txt
cat ${output}/TM.and.SP_TM.gene.list.txt | grep -P -o 'AT\d+\D\d+\.\d+' | sort | uniq > ${output}/Transmembrane.sort_uniq.genelist.txt


# listing 'Signal-peptide' having protein sequences
grep '| SP' ${scripts}/biolib_results/predicted_topologies.3line | grep -P -o 'AT\d+\D\d+\.\d+' > ${output}/Only.SP.list.out
grep 'SP+TM' ${scripts}/biolib_results/predicted_topologies.3line | grep -P -o 'AT\d+\D\d+\.\d+' >> ${output}/Only.SP.list.out
cat ${output}/phobius.out | grep '  Y' | grep -P -o 'AT\d+\D\d+\.\d+' >> ${output}/Only.SP.list.out




# running signalp-4.1 on the extracted rlp set after IPS analyzis.
# With a senstivity of 0.11, we extracted proteins with probable signal peptides.

conda deactivate
conda activate annotation
${signalp} -f short -u 0.11 -U 0.11 ${output}/newrlp_ips_list.fasta > ${output}/signalp.out 
grep 'Y' ${output}/signalp.out | grep -P -o 'AT\d+\D\d+\.\d+' | sort | uniq >> ${output}/Only.SP.list.out
cat ${output}/Only.SP.list.out | sort | uniq > ${output}/Signalp.list.txt



# extracting strict set of proteins
comm -12 <(sort ${output}/Transmembrane.sort_uniq.genelist.txt) <(sort ${output}/Signalp.list.txt) > ${output}/strict.list.txt


#extracting relaxed set of genes
cat ${output}/Transmembrane.sort_uniq.genelist.txt > ${output}/relaxed.list.out
cat ${output}/Signalp.list.txt >> ${output}/relaxed.list.out
cat ${output}/relaxed.list.out | sort | uniq > ${output}/relaxed.list.txt


rm -r ${output}/deep.TM.and.SP_TM.list.out
rm -r ${output}/deeptmhmm_TM.and.SP_TM.genelist.txt
rm -r ${output}/phobius.TM.and.SP_TM.gene.list.txt
rm -r ${output}/TM.and.SP_TM.gene.list.txt
rm -r ${output}/Transmembrane.sort_uniq.genelist.txt
rm -r ${output}/Only.SP.list.out
rm -r ${output}/signalp.out
rm -r ${output}/Signalp.list.txt
rm -r ${scripts}/biolib_results
rm -r ${output}/phobius.out
rm -r ${output}/relaxed.list.out




# In the next step we extract fasta sequences for the respective strict and relaxed set of RLPs.  
conda deactivate
conda activate HMM
seqtk subseq ${genome} ${output}/strict.list.txt > ${output}/strict.list.fasta
seqtk subseq ${genome} ${output}/relaxed.list.txt > ${output}/relaxed.list.fasta



conda deactivate
conda activate samtools

# extracting nucleotide sequence information newrlp_ips_list.txt 
awk '$3=="mRNA"' ${gff} > ${output}/mRNA.gff3
grep -wFf ${output}/newrlp_ips_list.txt ${output}/mRNA.gff3 > ${output}/newrlp_ips.mRNA.gff3
${gffread} -w ${output}/newrlp_ips.nu_genome.fasta -g ${nu_genome} -gff3 ${output}/newrlp_ips.mRNA.gff3

# extracting tm.sp_tm.signalp.genelist.txt sequences
awk '$3=="mRNA"' ${gff} > ${output}/mRNA.gff3
grep -wFf ${output}/strict.list.txt ${output}/mRNA.gff3 > ${output}/strict.mRNA.gff3
${gffread} -w ${output}/strict.genome.fasta -g ${nu_genome} -gff3 ${output}/strict.mRNA.gff3

# extracting sp.tm.sp_tm.signalp.genelist.txt sequences
awk '$3=="mRNA"' ${gff} > ${output}/mRNA.gff3
grep -wFf ${output}/relaxed.list.txt ${output}/mRNA.gff3 > ${output}/relaxed.mRNA.gff3
${gffread} -w ${output}/relaxed.genome.fasta -g ${nu_genome} -gff3 ${output}/relaxed.mRNA.gff3


# adding extra names to original rlpslist fasta with strict and relaxed geneids.
awk -v extra_name="relaxed" 'FNR==NR{gene_ids[$1]; next} /^>/ { gene_id=substr($1, 2); if(gene_id in gene_ids) print $1 " " extra_name; else print $0; next } 1' ${output}/relaxed.list.txt ${output}/newrlp_ips_list.fasta > ${output}/LRR_relaxed_set.fasta
awk -v new_extra_name="strict" 'FNR==NR{replace_gene_ids[$1]; next} /^>/ { gene_id=substr($1, 2); if(gene_id in replace_gene_ids) print $1 " " new_extra_name; else print $0; next } 1' ${output}/strict.list.txt ${output}/LRR_relaxed_set.fasta > ${output}/LRR_relaxed_srtict_set.fasta





