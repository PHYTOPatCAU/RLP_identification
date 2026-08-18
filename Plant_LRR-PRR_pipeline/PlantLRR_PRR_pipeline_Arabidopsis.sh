#!/bin/bash
#SBATCH --job-name=PlantLRR-PRR
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=PlantLRR-PRR.out
#SBATCH --error=PlantLRR-PRR.err
#SBATCH --partition=base
#SBATCH --mail-user=user_email_Id
#SBATCH --mail-type=ALL



##### Loading "gcc*" and "miniconda3" modules
module load gcc12-env
module load miniconda3



# Path variables ditected towards different directories, datasets, and software.
output='/path/to/output' # path to your output file
genome='/path/to/proteome.fasta' # path to your input proteome file
c3c2f='/path/to/Pfam/NC2C3.hmm' # path to your C2-C3F constructed HMMfile
pfam='/path/to/Pfam/' # path to downloaded Pfam database.
signalp='/path/to/signalp-4.1/signalp' # path to your localy installed signalp4.1 software
phobius='/path/to/phobius/phobius.pl' # path to your localy installed phobius software
ips='/path/to/interproscan-5.65-97.0/interproscan.sh' # path to your localy installed interproscan



# Path variables directed towards files containing Pfam-33.0 domain IDs for NLRs, RLK and RLP resistant genes.
NLR=$(cat ~/Pfam/NLR_DM_list.txt) # path to the file containing Pfam ids representing NLR domians
LRR=$(cat ~/Pfam/Lrr_DM_list.txt) # path to file containing Pfam ids representing LRR domians
RLK=$(cat ~/Pfam/Pkinase_DM_list.txt) # path to file containing Pfam ids representing RLK domians




### for loop to download and store Pfam database once and for all

for i in ${pfam}
        do
                if [ -f ${i}/Pfam-A.hmm ]
                        then
                                echo "Pfam 33.0 is downloaded"
                else
                                wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam33.0/Pfam-A.hmm.gz
                                mv Pfam-A.hmm.gz ${i}
                                cd ${i}
                                gunzip Pfam-A.hmm.gz
                fi
        done


conda activate HMM


for i in ${pfam}
        do
                if [ -f ${i}/NC2C3.hmm ]
                        then
                                echo "NC2C3F_HMM is already created"
                else
                                hmmbuild NC2C3.hmm ${i}/C2-C3F_domain.fasta
                fi
        done



# Flow of pipeline:  hmmer(c2-c3f) >> hmmer(to remove nlrs and rlks --> LRRs-NLRs-RLKs) >> ips >> deeptmhmm+phobius >> signalp >> relaxed + strict RLPs_set >> extracting nucleotide and protein fasta for 2 types of RLP set 
# At the end of the pipeline, we will have two sets of RLPs, relaxed and strict (complete and incomplete RLPs based on RLPs domains that they carry). 
module load gcc12-env
module load miniconda3



# Using miniconda we created a HMM environemnt to include 'hmmer-3.1b2' software and relared dependencies. 
# With the help of stream line editior (sed) '*' sign from the protein sequences are removed. 
# sequences with asterisk '*' sign will not be analyzed by interproscan software.
# performed hmmsearch using genome protein data against c2-c3f-HMM/RLP-HMM database (E value of 1e-3) to filter proteins with C2-C3F domain.
# Again, 'sed' command is used to extract short result from hmmsearch output (untill the pattern 'hmmsearch' found in thge original output).
# Filtered proteins were extracted both as list of protein IDs and fasta sequences by using grep command and seqtk software.
# seqtk software is also embeded inside HMM environment with related dependencies.
# Later, removing the extra files that are not useful for future steps.


conda activate HMM
sed -i 's/\*//g' ${genome}
hmmsearch --tblout ${output}/hmmc2-c3F.out -E 1e-3 --cpu 10 ${c3c2f} ${genome} > ${output}/hmmc2-c3F.out
sed -n '/hmmsearch/q; p' ${output}/hmmc2-c3F.out > ${output}/hmmc2-c3F.sed.out
grep -v '#' ${output}/hmmc2-c3F.sed.out | awk '{print $1}' | sort | uniq > ${output}/hmmc2-c3F.gene_list.txt
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
grep -v '#' ${output}/hmm.lrr.sed.out | awk '{print $1}' | sort | uniq > ${output}/hmm.lrr.genelist.txt
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
# Here, hmmerscan command is used to collect all the hit sequences against nlr.hmm database. hmmerscan command with an E value of 1e-3.
# Later, sed command is used to limit the result session. Further, grep command is used in the same fashion as before to extract the protein sequence list.
# At the end of this step, we removed the files with no future importance.

hmmscan --tblout ${output}/hmm.nlr.out -E 1e-3 --cpu 15 ${pfam}/nlr.hmm ${output}/hmm.lrr.genelist.fasta > ${output}/hmm.nlr.out
sed -n '/hmmscan/q; p' ${output}/hmm.nlr.out > ${output}/hmm.nlr.sed.out
grep -v '#' ${output}/hmm.nlr.sed.out | awk '{print $3}' | sort | uniq > ${output}/hmm.nlr-genelist.txt

rm -r ${output}/hmm.nlr.out
rm -r ${output}/hmm.nlr.sed.out








# Similar to LRR.HMM database, Pkinase.hmm database is constructed to extract all RLK receptors (they share LRR domains). But, the filtered proteins from this step are excluded from the main set of proteins.
# pkinase.hmm profile is developed using the Pfam domains unique to RLKs genes but excluding LRR domains.

for IDs in ${RLK}; do
grep -B 2 -A 1300 "ACC[[:space:]]*${IDs}$" ${pfam}/Pfam-A.hmm | awk 'BEGIN{doPrint=1}{if (doPrint==1) {print $0}; if ($0=="//") {doPrint=0}; if ($0=="--") {doPrint=0}}' >> ${pfam}/pkinase.hmm
done
hmmpress ${pfam}/pkinase.hmm


# Similar to nlr.hmm search, we used hmmscan for filtering piknase domain carrying protiens by passing via pkinase.hmm database. The input protein set are the entires extracted from lrr.filtering step.
# Here, hmmerscan command is used to collect all the hit sequences against pkinase.hmm database. hmmerscan command with an E value of 1e-3.
# Later, sed command is used to limit the result session. Further, grep command is used in the same fashion as before to extract the protein sequence list.
# At the end of this step, we removed the files with no future importance.

hmmscan --tblout ${output}/hmm.pkinase.out -E 1e-3 --cpu 10 ${pfam}/pkinase.hmm ${output}/hmm.lrr.genelist.fasta > ${output}/hmm.pkinase.out
sed -n '/hmmscan/q; p' ${output}/hmm.pkinase.out > ${output}/hmm.pkinase.sed.out
grep -v '#' ${output}/hmm.pkinase.sed.out | awk '{print $3}' | sort | uniq > ${output}/hmm.pkinase.genelist.txt

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
rm -r ${output}/hmm.nlr-genelist.txt
rm -r ${output}/hmm.lrr.genelist.fasta
rm -r ${output}/hmm.lrr.genelist.txt






# Performing Interproscan and removing probable protein like kinases and others.
# IPS software is decoded into a new conda environment .i.e, IPS_annotation with all suppoorting dependencies.
# IPS is mainly directed to extract gene families like, 'SSF52058' and 'SSF52047' and remove all other protein families.
# At the end of the step, we remove files with no future help.
conda deactivate
conda activate IPS_annotation
${ips} -f TSV,GFF3 --appl SUPERFAMILY,Gene3D -i ${output}/rlpslist.fasta -b ${output}/rlpsips.out
grep -e 'SSF56112' ${output}/rlpsips.out.tsv | awk '{print $1}' > ${output}/RMips.list.txt
grep -e 'G3DSA:1.10.510.10' ${output}/rlpsips.out.tsv | awk '{print $1}' >> ${output}/RMips.list.txt
comm -23 <(sort ${output}/rlpslist.txt) <(sort ${output}/RMips.list.txt) > ${output}/newrlp_ips_list.txt

conda deactivate
conda activate HMM
seqtk subseq ${genome} ${output}/newrlp_ips_list.txt > ${output}/newrlp_ips_list.fasta


# Extracting all the LRR-KInase with C2-C3f domain
grep 'SSF56112' ${output}/rlpsips.out.tsv | awk '{print $1}' >> ${output}/hmm.pkinase.genelist.txt
grep 'G3DSA:1.10.510.10' ${output}/rlpsips.out.tsv | awk '{print $1}' >> ${output}/hmm.pkinase.genelist.txt
cat ${output}/hmm.pkinase.genelist.txt | sort | uniq > ${output}/pkinase.genelist.txt
seqtk subseq ${genome} ${output}/pkinase.genelist.txt > ${output}/LRR_KINASE.fasta



######
######
###### After running the pipeline till here, extract the complete proteins from "newrlp_ips_list.fasta" and "LRR_KINASE.fasta" and run DeepTMHMM and NetGHPI web-based softwares.
###### Extract protein.topologies.3line results from the DeepTMHMM results and rename RLP and RLK results into RLK_predicted_topologies.3line and RLP_predicted_topologies.3line and store them under the folder "biolib_results" in "${output}" folder.
###### Extract the results from NetGPI and rename RLP and RLK outputs into RLP_netgpi_output.gff3 and RLK_netgpi_output.gff3, respectively and store these file in the "${output}" folder

### you can rerun the whole pipeline after adding the results and results will be generated in the end.

cat ${output}/rlpsips.out.tsv > ${output}/interproscan.output.txt


rm -r ${output}/RMips.list.txt
rm -r ${output}/rlpslist.fasta
rm -r ${output}/rlpslist.txt
rm -r ${output}/rlpsips.out.tsv
rm -r ${output}/rlpsips.out.gff3




# producing fasta file by removing first 70 sequences fromt each sequences
sed '/^>/!s/^.\{70\}//' ${output}/LRR_KINASE.fasta > ${output}/phobius_LRR_kinase_no-50AA.fasta
sed '/^>/!s/^.\{70\}//' ${output}/newrlp_ips_list.fasta > ${output}/phobius_LRR_RLP_no-50AA.fasta





# Till here, the pipeline helped in extracting rlps (complete and incomplete proteins). Further, we segreated protiens into complete (SP+LRR+TM+SCT) and incomplete protiens (missing one SP/TM) using sofwares to identify Transmembrane and signal peptide domains.

# running transmembrane analysis using DEEPTMHMM, NetGPI and PHOBIUS
# For both softwares, different conda environment has created to keep them separated with respective dependencies.
# Later, grep and cat commands are used to extract rlps with either SP or TM domains.
###### After running the pipeline till here, extract the complete proteins from "newrlp_ips_list.fasta" and "LRR_KINASE.fasta" and run DeepTMHMM and NetGHPI web-based softwares.
###### Extract protein.topologies.3line results from the DeepTMHMM results and rename RLP and RLK results into RLK_predicted_topologies.3line and RLP_predicted_topologies.3line and store them under the folder "biolib_results" in "${output}" folder.
###### Extract the results from NetGPI and rename RLP and RLK outputs into RLP_netgpi_output.gff3 and RLK_netgpi_output.gff3, respectively and store these file in the "${output}" folder

grep 'SP+TM' ${output}/biolib_results/RLP_predicted_topologies.3line > ${output}/deep.TM.and.SP_TM.list.out
grep '| TM' ${output}/biolib_results/RLP_predicted_topologies.3line >> ${output}/deep.TM.and.SP_TM.list.out
cat ${output}/deep.TM.and.SP_TM.list.out | awk '{print $1}' | sed 's/>//g' | sort | uniq > ${output}/deeptmhmm_TM.and.SP_TM.genelist.txt



${phobius} -short ${output}/phobius_LRR_RLP_no-50AA.fasta > ${output}/phobius.out
cat ${output}/phobius.out | grep -v '0  Y' | grep -v '0  0' | tail -n +2 | awk '{print $1}' | sort | uniq > ${output}/phobius.TM.and.SP_TM.gene.list.txt


# listing transmembrane having protein sequences
cat ${output}/deeptmhmm_TM.and.SP_TM.genelist.txt > ${output}/TM.and.SP_TM.gene.list.txt
cat ${output}/phobius.TM.and.SP_TM.gene.list.txt >> ${output}/TM.and.SP_TM.gene.list.txt
cat ${output}/RLP_netgpi_output.gff3 | grep -v '#' | awk '{print $1}' >> ${output}/TM.and.SP_TM.gene.list.txt
cat ${output}/TM.and.SP_TM.gene.list.txt | sort | uniq > ${output}/Transmembrane.sort_uniq.genelist.txt


# listing 'Signal-peptide' having protein sequences
grep '| SP' ${output}/biolib_results/RLP_predicted_topologies.3line | awk '{print $1}' | sed 's/>//g' > ${output}/Only.SP.list.out
grep 'SP+TM' ${output}/biolib_results/RLP_predicted_topologies.3line | awk '{print $1}' | sed 's/>//g' >> ${output}/Only.SP.list.out




# running signalp-4.1 on the extracted rlp set after IPS analyzis.
# With a senstivity of 0.11, we extracted proteins with probable signal peptides.


${signalp} -f short -u 0.11 -U 0.11 ${output}/newrlp_ips_list.fasta > ${output}/signalp.out
grep 'Y' ${output}/signalp.out | grep -v "#" | awk '{print $1}' | sort | uniq >> ${output}/Only.SP.list.out
cat ${output}/Only.SP.list.out | sort | uniq > ${output}/Signalp.list.txt



# extracting strict set of proteins
comm -12 <(sort ${output}/Transmembrane.sort_uniq.genelist.txt) <(sort ${output}/Signalp.list.txt) > ${output}/strict.list.txt


#extracting relaxed set of genes
cat ${output}/Transmembrane.sort_uniq.genelist.txt > ${output}/relaxed.list.out
cat ${output}/Signalp.list.txt >> ${output}/relaxed.list.out
cat ${output}/relaxed.list.out | sort | uniq > ${output}/relaxed.list.txt

#Extracting for LRR_Kinases
# Till here, the pipeline helped in extracting rlps (complete and incomplete proteins). Further, we segreated protiens into complete (SP+LRR+TM+SCT) and incomplete protiens (missing one SP/TM) using sofwares to identify Transmembrane an>
# running transmembrane analysis using DEEPTMHMM, NetGPI and PHOBIUS
# For both softwares, different conda environment has created to keep them separated with respective dependencies.
# Later, grep and cat commands are used to extract rlps with either SP or TM domains.
###### After running the pipeline till here, extract the complete proteins from "newrlp_ips_list.fasta" and "LRR_KINASE.fasta" and run DeepTMHMM and NetGHPI web-based softwares.
###### Extract protein.topologies.3line results from the DeepTMHMM results and rename RLP and RLK results into RLK_predicted_topologies.3line and RLP_predicted_topologies.3line and store them under the folder "biolib_results" in "${output}" folder.
###### Extract the results from NetGPI and rename RLP and RLK outputs into RLP_netgpi_output.gff3 and RLK_netgpi_output.gff3, respectively and store these file in the "${output}" folder


grep 'SP+TM' ${output}/biolib_results/RLK_predicted_topologies.3line > ${output}/kinase.deep.TM.and.SP_TM.list.out
grep '| TM' ${output}/biolib_results/RLK_predicted_topologies.3line >> ${output}/kinase.deep.TM.and.SP_TM.list.out
cat ${output}/kinase.deep.TM.and.SP_TM.list.out | awk '{print $1}' | sed 's/>//g' | sort | uniq > ${output}/kinase.deeptmhmm_TM.and.SP_TM.genelist.txt



${phobius} -short ${output}/phobius_LRR_kinase_no-50AA.fasta > ${output}/kinase.phobius.out
cat ${output}/kinase.phobius.out | grep -v '0  Y' | grep -v '0  0' | tail -n +2 | awk '{print $1}' | sort | uniq > ${output}/kinase.phobius.TM.and.SP_TM.gene.list.txt


# listing transmembrane having protein sequences
cat ${output}/kinase.deeptmhmm_TM.and.SP_TM.genelist.txt > ${output}/kinase.TM.and.SP_TM.gene.list.txt
cat ${output}/kinase.phobius.TM.and.SP_TM.gene.list.txt >> ${output}/kinase.TM.and.SP_TM.gene.list.txt
cat ${output}/RLK_netgpi_output.gff3 | grep -v '#' | awk '{print $1}' >> ${output}/kinase.TM.and.SP_TM.gene.list.txt
cat ${output}/kinase.TM.and.SP_TM.gene.list.txt | sort | uniq > ${output}/kinase.Transmembrane.sort_uniq.genelist.txt



# listing 'Signal-peptide' having protein sequences
grep '| SP' ${output}/biolib_results/RLK_predicted_topologies.3line | awk '{print $1}' | sed 's/>//g' > ${output}/kinase.Only.SP.list.out
grep 'SP+TM' ${output}/biolib_results/RLK_predicted_topologies.3line | awk '{print $1}' | sed 's/>//g' >> ${output}/kinase.Only.SP.list.out



# running signalp-4.1 on the extracted rlp set after IPS analyzis.
# With a senstivity of 0.11, we extracted proteins with probable signal peptides.


${signalp} -f short -u 0.11 -U 0.11 ${output}/LRR_KINASE.fasta > ${output}/kinase.signalp.out
grep 'Y' ${output}/kinase.signalp.out | grep -v "#" | awk '{print $1}' | sort | uniq >> ${output}/kinase.Only.SP.list.out
cat ${output}/kinase.Only.SP.list.out | sort | uniq > ${output}/kinase.Signalp.list.txt



# extracting strict set of proteins
comm -12 <(sort ${output}/kinase.Transmembrane.sort_uniq.genelist.txt) <(sort ${output}/kinase.Signalp.list.txt) > ${output}/kinase.strict.list.txt


#extracting relaxed set of genes
cat ${output}/kinase.Transmembrane.sort_uniq.genelist.txt > ${output}/kinase.relaxed.list.out
cat ${output}/kinase.Signalp.list.txt >> ${output}/kinase.relaxed.list.out
cat ${output}/kinase.relaxed.list.out | sort | uniq > ${output}/kinase.relaxed.list.txt







rm -r ${output}/deep.TM.and.SP_TM.list.out
rm -r ${output}/deeptmhmm_TM.and.SP_TM.genelist.txt
rm -r ${output}/phobius.TM.and.SP_TM.gene.list.txt
rm -r ${output}/TM.and.SP_TM.gene.list.txt
rm -r ${output}/Transmembrane.sort_uniq.genelist.txt
rm -r ${output}/Only.SP.list.out
rm -r ${output}/signalp.out
rm -r ${output}/Signalp.list.txt
rm -r ${output}/phobius.out
rm -r ${output}/relaxed.list.out




