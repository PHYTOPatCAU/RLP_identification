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
#SBATCH --mail-user=your@email.id
#SBATCH --mail-type=ALL


# Path variables ditected towards different directories, datasets, and softwares.

output='/work_beegfs/suaph302/output/solpeni_0716/TESTING'                                                             # path to your output file

genome='/work_beegfs/suaph302/solanum_genome/solpeni_0716/primarysopeni0716.protein.fasta'                     # path to your genome.fasta file

c3c2f='/work_beegfs/suaph302/Pfam/NC2C3.hmm'                                                                   # path to c2_c3F_hmm

pfam='/work_beegfs/suaph302/Pfam'                                                                              # Main path to C2_C3F, Pkinase, LRR, and NLR HMMs

signalp='/work_beegfs/suaph302/scripts/signalp-4.1/signalp'                                                    # path to SignalP4.1, for signal-peptide identification

phobius='/work_beegfs/suaph302/phobius/phobius.pl'                                                            # path to Phobius1.1, for transmembrane identification

ips='/work_beegfs/suaph302/IPS/interproscan-5.65-97.0/interproscan.sh'                                        # path to interproscan-5, to identify discountinous kinase domains



# this is the geneId we used to extract all the protein hits from the output of all the analysis.
gene_ID='Sopen\d+\D+\d+\.\d+'



# Flow of pipeline:  hmmer(c2-c3f) >> hmmer(to remove nlrs and rlks --> LRRs-NLRs-RLKs) >> ips >> deeptmhmm+phobius >> signalp >> relaxed + strict RLPs_set >> extracting nucleotide and protein fasta for 2 types of RLP set
# At the end of the pipeline, we will have two sets of RLPs, relaxed and strict (complete and incomplete RLPs based on RLPs domains that they carry).
module load gcc12-env/12.1.0
module load nano
module load miniconda3/4.12.0





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
hmmsearch --tblout ${output}/hmmc2-c3F.out -E 1e-3 --cpu 10 ${c3c2f} ${genome} > ${output}/hmmc2-c3F.out
sed -n '/hmmsearch/q; p' ${output}/hmmc2-c3F.out > ${output}/hmmc2-c3F.sed.out
grep -P -o ${gene_ID} ${output}/hmmc2-c3F.sed.out | sort | uniq > ${output}/hmmc2-c3F.gene_list.txt
seqtk subseq ${genome} ${output}/hmmc2-c3F.gene_list.txt  > ${output}/hmmc2-c3F.gene_list.fasta



# Next step in the pipeline is to extract the proteins with LRR domains as we focused on extracting resistant gene family having LRR domains.
# In this step, constructed LRR database is used to filter the proteins those initially passed C2-C3F domain prediction.
# Again, hmmersearch command is used to collect all the hit sequences against LRR database. hmmersearch command with an E value of 1e-3. The e value is taken from Nagou et al paper.
# Later, sed command is used to limit the result session. Further, grep and seqtk command is used in the same fashion as before to extract the protein sequence list and their respective fasta sequences.
# At the end of this step, we removed the files with no future importance.

hmmsearch --tblout ${output}/hmm.lrr.out -E 1e-3 --cpu 10 ${pfam}/lrr.hmm ${output}/hmmc2-c3F.gene_list.fasta > ${output}/hmm.lrr.out
sed -n '/hmmsearch/q; p' ${output}/hmm.lrr.out > ${output}/hmm.lrr.sed.out
grep -P -o ${gene_ID} ${output}/hmm.lrr.sed.out | sort | uniq > ${output}/hmm.lrr.genelist.txt
seqtk subseq ${genome} ${output}/hmm.lrr.genelist.txt > ${output}/hmm.lrr.genelist.fasta
rm -r ${output}/hmm.lrr.out
rm -r ${output}/hmm.lrr.sed.out




# Similar to nlr.hmm search, we used hmmscan for filtering piknase domain carrying protiens by passing via pkinase.hmm database. The input protein set are the entires extracted from lrr.filtering step.
# Here, hmmerscan command is used to collect all the hit sequences against pkinase.hmm database. hmmerscan command with an E value of 1e-5.
# Later, sed command is used to limit the result session. Further, grep command is used in the same fashion as before to extract the protein sequence list.
# At the end of this step, we removed the files with no future importance.

hmmscan --tblout ${output}/hmm.pkinase.out -E 1e-3 --cpu 10 ${pfam}/pkinase.hmm ${output}/hmm.lrr.genelist.fasta > ${output}/hmm.pkinase.out
sed -n '/hmmscan/q; p' ${output}/hmm.pkinase.out > ${output}/hmm.pkinase.sed.out
grep -P -o ${gene_ID} ${output}/hmm.pkinase.sed.out | sort | uniq > ${output}/hmm.pkinase.genelist.txt

rm -r ${output}/hmm.pkinase.sed.out
rm -r ${output}/hmm.pkinase.out



# Similar to lrr.hmm search, we used hmmscan for filtering nlr domain carrying protiens by passing via nlr.hmm database. The input protein set are the entires extracted in the previous lrr.filtering step.
# Here, hmmerscan command is used to collect all the hit sequences against nlr.hmm database. hmmerscan command with an E value of 1e-5.
# Later, sed command is used to limit the result session. Further, grep command is used in the same fashion as before to extract the protein sequence list.
# At the end of this step, we removed the files with no future importance.

hmmscan --tblout ${output}/hmm.nlr.out -E 1e-3 --cpu 15 ${pfam}/nlr.hmm ${output}/hmm.lrr.genelist.fasta > ${output}/hmm.nlr.out
sed -n '/hmmscan/q; p' ${output}/hmm.nlr.out > ${output}/hmm.nlr.sed.out
grep -P -o ${gene_ID} ${output}/hmm.nlr.sed.out | sort | uniq > ${output}/hmm.nlr-genelist.txt

rm -r ${output}/hmm.nlr.out
rm -r ${output}/hmm.nlr.sed.out




# Concatenating protein list having nlr and pkinase entires. Later, these proteins are deleted from the main set of LRR proteins to filter all probable rlps like proteins.
# cat, comm and seqtk commands are used to extract the rlps after deleting nlrs and rlks from the main set of lrr prteins.
# At the end, we have rlps that needs to be filtered for specific gene families via interproscan analyzis.


cat ${output}/hmm.pkinase.genelist.txt > ${output}/non-rlp.gene.list.txt
cat ${output}/hmm.nlr-genelist.txt   >> ${output}/non-rlp.gene.list.txt
comm -13 <(sort ${output}/non-rlp.gene.list.txt) <(sort ${output}/hmm.lrr.genelist.txt) > ${output}/rlpslist.txt
seqtk subseq ${genome} ${output}/rlpslist.txt > ${output}/rlpslist.fasta

rm -r ${output}/non-rlp.gene.list.txt
rm -r ${output}/hmm.lrr.genelist.fasta
rm -r ${output}/hmm.lrr.genelist.txt





# Performing Interproscan and removing probable protein like kinases and others.
# IPS software is decoded into a new conda environment .i.e, annotation with all suppoorting dependencies.
# IPS is mainly directed to extract gene families like, 'SSF52058' and 'SSF52047' and remove all other protein families.
# At the end of the step, we remove files with no future help.
#conda deactivate
#conda activate annotation
${ips} -f TSV,GFF3 --appl SUPERFAMILY,Gene3D -i ${output}/rlpslist.fasta -b ${output}/rlpsips.out
grep -e 'SSF56112' ${output}/rlpsips.out.tsv | grep -P -o ${gene_ID} > ${output}/RMips.list
grep -e 'G3DSA:1.10.510.10' ${output}/rlpsips.out.tsv | grep -P -o ${gene_ID} >> ${output}/RMips.list
cat ${output}/RMips.list | sort | uniq > ${output}/RMips.list.txt
comm -23 <(sort ${output}/rlpslist.txt) <(sort ${output}/RMips.list.txt) > ${output}/newrlp_ips_list.txt

conda deactivate
conda activate HMM
seqtk subseq ${genome} ${output}/newrlp_ips_list.txt > ${output}/newrlp_ips_list.fasta


# Extracting all the LRR-KInase with C2-C3f domain
grep 'SSF56112' ${output}/rlpsips.out.tsv | grep -P -o ${gene_ID} >> ${output}/hmm.pkinase.genelist.txt
grep 'G3DSA:1.10.510.10' ${output}/rlpsips.out.tsv | grep -P -o ${gene_ID} >> ${output}/hmm.pkinase.genelist.txt
grep -P -o ${gene_ID} ${output}/hmm.pkinase.genelist.txt | sort | uniq > ${output}/pkinase.genelist.txt
seqtk subseq ${genome} ${output}/pkinase.genelist.txt > ${output}/LRR_KINASE.fasta


cat ${output}/rlpsips.out.tsv > ${output}/SSF56112_results.tsv
rm -r ${output}/RMips.list.txt
rm -r ${output}/rlpslist.fasta
rm -r ${output}/rlpslist.txt
rm -r ${output}/rlpsips.out*



# producing fasta file by removing first 50 sequences fromt each sequences
sed '/^>/!s/^.\{70\}//' ${output}/LRR_KINASE.fasta > ${output}/phobius_LRR_kinase_no-50AA.fasta
sed '/^>/!s/^.\{70\}//' ${output}/newrlp_ips_list.fasta > ${output}/phobius_LRR_RLP_no-50AA.fasta









# Till here, the pipeline helped in extracting rlps (complete and incomplete proteins). Further, we segreated protiens into complete (SP+LRR+TM+SCT) and incomplete protiens (missing one SP/TM) using sofwares to identify Transmembrane and signal peptide domains.

# running transmembrane analysis using DEEPTMHMM and PHOBIUS
# For both softwares, different conda environment has created to keep them separated with respective dependencies.
# Later, grep and cat commands are used to extract rlps with either SP or TM domains.

conda deactivate
conda activate deep

#biolib run DTU/DeepTMHMM --fasta ${output}/newrlp_ips_list.fasta
grep 'SP+TM' ${scripts}/biolib_results/RLP_predicted_topologies.3line > ${output}/deep.TM.and.SP_TM.list.out
grep '| TM' ${scripts}/biolib_results/RLP_predicted_topologies.3line >> ${output}/deep.TM.and.SP_TM.list.out
cat ${output}/deep.TM.and.SP_TM.list.out | grep -P -o ${gene_ID} | sort | uniq > ${output}/deeptmhmm_TM.and.SP_TM.genelist.txt


${phobius} -short ${output}/phobius_LRR_RLP_no-50AA.fasta > ${output}/phobius.out
cat ${output}/phobius.out | grep -v '0  Y' | grep -v '0  0' | grep -P -o ${gene_ID} | sort | uniq > ${output}/phobius.TM.and.SP_TM.gene.list.txt


# listing transmembrane having protein sequences
cat ${output}/deeptmhmm_TM.and.SP_TM.genelist.txt > ${output}/TM.and.SP_TM.gene.list.txt
cat ${output}/phobius.TM.and.SP_TM.gene.list.txt >> ${output}/TM.and.SP_TM.gene.list.txt
cat ${output}/output.gff3 | grep -P -o ${gene_ID} >> ${output}/TM.and.SP_TM.gene.list.txt
cat ${output}/TM.and.SP_TM.gene.list.txt | grep -P -o ${gene_ID} | sort | uniq > ${output}/Transmembrane.sort_uniq.genelist.txt



# listing 'Signal-peptide' having protein sequences
grep '| SP' ${scripts}/biolib_results/RLP_predicted_topologies.3line | grep -P -o ${gene_ID} > ${output}/Only.SP.list.out
grep 'SP+TM' ${scripts}/biolib_results/RLP_predicted_topologies.3line | grep -P -o ${gene_ID} >> ${output}/Only.SP.list.out
cat ${output}/phobius.out | grep '  Y' | grep -P -o ${gene_ID} >> ${output}/Only.SP.list.out




# running signalp-4.1 on the extracted rlp set after IPS analyzis.
# With a senstivity of 0.11, we extracted proteins with probable signal peptides.

conda deactivate
conda activate annotation
${signalp} -f short -u 0.11 -U 0.11 ${output}/newrlp_ips_list.fasta > ${output}/signalp.out
grep 'Y' ${output}/signalp.out | grep -P -o ${gene_ID} | sort | uniq >> ${output}/Only.SP.list.out
cat ${output}/Only.SP.list.out | sort | uniq > ${output}/Signalp.list.txt



# extracting strict set of proteins
comm -12 <(sort ${output}/Transmembrane.sort_uniq.genelist.txt) <(sort ${output}/Signalp.list.txt) > ${output}/strict.list.txt


#extracting relaxed set of genes
cat ${output}/Transmembrane.sort_uniq.genelist.txt > ${output}/relaxed.list.out
cat ${output}/Signalp.list.txt >> ${output}/relaxed.list.out
cat ${output}/relaxed.list.out | sort | uniq > ${output}/relaxed.list.txt






#Extracting for LRR_Kinases
# Till here, the pipeline helped in extracting rlps (complete and incomplete proteins). Further, we segreated protiens into complete (SP+LRR+TM+SCT) and incomplete protiens (missing one SP/TM) using sofwares to identify Transmembrane an>
# running transmembrane analysis using DEEPTMHMM and PHOBIUS
# For both softwares, different conda environment has created to keep them separated with respective dependencies.
# Later, grep and cat commands are used to extract rlps with either SP or TM domains.

conda deactivate
conda activate deep

#biolib run DTU/DeepTMHMM --fasta ${output}/LRR_KINASE.fasta
grep 'SP+TM' ${scripts}/biolib_results/RLK_predicted_topologies.3line > ${output}/kinase.deep.TM.and.SP_TM.list.out
grep '| TM' ${scripts}/biolib_results/RLK_predicted_topologies.3line >> ${output}/kinase.deep.TM.and.SP_TM.list.out
cat ${output}/kinase.deep.TM.and.SP_TM.list.out | grep -P -o ${gene_ID} | sort | uniq > ${output}/kinase.deeptmhmm_TM.and.SP_TM.genelist.txt


${phobius} -short ${output}/phobius_LRR_kinase_no-50AA.fasta > ${output}/kinase.phobius.out
cat ${output}/kinase.phobius.out | grep -v '0  Y' | grep -v '0  0' | grep -P -o ${gene_ID} | sort | uniq > ${output}/kinase.phobius.TM.and.SP_TM.gene.list.txt


# listing transmembrane having protein sequences
cat ${output}/kinase.deeptmhmm_TM.and.SP_TM.genelist.txt > ${output}/kinase.TM.and.SP_TM.gene.list.txt
cat ${output}/kinase.phobius.TM.and.SP_TM.gene.list.txt >> ${output}/kinase.TM.and.SP_TM.gene.list.txt
cat ${output}/kinase.output.gff3 | grep -P -o ${gene_ID} >> ${output}/kinase.TM.and.SP_TM.gene.list.txt
cat ${output}/kinase.TM.and.SP_TM.gene.list.txt | grep -P -o ${gene_ID} | sort | uniq > ${output}/kinase.Transmembrane.sort_uniq.genelist.txt



# listing 'Signal-peptide' having protein sequences
grep '| SP' ${scripts}/biolib_results/RLK_predicted_topologies.3line | grep -P -o ${gene_ID} > ${output}/kinase.Only.SP.list.out

grep 'SP+TM' ${scripts}/biolib_results/RLK_predicted_topologies.3line | grep -P -o ${gene_ID} >> ${output}/kinase.Only.SP.list.out
cat ${output}/kinase.phobius.out | grep '  Y' | grep -P -o ${gene_ID} >> ${output}/kinase.Only.SP.list.out




# running signalp-4.1 on the extracted rlp set after IPS analyzis.
# With a senstivity of 0.11, we extracted proteins with probable signal peptides.

conda deactivate
conda activate annotation
${signalp} -f short -u 0.11 -U 0.11 ${output}/LRR_KINASE.fasta > ${output}/kinase.signalp.out
grep 'Y' ${output}/kinase.signalp.out | grep -P -o ${gene_ID} | sort | uniq >> ${output}/kinase.Only.SP.list.out
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







