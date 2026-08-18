
**###### The PlantLRR-PRR pipeline is constructed to extract genome-wide LRR-RLP and LRR-RLK more precisely than any other previous pipeline.**

**###### 5prime\_extension is constructed to identify the correct start codon of all the annotated genes in your input data.**







\#### First, clone the complete git directory

**git clone https://github.com/PHYTOPatCAU/RLP\_identification.git**







\##### Install all the conda environments included in the "envs" directory

1. **HMM                   ##### used in PlantLRR-PRR pipeline** 
2. **annotation            ##### used in PlantLRR-PRR pipeline**
3. **5_extension_env       ##### used in 5prime-extension**







**##### locally install** 

1\. SignalP4.1

2\. Phobius-1.1

3\. Interproscan-5











\##### To successfully run the PlantLRR-PRR pipeline, write the correct paths to the following inside the original PlantLRR-PRR-pipeline.sh

output='/path/to/output'                                 # path to your output file

genome='/path/to/proteome.fasta'                         # path to your input proteome file

c3c2f='/path/to/Pfam/NC2C3.hmm'                          # path to your C2-C3F constructed HMMfile (Construction is included inside the original script )

pfam='/path/to/Pfam/'                                    # path to downloaded Pfam database.

signalp='/path/to/signalp-4.1/signalp'                   # path to your localy installed signalp4.1 software

phobius='/path/to/phobius/phobius.pl'                    # path to your localy installed phobius software

ips='/path/to/interproscan-5.65-97.0/interproscan.sh'    # path to your localy installed interproscan









\##### Since you have to run DEEPTMHMM and NetGPI on the web browser, extract "**LRR\_KINASE.fasta**" and **"newrlp\_ips\_list. fasta"** files after you run the PlantLRR-PRR pipeline for the first time

1. **${output}/LRR\_KINASE.fasta**
2. **${output}/newrlp\_ips\_list.fasta**







\#### Run both  "**LRR\_KINASE.fasta**" and **"newrlp\_ips\_list.fasta"** in DEEPTMHMM and extract the results in the form of  files and rename as,

1. RLP\_predicted\_topologies.3line      ### for  **"newrlp\_ips\_list.fasta"**   # output for RLPs
2. RLK\_predicted\_topologies.3line      ### for  "**LRR\_KINASE.fasta**"        # output for RLKs

\##### Store these two outputs inside the folder "**${output}/biolib\_results**"









\#### Again run both  "**LRR\_KINASE.fasta**" and **"newrlp\_ips\_list.fasta"** in NetGPI and extract the results in the form of  files and rename as,

1\. ${output}/output.gff3                   ### for  **"newrlp\_ips\_list.fasta"**      # output for RLPs

2\. ${output}/kinase.output.gff3            ### for  "**LRR\_KINASE.fasta**"           # output for RLKs





\##### At last, rerun the entire pipeline. The final complete RLPs and RLKs are stored under the filename,

1. strict.list.txt                        #### complete RLPs list
2. kinase.strict.list.txt                 #### complete RLKs list
3. strict.list.fasta                      #### complete RLPs fasta
4. kinase.strict.list.fasta               #### complete RLPs fasta







**####### Final \_output**



\# FINAL LRR-RLK outputs

1\. kinase.strict.list.txt                                          # output file containing all complete LRR-RLKs  

2\. kinase.relaxed.list.txt                                         # output file containing all incomplete LRR-RLKs 

3\. pkinase.genelist.txt                                            # output file containing all crude LRR-RLKs



\# FINAL LRR-RLK outputs

1\. strict.list.txt                                                 # output file containing all complete LRR-RLPs

2\. relaxed.list.txt                                                # output file containing all incomplete LRR-RLPs

3\. newrlp\_ips\_list.txt                                             # output file containing all crude LRR-RLPs





