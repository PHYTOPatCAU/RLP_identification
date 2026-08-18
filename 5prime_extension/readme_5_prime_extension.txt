------------------------------# 5_extension


Include the files and scripts mentioned below in the same directory
1. Checking_gff_and_finding_right_protein_start.sh
2. First_gff_extension.py
3. Comparing_original_and_extended_proteins.py
4. exteracting-extended_justified_protein_txt_with_numberofbp_extension.sh
5. Final_justified_genes_uppdating_gff.py
6. original_protein.fasta (rename your input_protein.fasta as original_protein.fasta)
7. input_original.gff
8. protein_ids.txt (identical to headers in Input_protein.fasta): Make sure you have the same file name "protein_ids.txt"



####Install gffread and create a variable gffread to show the path to the gffread software
####Install HMM conda environment 
####Install 5_extension_env conda environment 
####Install blast conda environment 

############ All used conda environments are uploaded as .yml files, which can be used to create specific conda environments.
_________________________________________________
1. Checking_gff_and_finding_right_protein_start.sh

########This is the main script, and go through the script to add the path to your genome.fasta, input_gff, non-reductant proteins from a specific family like brassicaceae for Arabidopsis, and  create the gffread variable and add the path of gffread software.


_________________________________________________
2. In the file Comparing_original_and_extended_proteins.py


In line 57: Add the path of the same directory which you used to store all your input files and scripts  
57 line   base_dir = "~/IPS/genome.proteins/Athaliana/Correcting_the_gff_Ath/Extesnion_results"

           # Input files
60 line     original_fasta = os.path.join(base_dir, "original_protein.fasta")
           extended_fasta = os.path.join(base_dir, "Extracting_proteins_added_with_ext_name.fasta")
62 line    protein_ids_file = os.path.join(base_dir, "protein_ids.txt")
           output_fasta = os.path.join(base_dir, "Extracting_the_compared_and_justified_sequences.fasta")
           alignment_dir = os.path.join(base_dir, "alignments")
           csv_file = os.path.join(base_dir, "Extracting_extended_compared_justified_gff_updates.csv")

____________________________________________________
3. In the file exteracting-extended_justified_protein_txt_with_numberofbp_extension.sh

19th line: output="~/IPS/genome.proteins/Athaliana/Correcting_the_gff_Ath/Extesnion_results"
           rm -r ${output}/actual__extended_homlog_jusitified__genes.csv
           rm -r ${output}/Final_justified_genes_and_coordinateds.csv


## create the output variable by adding the path of your output directory; it is the same as the one used to store all your input files and scripts.



