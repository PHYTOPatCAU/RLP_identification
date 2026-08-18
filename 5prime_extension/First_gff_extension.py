import sys

def extend_coordinates(gff_file, output_file):
    genes = {}  # Store gene information to track first/last exons and CDS
    chromosome_lengths = {}  # Optional: Store chromosome lengths if available
    updated_lines = 0

    with open(gff_file, 'r') as infile:
        lines = infile.readlines()

    # First pass: Identify first exon/CDS for + strand and last exon/CDS for - strand
    for line in lines:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue

        chrom, _, feature_type, start, end, _, strand, _, attributes = fields
        start, end = int(start), int(end)
        
        # Extract gene ID
        gene_id = None
        for attr in attributes.split(";"):
            if attr.startswith("ID="):
                gene_id = attr.split("=")[1]
            elif attr.startswith("Parent="):
                gene_id = attr.split("=")[1]  # For exons/CDSs, use Parent ID
        if not gene_id:
            continue

        if gene_id not in genes:
            genes[gene_id] = {"strand": strand, "chrom": chrom, "gene": None, "mRNA": None, "first_exon": None, "last_exon": None, "first_CDS": None, "last_CDS": None}

        if feature_type == "gene":
            genes[gene_id]["gene"] = (start, end)
        elif feature_type == "mRNA":
            genes[gene_id]["mRNA"] = (start, end)
        elif feature_type == "exon":
            if strand == "+":
                if not genes[gene_id]["first_exon"] or start < genes[gene_id]["first_exon"][0]:
                    genes[gene_id]["first_exon"] = (start, end)
            else:
                if not genes[gene_id]["last_exon"] or end > genes[gene_id]["last_exon"][1]:
                    genes[gene_id]["last_exon"] = (start, end)
        elif feature_type == "CDS":
            if strand == "+":
                if not genes[gene_id]["first_CDS"] or start < genes[gene_id]["first_CDS"][0]:
                    genes[gene_id]["first_CDS"] = (start, end)
            else:
                if not genes[gene_id]["last_CDS"] or end > genes[gene_id]["last_CDS"][1]:
                    genes[gene_id]["last_CDS"] = (start, end)

    # Second pass: Modify start and end positions based on the rules
    with open(output_file, 'w') as outfile:
        for line in lines:
            if line.startswith("#"):
                outfile.write(line)
                continue

            fields = line.strip().split("\t")
            chrom, _, feature_type, start, end, _, strand, _, attributes = fields
            start, end = int(start), int(end)
            
            # Extract gene ID again
            gene_id = None
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    gene_id = attr.split("=")[1]
                elif attr.startswith("Parent="):
                    gene_id = attr.split("=")[1]
            if not gene_id or gene_id not in genes:
                outfile.write(line)
                continue

            gene_info = genes[gene_id]

            # Skip genes at chromosome boundaries
            if (strand == "+" and start == 1) or (strand == "-" and end >= chromosome_lengths.get(chrom, float('inf')) - 1):
                outfile.write(line)
                continue

            # Extend gene, mRNA, first exon, first CDS in + strand
            if strand == "+":
                if feature_type in ["gene", "mRNA"]:
                    start = max(1, start - 4200)
                elif feature_type == "exon" and (start, end) == gene_info["first_exon"]:
                    start = max(1, start - 4200)
                elif feature_type == "CDS" and (start, end) == gene_info["first_CDS"]:
                    start = max(1, start - 4200)
                updated_lines += 1
            
            # Extend gene, mRNA, last exon, last CDS in - strand
            elif strand == "-":
                if feature_type in ["gene", "mRNA"]:
                    end += 4200
                elif feature_type == "exon" and (start, end) == gene_info["last_exon"]:
                    end += 4200
                elif feature_type == "CDS" and (start, end) == gene_info["last_CDS"]:
                    end += 4200
                updated_lines += 1

            fields[3], fields[4] = str(start), str(end)
            outfile.write("\t".join(fields) + "\n")
    
    print(f"Modified {updated_lines} features in the GFF file.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extend_coordinates.py <input_gff> <output_gff>")
        sys.exit(1)

    extend_coordinates(sys.argv[1], sys.argv[2])
