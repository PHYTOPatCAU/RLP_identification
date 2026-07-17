import sys
import csv

def load_extension_lengths(csv_file):
    extension_dict = {}
    gene_id_map = {}
    with open(csv_file, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            full_id = row['Protein_ID']
            bp = int(row['BP_To_Add'])
            extension_dict[full_id] = bp
            base_id = full_id.split('.')[0]
            gene_id_map[base_id] = bp  # Map for gene-level match
    return extension_dict, gene_id_map

def extend_coordinates_custom(gff_file, csv_file, output_file):
    extensions, gene_level_extensions = load_extension_lengths(csv_file)
    genes = {}

    with open(gff_file, 'r') as infile:
        lines = infile.readlines()

    # First pass: find first/last features
    for line in lines:
        if line.startswith("#"): continue
        fields = line.strip().split("\t")
        if len(fields) < 9: continue

        chrom, _, feature, start, end, _, strand, _, attributes = fields
        start, end = int(start), int(end)

        gene_id = None
        for attr in attributes.split(";"):
            if attr.startswith("ID="):
                gene_id = attr.split("=")[1]
            elif attr.startswith("Parent="):
                gene_id = attr.split("=")[1]

        if not gene_id: continue

        if feature == "gene":
            # For genes, strip version (.1, .2)
            base_id = gene_id.split('.')[0]
            if base_id not in gene_level_extensions:
                continue
            extension_bp = gene_level_extensions[base_id]
        else:
            if gene_id not in extensions:
                continue
            extension_bp = extensions[gene_id]

        if gene_id not in genes:
            genes[gene_id] = {"strand": strand, "features": {}, "bp": extension_bp}

        key = None
        if feature in ["gene", "mRNA"]:
            key = feature
        elif feature == "exon":
            key = "first_exon" if strand == "+" else "last_exon"
        elif feature == "CDS":
            key = "first_CDS" if strand == "+" else "last_CDS"

        if key and key not in genes[gene_id]["features"]:
            genes[gene_id]["features"][key] = (start, end)
        elif key:
            if strand == "+" and start < genes[gene_id]["features"][key][0]:
                genes[gene_id]["features"][key] = (start, end)
            elif strand == "-" and end > genes[gene_id]["features"][key][1]:
                genes[gene_id]["features"][key] = (start, end)

    # Second pass: modify and write
    with open(output_file, 'w') as out:
        for line in lines:
            if line.startswith("#"):
                out.write(line)
                continue

            fields = line.strip().split("\t")
            chrom, _, feature, start, end, _, strand, _, attributes = fields
            start, end = int(start), int(end)

            gene_id = None
            for attr in attributes.split(";"):
                if attr.startswith("ID="):
                    gene_id = attr.split("=")[1]
                elif attr.startswith("Parent="):
                    gene_id = attr.split("=")[1]

            if not gene_id:
                out.write(line)
                continue

            # Decide extension
            ext = None
            if feature in ["gene", "mRNA"]:
                base_id = gene_id.split('.')[0]
                if base_id in gene_level_extensions:
                    ext = gene_level_extensions[base_id]
            elif gene_id in extensions:
                ext = extensions[gene_id]

            if not ext or gene_id not in genes:
                out.write(line)
                continue

            gene_info = genes[gene_id]

            if strand == "+":
                if feature in ["gene", "mRNA"]:
                    start = max(1, start - ext)
                elif feature == "exon" and (start, end) == gene_info["features"].get("first_exon"):
                    start = max(1, start - ext)
                elif feature == "CDS" and (start, end) == gene_info["features"].get("first_CDS"):
                    start = max(1, start - ext)
            elif strand == "-":
                if feature in ["gene", "mRNA"]:
                    end += ext
                elif feature == "exon" and (start, end) == gene_info["features"].get("last_exon"):
                    end += ext
                elif feature == "CDS" and (start, end) == gene_info["features"].get("last_CDS"):
                    end += ext

            fields[3], fields[4] = str(start), str(end)
            out.write("\t".join(fields) + "\n")

    print(f"✅ Extension complete. Output written to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extend_by_custom_bp.py <input.gff> <extensions.csv> <output.gff>")
        sys.exit(1)
    
    extend_coordinates_custom(sys.argv[1], sys.argv[2], sys.argv[3])
