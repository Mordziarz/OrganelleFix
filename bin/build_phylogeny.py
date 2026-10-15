#!/usr/bin/env python3
import argparse
import os
import subprocess
import shutil
import re
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="Build concatenated alignment and run IQ-TREE for organelle genomes.")
    parser.add_argument("--organelle", required=True, help="Organelle name (plastome/mitogenome)")
    parser.add_argument("--samples", required=True, nargs='+', help="List of sample IDs")
    parser.add_argument("--cds_files", required=True, nargs='+', help="List of CDS fasta files corresponding to samples")
    parser.add_argument("--outgroup", default=None, help="Outgroup sample ID for IQ-TREE rooting")
    return parser.parse_args()

def normalize_name(header):
    m = re.search(r'gene=([a-zA-Z0-9_-]+)', header)
    if m:
        gene_name = m.group(1)
    else:
        parts = header.replace('>', '').split('_')
        gene_name = "unnamed"
        for p in parts:
            if len(p) >= 2 and p[0].isalpha() and p[1].islower():
                gene_name = p
                break
        if gene_name == "unnamed" and len(parts) > 1:
            gene_name = parts[1]

    cleaned = re.sub(r'[^a-zA-Z0-9]', '', gene_name)
    return cleaned.lower()

def parse_cds(fasta_path):
    genes = {}
    if not os.path.exists(fasta_path):
        return genes
    current_gene, current_seq = None, []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_gene:
                    genes[current_gene] = "".join(current_seq)
                current_gene = normalize_name(line[1:])
                current_seq = []
            else:
                current_seq.append(line)
        if current_gene:
            genes[current_gene] = "".join(current_seq)
    return genes

def main():
    args = parse_arguments()
    if len(args.samples) != len(args.cds_files):
        print("[ERROR] The number of samples does not match the number of CDS files!")
        sys.exit(1)

    sample_data = {}
    for sample, cds_path in zip(args.samples, args.cds_files):
        parsed = parse_cds(cds_path)
        sample_data[sample] = parsed
        print(f"[DEBUG] Sample {sample}: found {len(parsed)} genes in {cds_path}")

    if not sample_data:
        print("[ERROR] No input data for phylogenetics.")
        sys.exit(1)

    gene_sets = [set(d.keys()) for d in sample_data.values()]
    common_genes = sorted(list(set.intersection(*gene_sets)))

    if not common_genes:
        print(f"[ERROR] No common genes found for organelle {args.organelle} across all samples! Check gene names in GBK/CDS files.")
        first_sample = list(sample_data.keys())[0]
        print(f"[DEBUG] Sample '{first_sample}' sample genes sample: {list(sample_data[first_sample].keys())[:10]}")
        sys.exit(1)

    print(f"[INFO] Found {len(common_genes)} common genes for {args.organelle}: {common_genes}")

    aligned_sequences = {sample: "" for sample in args.samples}

    for gene in common_genes:
        tmp_fa = f"temp_{gene}.fasta"
        with open(tmp_fa, 'w') as out:
            for sample in args.samples:
                seq = sample_data[sample][gene]
                out.write(f">{sample}\n{seq}\n")

        aligned_fa = f"aligned_{gene}.fasta"
        mafft_cmd = ["mafft", "--auto", "--quiet", tmp_fa]
        with open(aligned_fa, 'w') as out:
            subprocess.run(mafft_cmd, stdout=out, check=True)

        current_sample = None
        current_seq_parts = []
        with open(aligned_fa, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_sample:
                        aligned_sequences[current_sample] += "".join(current_seq_parts)
                    current_sample = line[1:]
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line)
            if current_sample:
                aligned_sequences[current_sample] += "".join(current_seq_parts)

        if os.path.exists(tmp_fa): os.remove(tmp_fa)
        if os.path.exists(aligned_fa): os.remove(aligned_fa)

    supermatrix_file = f"{args.organelle}_supermatrix.fasta"
    with open(supermatrix_file, 'w') as f:
        for sample, seq in aligned_sequences.items():
            f.write(f">{sample}\n{seq}\n")

    print(f"[INFO] Supermatrix created: {supermatrix_file}.")

    if len(args.samples) < 3:
        print(f"[WARNING] IQ-TREE requires at least 3 sequences, but only {len(args.samples)} provided. Skipping tree inference.")
        dummy_tree = f"{args.organelle}_skipped.treefile"
        with open(dummy_tree, 'w') as f:
            f.write("# Tree generation skipped: fewer than 3 samples provided.\n")
        sys.exit(0)

    iqtree_bin = "iqtree2" if shutil.which("iqtree2") else "iqtree"
    if not shutil.which(iqtree_bin):
        print("[ERROR] Neither 'iqtree2' nor 'iqtree' binary found in PATH!")
        sys.exit(1)

    iqtree_cmd = [iqtree_bin, "-s", supermatrix_file, "-st", "DNA", "-m", "MFP", "-bb", "1000", "-nt", "AUTO", "-redo"]
    
    if args.outgroup:
        if args.outgroup not in args.samples:
            print(f"[ERROR] Specified outgroup '{args.outgroup}' is not present in the samples list!")
            sys.exit(1)
        iqtree_cmd.extend(["-o", args.outgroup])
        print(f"[INFO] Using outgroup for IQ-TREE: {args.outgroup}")

    res = subprocess.run(iqtree_cmd)

    if res.returncode == 0:
        print(f"[INFO] Phylogenetic tree for {args.organelle} successfully generated!")
    else:
        print(f"[WARNING] IQ-TREE finished with exit code {res.returncode}")

if __name__ == '__main__':
    main()
