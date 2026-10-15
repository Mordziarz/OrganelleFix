#!/usr/bin/env python3
import argparse
import os
import csv
import re

GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

HYDROPHOBIC_AAS = {'I', 'L', 'V', 'F', 'M', 'W'}

def parse_arguments():
    parser = argparse.ArgumentParser(description="Universal biochemical RNA editing prediction for plant organelle CDS.")
    parser.add_argument("--sample", required=True, help="Sample name")
    parser.add_argument("--organelle", required=True, help="Organelle type (plastome/mitogenome)")
    parser.add_argument("--cds", required=True, help="CDS fasta file")
    return parser.parse_args()

def parse_fasta(fasta_path):
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
                header = line[1:]
                m = re.search(r'gene=([A-Za-z0-9_-]+)', header)
                current_gene = m.group(1) if m else header.split()[0].split('_')[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_gene:
            genes[current_gene] = "".join(current_seq)
    return genes

def universal_smart_filter(gene_name, seq):
    """Universal and rigorous biochemical filter (C-to-U with hydrophobicity preference at the 1st or 2nd codon position)"""
    editing_events = []
    seq = seq.upper()
    codon_len = (len(seq) // 3) * 3
    
    for i in range(0, codon_len, 3):
        codon = seq[i:i+3]
        if 'N' in codon or len(codon) < 3:
            continue
        orig_aa = GENETIC_CODE.get(codon, 'X')
        
        for pos_in_codon in range(3):
            if codon[pos_in_codon] == 'C':
                mut_codon_list = list(codon)
                mut_codon_list[pos_in_codon] = 'T'
                mut_codon = "".join(mut_codon_list)
                mut_aa = GENETIC_CODE.get(mut_codon, 'X')
                
                # Rygorystyczny filtr biologiczny dla roślin organellarnych
                if orig_aa != mut_aa and mut_aa != '*' and mut_aa in HYDROPHOBIC_AAS and pos_in_codon in [0, 1]:
                    editing_events.append({
                        'gene': gene_name,
                        'codon_position': (i // 3) + 1,
                        'nucleotide_position': i + pos_in_codon + 1,
                        'original_codon': codon,
                        'edited_codon': mut_codon,
                        'original_aa': orig_aa,
                        'edited_aa': mut_aa,
                        'type': 'C-to-U',
                        'confidence_score': 'High'
                    })
    return editing_events

def main():
    args = parse_arguments()
    out_file = f"{args.sample}_{args.organelle}_rna_editing_predictions.csv"
    
    if not os.path.exists(args.cds):
        print(f"[ERROR] Plik CDS nie istnieje: {args.cds}")
        return

    genes = parse_fasta(args.cds)
    all_predictions = []
    
    for gene, seq in genes.items():
        if gene.startswith("NC_") and len(gene) > 10:
            continue
        all_predictions.extend(universal_smart_filter(gene, seq))
        
    fieldnames = ['gene', 'codon_position', 'nucleotide_position', 'original_codon', 'edited_codon', 'original_aa', 'edited_aa', 'type', 'confidence_score']
    with open(out_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=';')
        writer.writeheader()
        for row in all_predictions:
            writer.writerow(row)
            
    print(f"[INFO] Universal filter scan: Found {len(all_predictions)} high-confidence editing sites for {args.sample} ({args.organelle}). Saved to {out_file}")

if __name__ == '__main__':
    main()
