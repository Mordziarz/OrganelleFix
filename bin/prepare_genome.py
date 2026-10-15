#!/usr/bin/env python3
import os
import sys
import argparse
import shutil
import re
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def clean_gene_name(name):
    """Single copy."""
    if not name or name == "unnamed" or name.startswith("unnamed_"):
        return "unnamed"
    cleaned = re.sub(r'[^a-zA-Z0-9]', '', str(name)).lower()
    return cleaned if cleaned else "unnamed"

def is_valid_file(path):
    return path and path != "None" and os.path.exists(path) and "NO_" not in os.path.basename(path) and os.path.getsize(path) > 0

def extract_circular_seq(seq_rec, start, end, strand):
    """CDS extraction"""
    try:
        if start <= end:
            sub_seq = seq_rec.seq[start:end]
        else:
            sub_seq = seq_rec.seq[start:] + seq_rec.seq[:end]
            
        if strand == "-":
            sub_seq = sub_seq.reverse_complement()
        return sub_seq
    except Exception:
        return Seq("")

def main():
    parser = argparse.ArgumentParser(description="Genome preparation supporting GenBank and FASTA+BED.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--fasta", default=None)
    parser.add_argument("--gbk", default=None)
    parser.add_argument("--bed", default=None)
    args = parser.parse_args()

    sample_id = args.sample_id
    has_fasta = is_valid_file(args.fasta)
    has_gbk = is_valid_file(args.gbk)
    has_bed = is_valid_file(args.bed)

    if not has_fasta and not has_gbk and not has_bed:
        print(f"Error for sample {sample_id}: No valid input files provided!", file=sys.stderr)
        sys.exit(1)

    bed_rows = []

    if has_gbk:
        records = list(SeqIO.parse(args.gbk, "genbank"))
        SeqIO.write(records, "genome.fasta", "fasta")
        
        for record in records:
            chrom = record.id
            for feat in record.features:
                if feat.type in ["gene", "CDS", "tRNA", "rRNA"]:
                    if 'pseudo' in feat.qualifiers or 'pseudogene' in feat.qualifiers:
                        continue
                    raw_name = feat.qualifiers.get("gene", feat.qualifiers.get("locus_tag", ["unnamed"]))[0]
                    strand = "+" if (feat.location.strand is not None and feat.location.strand >= 0) else "-"
                    
                    parts = feat.location.parts if hasattr(feat.location, 'parts') else [feat.location]
                    for part in parts:
                        bed_rows.append((chrom, int(part.start), int(part.end), raw_name, feat.type, strand))

    elif has_fasta and has_bed:
        shutil.copy(args.fasta, "genome.fasta")
        with open(args.bed, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 6:
                    continue
                chrom, s_str, e_str, name, ftype, strand = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
                try:
                    bed_rows.append((chrom, int(s_str), int(e_str), name, ftype, strand))
                except ValueError:
                    continue
    else:
        print(f"Error for sample {sample_id}: Invalid input combination. Provide either gbk OR (fasta + bed).", file=sys.stderr)
        sys.exit(1)

    bed_rows.sort(key=lambda x: (x[0], x[1], -x[2]))
    with open("annotations.bed", "w") as bed_file:
        for chrom, start, end, name, ftype, strand in bed_rows:
            bed_file.write(f"{chrom}\t{start}\t{end}\t{name}\t{ftype}\t{strand}\n")

    genome_dict = SeqIO.to_dict(SeqIO.parse("genome.fasta", "fasta"))
    cds_entries = []

    for chrom, start, end, name, ftype, strand in bed_rows:
        if ftype.upper() == "CDS":
            if chrom in genome_dict:
                gkey = clean_gene_name(name)
                if gkey != "unnamed":
                    cds_entries.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'name': name,
                        'gkey': gkey,
                        'strand': strand
                    })

    groups = defaultdict(list)
    for entry in cds_entries:
        groups[(entry['chrom'], entry['gkey'], entry['strand'])].append(entry)

    cds_records = []
    seen_gene_keys = set()

    for (chrom, gkey, strand), entries in groups.items():
        if gkey in seen_gene_keys:
            continue
        seen_gene_keys.add(gkey)

        if strand == "-":
            entries.sort(key=lambda x: x['start'], reverse=True)
        else:
            entries.sort(key=lambda x: x['start'])

        seq_rec = genome_dict[chrom]
        full_cds_seq = Seq("")
        min_start = min(e['start'] for e in entries)
        max_end = max(e['end'] for e in entries)
        orig_name = entries[0]['name']

        for entry in entries:
            sub_exon = extract_circular_seq(seq_rec, entry['start'], entry['end'], strand)
            full_cds_seq += sub_exon

        if len(full_cds_seq) > 0:
            fid = f"{chrom}_{orig_name}_CDS_{min_start}_{max_end}_{strand}"
            desc = f"gene={orig_name} type=CDS exons={len(entries)} location=[{min_start}:{max_end}]"
            cds_records.append(SeqRecord(full_cds_seq, id=fid, description=desc))

    SeqIO.write(cds_records, "CDS.fasta", "fasta")

if __name__ == "__main__":
    main()
