# OrganelleFix
A Nextflow pipeline for organellar genome analysis.

# Usage

1. Install Nextflow via Conda.
2. Ensure Docker is installed and running on your system.
3. Build `samples_table`

## Samples Configuration File (-samples_table)

The `-samples_table` input file controls the execution of the **OrganelleFix** pipeline. It allows defining organellar samples and supports two alternative data delivery methods (GenBank **or** FASTA + BED).

## Column Description:
* **`sample_id`**: Unique identifier for the sample or species (used in output folder names).
* **`organelle`**: Type of organellar genome (e.g., `plastome` or `mitogenome`).
* **`comparison_group`**: Comparative or taxonomic group (e.g., `Brassicaceae`), useful for downstream bulk analyses.
* **`fasta_path`**: Path to the genome FASTA file (required in FASTA+BED mode; leave empty if using a GenBank file).
* **`gbk_path`**: Path to the GenBank `.gb` file (required in GBK mode; leave empty if using the FASTA + BED pair).
* **`bed_path`**: Path to the annotations BED file (required in FASTA+BED mode; leave empty if using a GenBank file).

---

## Example Content of `-samples_table`:

```csv
sample_id,organelle,comparison_group,fasta_path,gbk_path,bed_path
Species_1,plastome,Group_A,,data_test/Species_1_plastome.gb,
Species_1,mitogenome,Group_A,,data_test/Species_1_mitogenome.gb,
Species_2,plastome,Group_A,,data_test/Species_2_plastome.gb,
Species_2,mitogenome,Group_A,,data_test/Species_2_mitogenome.gb,
Species_3,plastome,Group_B,/path/to/Species_3_plastome.fasta,,/path/to/Species_3_plastome.bed
Species_3,mitogenome,Group_B,/path/to/Species_3_mitogenome.fasta,,/path/to/Species_3_mitogenome.bed
```

4. Run pipeline

## Pipeline Parameters

- **`--samples_table`** (`file`, required): Path to input CSV samples sheet.
- **`--ks_ka_reference_plastome`** (`file`, optional): Reference plastome sequence for Ka/Ks analysis.
- **`--ks_ka_reference_mitogenome`** (`file`, optional): Reference mitogenome sequence for Ka/Ks analysis.
- **`--outgroup_plastome`** (`file`, optional): Outgroup sequence for plastome phylogenetic tree building.
- **`--outgroup_mitogenome`** (`file`, optional): Outgroup sequence for mitogenome phylogenetic tree building.

``` bash
nextflow run mordziarz/OrganelleFix -profile standard  --samples_table samples_table.csv -profile standard --ks_ka_reference_plastome "species" --ks_ka_reference_mitogenome "species" --outgroup_plastome "species" --outgroup_mitogenome "species"
```
# WARNING
The pipeline does **not** automatically validate the number of CDS features, nor whether they start and end correctly. It is **strongly recommended** to inspect all output files in the `results/prepared_genomes` folder. While downstream analyses based on nucleotide sequences are generally robust, verifying the correct CDS boundaries for your dataset is essential.

# Results

```text
results/
├── kaks/
├── phylogeny/
├── prepared_genomes/
├── repeats/
├── rna_editing/
├── rscu_results/
├── ssr/
└── transgenomes_results/
```

## RSCU Analysis

## SSR Analysis

## Repeat Analysis 

## Ka/Ks Analysis

## Phylogenetic Tree (CDS-based)

## RNA Editing Prediction

## Organellar Transfers (Mitogenome vs. Plastome)

# Citation

If you use **OrganelleFix** in your research, please cite the underlying tools and software packages integrated into this pipeline:

- **RSCU Analysis:** RSCUcaller ([10.1186/s12859-025-06166-5](https://doi.org/10.1186/s12859-025-06166-5))
- **Transfer Analysis:** transgenomes ([10.1007/s10265-026-01709-0](https://doi.org/10.1007/s10265-026-01709-0))
- **Phylogenetic Analysis:** IQ-TREE ([10.1093/molbev/msag117](https://doi.org/10.1093/molbev/msag117))
- **Alignment:** MAFFT ([10.1093/nar/gkf436](https://doi.org/10.1093/nar/gkf436))
- **Visualization:** ggplot2 ([10.1007/978-3-319-24277-4](https://doi.org/10.1007/978-3-319-24277-4))

> *Note:* A dedicated scientific publication for **OrganelleFix** is currently in preparation.