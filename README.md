# Bacterial AI-2 Quorum Sensing, *lsr* Operon, and Active Methyl (SAM) Cycle Pipeline

Modular Python workflows for screening bacterial genomes for Autoinducer-2 (AI-2) biosynthesis (`luxS`), active methyl (SAM) cycle enzymes (`metK`, `mtnN`, `metE`, `metH`), and *Escherichia coli*-type *lsr* uptake/regulatory loci (`lsrACDBFG`, `lsrR`, `lsrK`), with separate entry points for **Gram-positive** and **Gram-negative** genomes.

---

## Overview

Bacterial AI-2 signaling and SAM metabolism differ across taxonomic groups:

* **Gram-Positive Workflow (`run_gram_positive.py`):** Targets `luxS` and active methyl cycle genes (`metK`, `mtnN`, `metE`, `metH`) across Gram-positive reference taxa (*Bacillus*, *Staphylococcus*, *Streptococcus*, *Listeria*, *Enterococcus*), extracts flanking ORFs around each `luxS` locus (`window_bp=10000` by default), and generates a 300 dpi neighborhood schematic map with non-overlapping tracks. Many Gram-positive organisms (e.g., *Bacillus velezensis*) carry `luxS` and SAM-cycle genes without an *E. coli*-style `lsrACDBFG` operon.
* **Gram-Negative Workflow (`run_gram_negative.py`):** Targets *E. coli* K-12 MG1655 `luxS`, `lsrACDBFG`, and `lsrR`/`lsrK` while filtering cross-homologous sugar ABC transporters (`rbsA/B/C`, `mglA/B/C`, `xylF/G/H`) and FGGY carbohydrate kinases (`xylB`, `glpK`, `rhaB`, `araB`) using tagged `DECOY` reference entries. It also checks multi-gene cluster synteny (`annotate_synteny_status`), computes intergenic distances and strand orientations, scans upstream regions for a **candidate/heuristic** palindromic motif (`AATGAA.{6,8}TTCATT`) on a single strand, and renders 300 dpi operon and spacing plots.

---

## Requirements & Setup

### Python & External Tools
* **Python:** 3.10+
* **Python packages:** `biopython`, `pandas`, `numpy`, `matplotlib` (optional: `pyrodigal` for cleaner prokaryotic gene calling)
* **NCBI BLAST+ CLI:** `makeblastdb` and `blastp` must be available in your system `PATH`.
* **Network access (first run / cache refresh):** Required only when building or refreshing NCBI Entrez reference FASTA files.

### Environment Setup (Windows PowerShell)

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install biopython pandas numpy matplotlib

makeblastdb -version
blastp -version
python -c "import Bio, pandas, numpy, matplotlib; print('Environment ready.')"
python run_gram_positive.py
Key generated files:

gram_neg_ai2_k12_decoys_ref.fasta — Cached, gene=-tagged E. coli K-12 target and DECOY reference sequences

gram_neg_ai2_blast.tsv — Raw tabular BLASTP output

gram_neg_ai2_annotations.csv — Filtered, non-decoy AI-2 / lsr homologs

gram_neg_extracted_ai2_proteins.faa — Extracted protein sequences for matched ORFs

gram_neg_lsr_intergenic_distances.csv — Neighboring gene spacing (bp) and orientation (Tandem, Convergent, Divergent)

gram_neg_lsr_upstream_promoters.csv & gram_neg_lsr_candidate_motifs.csv — Upstream sequences and single-strand heuristic motif matches

lsr_operon_genomic_map.png & lsr_intergenic_distances_bar_chart.png — 300 dpi visualizations

Key Safeguards & Notes
Cache validation (_cache_ok): Reuses existing reference FASTA files when force_refresh=False only if files exist, are non-empty, carry gene= tags, and include all expected gene tags; otherwise re-fetches cleanly.

Coordinate & BLAST deduplication: ORF coordinates are parsed into forward genomic space once (all_orfs_df), and BLASTP hits are deduplicated by highest bitscore per ORF_ID before joining coordinates, preventing multi-HSP row multiplication and self-pair artifacts.

Biological scope: Homology-based detection indicates genetic potential rather than confirmed metabolite production or regulatory activity; metK, mtnN, metE, and metH are broadly distributed across bacteria and do not by themselves establish Gram staining.
