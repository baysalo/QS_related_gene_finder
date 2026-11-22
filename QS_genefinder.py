import json
import os
import re
import subprocess
import time
from typing import Optional
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ==============================================================================
# 1. GENOME & ORF COORDINATE PARSING
# ==============================================================================
def load_contig_lengths(genome_fasta: str) -> dict:
    lengths = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        lengths[record.id.split()[0]] = len(record.seq)
    return lengths


def _match_contig_name(orf_full_id: str, contig_lengths: dict) -> Optional[str]:
    if "_" in orf_full_id:
        candidate = orf_full_id.rsplit("_", 1)[0]
        if candidate in contig_lengths:
            return candidate
    matches = [cname for cname in contig_lengths if orf_full_id.startswith(cname)]
    return max(matches, key=len) if matches else None


def parse_orf_coords(orf_records_or_ids: list, contig_lengths: dict) -> pd.DataFrame:
    """
    Supports:
      1) <Contig>_ORF_<Idx>_<RawStart>_<RawEnd>_<fwd|rev>
      2) Prodigal/Pyrodigal description: <ID> # <Start> # <End> # <1|-1> # ...
    """
    pat_custom = re.compile(r"^(.*)_ORF_(\d+)_(\d+)_(\d+)_(fwd|rev)$")
    pat_prodigal = re.compile(r"^(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-1|1)\b")

    rows = []
    skipped = 0

    for item in orf_records_or_ids:
        if hasattr(item, "id") and hasattr(item, "description"):
            oid, desc = item.id, item.description
        else:
            oid, desc = str(item), str(item)

        contig = g_start = g_end = strand = None

        m1 = pat_custom.match(oid)
        if m1:
            contig = m1.group(1)
            raw_start = int(m1.group(3))
            raw_end = int(m1.group(4))
            strand = m1.group(5)
            c_len = contig_lengths.get(contig)
            if c_len is not None:
                if strand == "fwd":
                    g_start, g_end = raw_start, raw_end
                else:
                    g_start = c_len - raw_end + 1
                    g_end = c_len - raw_start + 1
        else:
            m2 = pat_prodigal.match(desc)
            if m2:
                contig = _match_contig_name(m2.group(1), contig_lengths)
                g_start = int(m2.group(2))
                g_end = int(m2.group(3))
                strand = "fwd" if m2.group(4) == "1" else "rev"

        if contig is None or g_start is None or g_end is None or contig not in contig_lengths:
            skipped += 1
            continue

        rows.append({
            "ORF_ID": oid,
            "Contig": contig,
            "Genomic_Start": min(g_start, g_end),
            "Genomic_End": max(g_start, g_end),
            "Strand": strand,
            "Forward": strand == "fwd",
        })

    if skipped > 0:
        print(f"Warning: {skipped} ORF records did not match expected coordinate headers and were skipped.")

    return pd.DataFrame(
        rows,
        columns=["ORF_ID", "Contig", "Genomic_Start", "Genomic_End", "Strand", "Forward"],
    )


def predict_bacterial_orfs_fallback(
    genome_fasta: str,
    out_faa: str = "predicted_orfs_clean.faa",
    min_aa_len: int = 80,
    table_id: int = 11,
    fmet_start: bool = True,
) -> int:
    start_codons = {"ATG", "GTG", "TTG"}
    stop_codons = {"TAA", "TAG", "TGA"}
    orf_records = []
    idx = 1

    for record in SeqIO.parse(genome_fasta, "fasta"):
        contig = record.id.split()[0]
        fwd_seq = str(record.seq).upper()
        rev_seq = str(record.seq.reverse_complement()).upper()

        for strand_label, seq_str in [("fwd", fwd_seq), ("rev", rev_seq)]:
            seq_len = len(seq_str)
            for frame in range(3):
                pos = frame
                current_start = None
                while pos + 2 < seq_len:
                    codon = seq_str[pos : pos + 3]
                    if current_start is None:
                        if codon in start_codons:
                            current_start = pos
                    elif codon in stop_codons:
                        aa_len = (pos - current_start) // 3
                        if aa_len >= min_aa_len:
                            raw_start = current_start + 1
                            raw_end = pos
                            aa_seq = str(
                                Seq(seq_str[current_start:pos]).translate(table=table_id, to_stop=True)
                            )
                            if fmet_start and len(aa_seq) > 0 and aa_seq[0] != "M":
                                aa_seq = "M" + aa_seq[1:]
                            orf_id = f"{contig}_ORF_{idx}_{raw_start}_{raw_end}_{strand_label}"
                            orf_records.append(
                                SeqRecord(Seq(aa_seq), id=orf_id, description=f"frame={frame + 1} len_aa={aa_len}")
                            )
                            idx += 1
                        current_start = None
                    pos += 3

    SeqIO.write(orf_records, out_faa, "fasta")
    print(f"Wrote {len(orf_records)} ORFs to '{out_faa}' (6-frame fallback; Pyrodigal is recommended when available).")
    return len(orf_records)


def ensure_orfs_faa(genome_fasta: str, out_faa: str = "predicted_orfs_clean.faa", force_refresh: bool = False) -> str:
    if not force_refresh and os.path.exists(out_faa) and os.path.getsize(out_faa) > 0:
        return out_faa
    try:
        import pyrodigal

        finder = pyrodigal.GeneFinder(meta=False)
        orf_records = []
        idx = 1
        for record in SeqIO.parse(genome_fasta, "fasta"):
            contig = record.id.split()[0]
            genes = finder.find_genes(str(record.seq).upper().encode("ascii"))
            for g in genes:
                strand_label = "fwd" if g.strand == 1 else "rev"
                orf_id = f"{contig}_ORF_{idx}_{int(g.begin)}_{int(g.end)}_{strand_label}"
                orf_records.append(SeqRecord(Seq(g.translate()), id=orf_id, description="pyrodigal"))
                idx += 1
        SeqIO.write(orf_records, out_faa, "fasta")
        print(f"Predicted {len(orf_records)} non-overlapping CDSs with Pyrodigal into '{out_faa}'.")
    except ImportError:
        predict_bacterial_orfs_fallback(genome_fasta, out_faa=out_faa)
    return out_faa


# ==============================================================================
# 2. GRAM-POSITIVE TAGGED REFERENCE BUILDER & BLASTP
# ==============================================================================
def _entrez_fetch(query: str, out_fasta: str, retmax: int = 100) -> int:
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urlencode({
        "db": "protein",
        "term": query,
        "usehistory": "y",
        "retmax": retmax,
        "retmode": "json",
    })
    with urlopen(esearch_url, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))

    esr = data.get("esearchresult", {})
    total_count = int(esr.get("count", 0))
    webenv = esr.get("webenv", "")
    query_key = esr.get("querykey", "")
    if total_count == 0 or not webenv or not query_key:
        raise RuntimeError(f"NCBI esearch returned no hits (count={total_count}).")

    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode({
        "db": "protein",
        "query_key": query_key,
        "WebEnv": webenv,
        "retmax": retmax,
        "rettype": "fasta",
        "retmode": "text",
    })
    with urlopen(efetch_url, timeout=60) as resp, open(out_fasta, "wb") as fh:
        fh.write(resp.read())

    return sum(1 for _ in SeqIO.parse(out_fasta, "fasta"))


def tag_from_title(stitle: str) -> str:
    m = re.search(r"gene=(\S+)", str(stitle))
    return m.group(1) if m else "?"


def _cache_ok(out_fasta: str, jobs: list) -> bool:
    if not (os.path.exists(out_fasta) and os.path.getsize(out_fasta) > 0):
        return False
    expected = {tag for _, gm in jobs for tag in gm.values()}
    found = set()
    for r in SeqIO.parse(out_fasta, "fasta"):
        tag = tag_from_title(r.description)
        if tag == "?":
            print(f"Cache '{out_fasta}' lacks gene= tags; re-downloading.")
            return False
        found.add(tag)
    missing = sorted(expected - found)
    if missing:
        print(f"Cache '{out_fasta}' is missing tags {missing}; re-downloading.")
        return False
    return True


def build_tagged_refs(jobs: list, out_fasta: str, retmax_per_gene: int = 100, force_download: bool = False) -> int:
    if not force_download and _cache_ok(out_fasta, jobs):
        count = sum(1 for _ in SeqIO.parse(out_fasta, "fasta"))
        print(f"Using validated cache '{out_fasta}' ({count} records).")
        return count

    records = []
    seen = set()
    for taxon_q, gene_map in jobs:
        for gene, tag in gene_map.items():
            tmp = f"_tmp_gp_{gene}_{os.getpid()}.fasta"
            try:
                _entrez_fetch(f"{taxon_q} AND {gene}[Gene] AND refseq[filter]", tmp, retmax=retmax_per_gene)
                for r in SeqIO.parse(tmp, "fasta"):
                    if r.id in seen:
                        continue
                    seen.add(r.id)
                    rest = r.description[len(r.id) :].strip()
                    r.description = f"{r.id} gene={tag} {rest}".strip()
                    records.append(r)
            except (RuntimeError, URLError, TimeoutError) as e:
                print(f"Warning: failed to fetch '{gene}' ({taxon_q}): {e}")
            finally:
                if os.path.exists(tmp):
                    os.remove(tmp)
                time.sleep(0.4)

    SeqIO.write(records, out_fasta, "fasta")
    print(f"Wrote {len(records)} tagged Gram-positive reference sequences to '{out_fasta}'.")
    return len(records)


def run_gram_pos_blastp(
    query_faa: str = "predicted_orfs_clean.faa",
    ref_fasta: str = "gram_pos_luxs_methyl_ref.fasta",
    blast_tsv: str = "gram_pos_ai2_methyl_blast.tsv",
    out_csv: str = "gram_pos_luxs_methyl_summary.csv",
    min_pident: float = 30.0,
    min_aln_len: int = 40,
    min_scov: float = 0.50,
) -> pd.DataFrame:
    subprocess.run(["makeblastdb", "-in", ref_fasta, "-dbtype", "prot"], check=True)
    subprocess.run(
        [
            "blastp",
            "-query", query_faa,
            "-db", ref_fasta,
            "-out", blast_tsv,
            "-evalue", "1e-5",
            "-max_target_seqs", "5",
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore qlen slen stitle",
        ],
        check=True,
    )

    cols = [
        "ORF_ID", "Subject_ID", "Percent_Identity",
        "Alignment_Length", "E_value", "Bitscore",
        "Query_Len", "Subject_Len", "Annotation",
    ]
    if not os.path.exists(blast_tsv) or os.path.getsize(blast_tsv) == 0:
        return pd.DataFrame(columns=cols + ["Subject_Cov", "Gene_Symbol"])

    df = pd.read_csv(blast_tsv, sep="\t", names=cols)
    df["Subject_Cov"] = df["Alignment_Length"] / df["Subject_Len"].replace(0, float("nan"))

    hits = df[
        (df["Percent_Identity"] >= min_pident)
        & (df["Alignment_Length"] >= min_aln_len)
        & (df["Subject_Cov"] >= min_scov)
    ].copy()

    hits = (
        hits.sort_values(["Bitscore", "E_value"], ascending=[False, True])
        .drop_duplicates(subset=["ORF_ID"], keep="first")
        .reset_index(drop=True)
    )
    hits["Gene_Symbol"] = hits["Annotation"].apply(tag_from_title)
    hits = hits[hits["Gene_Symbol"] != "?"].sort_values(["Gene_Symbol", "E_value"]).reset_index(drop=True)
    hits.to_csv(out_csv, index=False)

    print("\n=== Gram-Positive LuxS & Active Methyl Cycle Summary ===")
    if hits.empty:
        print("No luxS or methyl cycle homologs passed filters.")
    else:
        print("Counts per gene:", hits.groupby("Gene_Symbol")["ORF_ID"].nunique().to_dict())
        print(hits[["ORF_ID", "Gene_Symbol", "Percent_Identity", "Subject_Cov", "E_value"]].to_string(index=False))

    return hits


# ==============================================================================
# 3. LUXS NEIGHBORHOOD TABLE & 300 DPI SCHEMATIC MAP
# ==============================================================================
def summarize_luxs_neighborhood(
    all_orfs_df: pd.DataFrame,
    hits_df: pd.DataFrame,
    window_bp: int = 10000,
    out_csv: str = "luxs_flanking_neighborhood.csv",
) -> pd.DataFrame:
    if hits_df.empty or all_orfs_df.empty:
        return pd.DataFrame()

    all_unique = all_orfs_df.drop_duplicates(subset=["ORF_ID"]).copy()
    coords = all_unique.merge(
        hits_df[["ORF_ID", "Gene_Symbol", "Percent_Identity", "Subject_Cov", "E_value"]],
        on="ORF_ID",
        how="inner",
    )
    luxs_rows = coords[coords["Gene_Symbol"] == "luxS"]
    if luxs_rows.empty:
        print("No 'luxS' hit found for neighborhood summary.")
        return pd.DataFrame()

    sym_map = dict(zip(coords["ORF_ID"], coords["Gene_Symbol"]))
    frames = []

    for idx, (_, lrow) in enumerate(luxs_rows.iterrows(), start=1):
        contig = lrow["Contig"]
        w_start = max(1, int(lrow["Genomic_Start"]) - window_bp)
        w_end = int(lrow["Genomic_End"]) + window_bp

        sub = all_unique[
            (all_unique["Contig"] == contig)
            & (all_unique["Genomic_End"] >= w_start)
            & (all_unique["Genomic_Start"] <= w_end)
        ].sort_values(["Genomic_Start", "Genomic_End"]).copy()

        sub["LuxS_Focus_ID"] = lrow["ORF_ID"]
        sub["LuxS_Focus_Idx"] = idx
        sub["Pathway_Role"] = sub["ORF_ID"].map(sym_map).fillna("Flanking_ORF")
        sub["Same_Strand_As_LuxS"] = sub["Strand"] == lrow["Strand"]
        frames.append(sub)

    out_df = pd.concat(frames, ignore_index=True)
    out_df.to_csv(out_csv, index=False)
    print(f"Saved {len(out_df)} neighborhood ORF rows (+/- {window_bp} bp) to '{out_csv}'.")
    return out_df


def _assign_tracks(genes_df: pd.DataFrame, pad_bp: float = 80.0) -> pd.DataFrame:
    if genes_df.empty:
        return genes_df.assign(Track=pd.Series(dtype=int))
    ordered = genes_df.sort_values(["Genomic_Start", "Genomic_End"]).copy()
    track_ends = []
    tracks = []
    for _, row in ordered.iterrows():
        s, e = float(row["Genomic_Start"]), float(row["Genomic_End"])
        placed = None
        for t_idx, t_end in enumerate(track_ends):
            if s > t_end + pad_bp:
                track_ends[t_idx] = e
                placed = t_idx
                break
        if placed is None:
            placed = len(track_ends)
            track_ends.append(e)
        tracks.append(placed)
    ordered["Track"] = tracks
    return ordered


def plot_luxs_neighborhood_map(
    all_orfs_df: pd.DataFrame,
    hits_df: pd.DataFrame,
    window_bp: int = 10000,
    out_png: str = "gram_pos_luxs_neighborhood_map.png",
):
    if all_orfs_df.empty or hits_df.empty:
        return

    all_unique = all_orfs_df.drop_duplicates(subset=["ORF_ID"]).copy()
    coords = all_unique.merge(
        hits_df[["ORF_ID", "Gene_Symbol", "Percent_Identity"]],
        on="ORF_ID",
        how="inner",
    )
    luxs_rows = coords[coords["Gene_Symbol"] == "luxS"].sort_values(["Contig", "Genomic_Start"])
    if luxs_rows.empty:
        return

    sym_map = dict(zip(coords["ORF_ID"], coords["Gene_Symbol"]))
    pid_map = dict(zip(coords["ORF_ID"], coords["Percent_Identity"]))

    palette = {
        "luxS": "#F4A261",
        "metK": "#2A9D8F",
        "mtnN": "#457B9D",
        "metE": "#8AB17D",
        "metH": "#6D597A",
        "Flanking": "#E0E0E0",
    }

    n_panels = len(luxs_rows)
    fig, axes = plt.subplots(n_panels, 1, figsize=(13, max(3.6, 3.0 * n_panels)), squeeze=False)

    for p_idx, (_, lrow) in enumerate(luxs_rows.iterrows()):
        ax = axes[p_idx, 0]
        contig = lrow["Contig"]
        focus_id = lrow["ORF_ID"]
        w_start = max(1, int(lrow["Genomic_Start"]) - window_bp)
        w_end = int(lrow["Genomic_End"]) + window_bp
        span = max(1, w_end - w_start)

        sub = all_unique[
            (all_unique["Contig"] == contig)
            & (all_unique["Genomic_End"] >= w_start)
            & (all_unique["Genomic_Start"] <= w_end)
        ].copy()
        sub["Gene_Symbol"] = sub["ORF_ID"].map(sym_map).fillna("Flanking")
        sub["Percent_Identity"] = sub["ORF_ID"].map(pid_map)
        sub = _assign_tracks(sub, pad_bp=span * 0.006)

        n_tracks = max(1, int(sub["Track"].max()) + 1 if not sub.empty else 1)
        step = 0.42

        for _, gene in sub.iterrows():
            x0, x1 = float(gene["Genomic_Start"]), float(gene["Genomic_End"])
            width = max(1.0, x1 - x0)
            sym = gene["Gene_Symbol"]
            color = palette.get(sym, "#E0E0E0")
            is_focus = gene["ORF_ID"] == focus_id

            y_mid = 0.40 + int(gene["Track"]) * step
            half_h = 0.11
            head = min(width * 0.25, span * 0.015)

            if gene["Forward"]:
                body_end = max(x0, x1 - head)
                pts = [(x0, y_mid - half_h), (body_end, y_mid - half_h), (x1, y_mid), (body_end, y_mid + half_h), (x0, y_mid + half_h)]
            else:
                body_start = min(x1, x0 + head)
                pts = [(x1, y_mid - half_h), (body_start, y_mid - half_h), (x0, y_mid), (body_start, y_mid + half_h), (x1, y_mid + half_h)]

            ax.add_patch(
                Polygon(
                    pts,
                    closed=True,
                    facecolor=color,
                    edgecolor="#B22222" if is_focus else "black",
                    linewidth=1.5 if is_focus else 0.7,
                    linestyle="--" if is_focus else "-",
                )
            )
            if sym != "Flanking":
                pid_s = f" ({gene['Percent_Identity']:.0f}%)" if pd.notna(gene["Percent_Identity"]) else ""
                ax.text(
                    (x0 + x1) / 2.0,
                    y_mid + half_h + 0.04,
                    f"{sym}{pid_s}" + (" [focus]" if is_focus else ""),
                    ha="center",
                    va="bottom",
                    fontsize=8.5,
                    fontweight="bold",
                )

        ax.set_xlim(w_start, w_end)
        ax.set_ylim(0.15, 0.35 + n_tracks * step + 0.25)
        ax.set_yticks([])
        ax.set_xlabel("Genomic Coordinate (bp)")
        ax.set_title(f"Gram-Positive luxS (+/- {window_bp:,} bp) — {contig}", loc="left", fontsize=10.5, fontweight="bold")
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)

    legend_handles = [
        Patch(facecolor=palette["luxS"], edgecolor="#B22222", linestyle="--", label="luxS (AI-2 / DPD)"),
        Patch(facecolor=palette["metK"], edgecolor="black", label="metK (Met -> SAM)"),
        Patch(facecolor=palette["mtnN"], edgecolor="black", label="mtnN (SAH -> SRH)"),
        Patch(facecolor=palette["metE"], edgecolor="black", label="metE (B12-indep. Met)"),
        Patch(facecolor=palette["metH"], edgecolor="black", label="metH (B12-dep. Met)"),
        Patch(facecolor=palette["Flanking"], edgecolor="black", label="Neighboring ORF"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3, fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved 300 dpi map '{out_png}'.")


if __name__ == "__main__":
    genome_fa = "bacterial_genome.fasta"
    orfs_faa = "predicted_orfs_clean.faa"
    ref_fa = "gram_pos_luxs_methyl_ref.fasta"
    force_refresh = False

    contig_lengths = load_contig_lengths(genome_fa)
    ensure_orfs_faa(genome_fa, out_faa=orfs_faa, force_refresh=force_refresh)
    orf_records = list(SeqIO.parse(orfs_faa, "fasta"))
    all_orfs_df = parse_orf_coords(orf_records, contig_lengths)

    gp_taxon = (
        "(Bacillus[Organism] OR Staphylococcus[Organism] OR Streptococcus[Organism] "
        "OR Listeria[Organism] OR Enterococcus[Organism])"
    )
    gp_genes = {
        "luxS": "luxS",
        "metK": "metK",
        "mtnN": "mtnN",
        "metE": "metE",
        "metH": "metH",
    }

    build_tagged_refs(jobs=[(gp_taxon, gp_genes)], out_fasta=ref_fa, retmax_per_gene=100, force_download=force_refresh)
    hits_df = run_gram_pos_blastp(query_faa=orfs_faa, ref_fasta=ref_fa)
    if not hits_df.empty:
        summarize_luxs_neighborhood(all_orfs_df, hits_df, window_bp=10000)
        plot_luxs_neighborhood_map(all_orfs_df, hits_df, window_bp=10000)

  import json
import os
import re
import subprocess
import time
from typing import Optional
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# ==============================================================================
# 1. GENOME & ORF COORDINATE PARSING
# ==============================================================================
def load_contig_lengths(genome_fasta: str) -> dict:
    return {rec.id.split()[0]: len(rec.seq) for rec in SeqIO.parse(genome_fasta, "fasta")}


def load_genome_seqs(genome_fasta: str) -> dict:
    return {rec.id.split()[0]: str(rec.seq).upper() for rec in SeqIO.parse(genome_fasta, "fasta")}


def _match_contig_name(orf_full_id: str, contig_lengths: dict) -> Optional[str]:
    if "_" in orf_full_id:
        candidate = orf_full_id.rsplit("_", 1)[0]
        if candidate in contig_lengths:
            return candidate
    matches = [cname for cname in contig_lengths if orf_full_id.startswith(cname)]
    return max(matches, key=len) if matches else None


def parse_orf_coords(orf_records_or_ids: list, contig_lengths: dict) -> pd.DataFrame:
    pat_custom = re.compile(r"^(.*)_ORF_(\d+)_(\d+)_(\d+)_(fwd|rev)$")
    pat_prodigal = re.compile(r"^(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-1|1)\b")

    rows = []
    for item in orf_records_or_ids:
        if hasattr(item, "id") and hasattr(item, "description"):
            oid, desc = item.id, item.description
        else:
            oid, desc = str(item), str(item)

        contig = g_start = g_end = strand = None
        m1 = pat_custom.match(oid)
        if m1:
            contig = m1.group(1)
            raw_start = int(m1.group(3))
            raw_end = int(m1.group(4))
            strand = m1.group(5)
            c_len = contig_lengths.get(contig)
            if c_len is not None:
                if strand == "fwd":
                    g_start, g_end = raw_start, raw_end
                else:
                    g_start = c_len - raw_end + 1
                    g_end = c_len - raw_start + 1
        else:
            m2 = pat_prodigal.match(desc)
            if m2:
                contig = _match_contig_name(m2.group(1), contig_lengths)
                g_start = int(m2.group(2))
                g_end = int(m2.group(3))
                strand = "fwd" if m2.group(4) == "1" else "rev"

        if contig is None or g_start is None or g_end is None or contig not in contig_lengths:
            continue

        rows.append({
            "ORF_ID": oid,
            "Contig": contig,
            "Genomic_Start": min(g_start, g_end),
            "Genomic_End": max(g_start, g_end),
            "Strand": strand,
            "Forward": strand == "fwd",
        })

    return pd.DataFrame(rows, columns=["ORF_ID", "Contig", "Genomic_Start", "Genomic_End", "Strand", "Forward"])


def ensure_orfs_faa(genome_fasta: str, out_faa: str = "predicted_orfs_clean.faa", force_refresh: bool = False) -> str:
    if not force_refresh and os.path.exists(out_faa) and os.path.getsize(out_faa) > 0:
        return out_faa
    try:
        import pyrodigal
        from Bio.SeqRecord import SeqRecord

        finder = pyrodigal.GeneFinder(meta=False)
        orf_records = []
        idx = 1
        for record in SeqIO.parse(genome_fasta, "fasta"):
            contig = record.id.split()[0]
            for g in finder.find_genes(str(record.seq).upper().encode("ascii")):
                strand_label = "fwd" if g.strand == 1 else "rev"
                orf_id = f"{contig}_ORF_{idx}_{int(g.begin)}_{int(g.end)}_{strand_label}"
                orf_records.append(SeqRecord(Seq(g.translate()), id=orf_id, description="pyrodigal"))
                idx += 1
        SeqIO.write(orf_records, out_faa, "fasta")
        print(f"Predicted {len(orf_records)} CDSs with Pyrodigal into '{out_faa}'.")
    except ImportError:
        from Bio.SeqRecord import SeqRecord

        start_codons = {"ATG", "GTG", "TTG"}
        stop_codons = {"TAA", "TAG", "TGA"}
        orf_records = []
        idx = 1
        for record in SeqIO.parse(genome_fasta, "fasta"):
            contig = record.id.split()[0]
            for strand_label, seq_str in [("fwd", str(record.seq).upper()), ("rev", str(record.seq.reverse_complement()).upper())]:
                for frame in range(3):
                    pos = frame
                    cur_start = None
                    while pos + 2 < len(seq_str):
                        codon = seq_str[pos : pos + 3]
                        if cur_start is None:
                            if codon in start_codons:
                                cur_start = pos
                        elif codon in stop_codons:
                            if (pos - cur_start) // 3 >= 80:
                                aa_seq = str(Seq(seq_str[cur_start:pos]).translate(table=11, to_stop=True))
                                if aa_seq and aa_seq[0] != "M":
                                    aa_seq = "M" + aa_seq[1:]
                                orf_id = f"{contig}_ORF_{idx}_{cur_start + 1}_{pos}_{strand_label}"
                                orf_records.append(SeqRecord(Seq(aa_seq), id=orf_id, description=f"frame={frame + 1}"))
                                idx += 1
                            cur_start = None
                        pos += 3
        SeqIO.write(orf_records, out_faa, "fasta")
        print(f"Wrote {len(orf_records)} ORFs to '{out_faa}'.")
    return out_faa
import os
from Bio import SeqIO

from deneme1_1 import (
    build_tagged_refs,
    load_contig_lengths,
    parse_orf_coords,
    predict_bacterial_orfs,
    run_methyl_cycle_blastp,
    summarize_luxs_neighborhood,
)
from luxs_neighborhood_plot import plot_luxs_neighborhood_map


def run_gram_positive_workflow(
    genome_fa: str = "bacterial_genome.fasta",
    orfs_faa: str = "predicted_orfs_clean.faa",
    methyl_ref_fa: str = "gram_pos_methyl_cycle_ref.fasta",
    window_bp: int = 10000,
    force_refresh: bool = False,
) -> None:
    contig_lengths = load_contig_lengths(genome_fa)

    if force_refresh or not os.path.exists(orfs_faa):
        predict_bacterial_orfs(
            genome_fasta=genome_fa,
            out_faa=orfs_faa,
            min_aa_len=80,
            table_id=11,
            fmet_start=True,
        )

    orf_records = list(SeqIO.parse(orfs_faa, "fasta"))
    all_orfs_df = parse_orf_coords(orf_records, contig_lengths)
    if all_orfs_df.empty:
        print("No valid ORF coordinates parsed; check FASTA headers and contig names.")
        return

    gp_taxon = (
        "(Bacillus[Organism] OR Staphylococcus[Organism] OR Streptococcus[Organism] "
        "OR Listeria[Organism] OR Enterococcus[Organism])"
    )
    methyl_gene_map = {
        "luxS": "luxS",
        "metK": "metK",
        "mtnN": "mtnN",
        "metE": "metE",
        "metH": "metH",
    }

    build_tagged_refs(
        jobs=[(gp_taxon, methyl_gene_map)],
        out_fasta=methyl_ref_fa,
        retmax_per_gene=100,
        force_download=force_refresh,
    )

    methyl_hits = run_methyl_cycle_blastp(
        query_faa=orfs_faa,
        ref_fasta=methyl_ref_fa,
        blast_tsv="gram_pos_methyl_cycle_blast.tsv",
        out_csv="gram_pos_methyl_cycle_summary.csv",
    )

    if not methyl_hits.empty:
        summarize_luxs_neighborhood(
            all_orfs_df=all_orfs_df,
            methyl_hits_df=methyl_hits,
            window_bp=window_bp,
            out_csv="gram_pos_luxs_neighborhood.csv",
        )
        plot_luxs_neighborhood_map(
            all_orfs_df=all_orfs_df,
            methyl_hits_df=methyl_hits,
            window_bp=window_bp,
            out_png="gram_pos_luxs_neighborhood_map.png",
            out_pdf="gram_pos_luxs_neighborhood_map.pdf",
        )


if __name__ == "__main__":
    run_gram_positive_workflow(
        genome_fa="bacterial_genome.fasta",
        orfs_faa="predicted_orfs_clean.faa",
        methyl_ref_fa="gram_pos_methyl_cycle_ref.fasta",
        window_bp=10000,
        force_refresh=False,
    )

# ==============================================================================
# 2. GRAM-NEGATIVE LSR / AI-2 REFERENCE DB WITH SUGAR & KINASE DECOYS
# ==============================================================================
def _entrez_fetch(query: str, out_fasta: str, retmax: int = 100) -> int:
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urlencode({
        "db": "protein",
        "term": query,
        "usehistory": "y",
        "retmax": retmax,
        "retmode": "json",
    })
    with urlopen(esearch_url, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))

    esr = data.get("esearchresult", {})
    if int(esr.get("count", 0)) == 0 or not esr.get("webenv") or not esr.get("querykey"):
        raise RuntimeError("NCBI esearch returned no records.")

    efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode({
        "db": "protein",
        "query_key": esr["querykey"],
        "WebEnv": esr["webenv"],
        "retmax": retmax,
        "rettype": "fasta",
        "retmode": "text",
    })
    with urlopen(efetch_url, timeout=60) as resp, open(out_fasta, "wb") as fh:
        fh.write(resp.read())
    return sum(1 for _ in SeqIO.parse(out_fasta, "fasta"))


def tag_from_title(stitle: str) -> str:
    m = re.search(r"gene=(\S+)", str(stitle))
    return m.group(1) if m else "?"


def _cache_ok(out_fasta: str, jobs: list) -> bool:
    if not (os.path.exists(out_fasta) and os.path.getsize(out_fasta) > 0):
        return False
    expected = {tag for _, gm in jobs for tag in gm.values()}
    found = set()
    for r in SeqIO.parse(out_fasta, "fasta"):
        tag = tag_from_title(r.description)
        if tag == "?":
            return False
        found.add(tag)
    return not (expected - found)


def build_tagged_refs(jobs: list, out_fasta: str, retmax_per_gene: int = 100, force_download: bool = False) -> int:
    if not force_download and _cache_ok(out_fasta, jobs):
        count = sum(1 for _ in SeqIO.parse(out_fasta, "fasta"))
        print(f"Using validated cache '{out_fasta}' ({count} records).")
        return count

    records = []
    seen = set()
    for taxon_q, gene_map in jobs:
        for gene, tag in gene_map.items():
            tmp = f"_tmp_gn_{gene}_{os.getpid()}.fasta"
            try:
                _entrez_fetch(f"{taxon_q} AND {gene}[Gene] AND refseq[filter]", tmp, retmax=retmax_per_gene)
                for r in SeqIO.parse(tmp, "fasta"):
                    if r.id in seen:
                        continue
                    seen.add(r.id)
                    rest = r.description[len(r.id) :].strip()
                    r.description = f"{r.id} gene={tag} {rest}".strip()
                    records.append(r)
            except (RuntimeError, URLError, TimeoutError) as e:
                print(f"Warning: failed to fetch '{gene}': {e}")
            finally:
                if os.path.exists(tmp):
                    os.remove(tmp)
                time.sleep(0.4)

    SeqIO.write(records, out_fasta, "fasta")
    print(f"Wrote {len(records)} tagged reference sequences to '{out_fasta}'.")
    return len(records)


def run_gram_neg_blastp(
    query_faa: str = "predicted_orfs_clean.faa",
    ref_fasta: str = "gram_neg_decoys_ref.fasta",
    blast_tsv: str = "gram_neg_lsr_blast.tsv",
    out_csv: str = "gram_neg_lsr_annotations.csv",
    min_pident: float = 30.0,
    min_aln_len: int = 40,
    min_scov: float = 0.50,
) -> pd.DataFrame:
    subprocess.run(["makeblastdb", "-in", ref_fasta, "-dbtype", "prot"], check=True)
    subprocess.run(
        [
            "blastp",
            "-query", query_faa,
            "-db", ref_fasta,
            "-out", blast_tsv,
            "-evalue", "1e-5",
            "-max_target_seqs", "5",
            "-outfmt", "6 qseqid sseqid pident length evalue bitscore qlen slen stitle",
        ],
        check=True,
    )

    cols = [
        "ORF_ID", "Subject_ID", "Percent_Identity",
        "Alignment_Length", "E_value", "Bitscore",
        "Query_Len", "Subject_Len", "Annotation",
    ]
    if not os.path.exists(blast_tsv) or os.path.getsize(blast_tsv) == 0:
        return pd.DataFrame(columns=cols + ["Subject_Cov", "Gene_Symbol"])

    df = pd.read_csv(blast_tsv, sep="\t", names=cols)
    df["Subject_Cov"] = df["Alignment_Length"] / df["Subject_Len"].replace(0, float("nan"))

    hits = df[
        (df["Percent_Identity"] >= min_pident)
        & (df["Alignment_Length"] >= min_aln_len)
        & (df["Subject_Cov"] >= min_scov)
    ].copy()

    hits = (
        hits.sort_values(["Bitscore", "E_value"], ascending=[False, True])
        .drop_duplicates(subset=["ORF_ID"], keep="first")
        .reset_index(drop=True)
    )
    hits["Gene_Symbol"] = hits["Annotation"].apply(tag_from_title)
    hits = hits[~hits["Gene_Symbol"].isin(["DECOY", "?"])].sort_values(["Gene_Symbol", "E_value"]).reset_index(drop=True)
    hits.to_csv(out_csv, index=False)
    print(f"Saved {len(hits)} filtered Gram-negative AI-2 / lsr homologs to '{out_csv}'.")
    return hits


def annotate_lsr_synteny(coords_df: pd.DataFrame, max_gap_bp: int = 5000, min_cluster_genes: int = 2) -> pd.DataFrame:
    if coords_df.empty:
        return coords_df
    df = coords_df.drop_duplicates(subset=["ORF_ID"]).sort_values(["Contig", "Genomic_Start", "Genomic_End"]).reset_index(drop=True)
    cids = [1]
    cur_cid = 1
    cur_contig = df.loc[0, "Contig"]
    cur_max_end = int(df.loc[0, "Genomic_End"])

    for i in range(1, len(df)):
        r = df.loc[i]
        if r["Contig"] == cur_contig and (int(r["Genomic_Start"]) - cur_max_end) <= max_gap_bp:
            cids.append(cur_cid)
            cur_max_end = max(cur_max_end, int(r["Genomic_End"]))
        else:
            cur_cid += 1
            cur_contig = r["Contig"]
            cur_max_end = int(r["Genomic_End"])
            cids.append(cur_cid)

    df["Cluster_ID"] = cids
    lsr_genes = {"lsrA", "lsrB", "lsrC", "lsrD", "lsrF", "lsrG", "lsrR", "lsrK"}
    counts = df[df["Gene_Symbol"].isin(lsr_genes)].groupby("Cluster_ID")["Gene_Symbol"].nunique().to_dict()

    statuses = []
    for _, r in df.iterrows():
        sym = r["Gene_Symbol"]
        if sym in lsr_genes and counts.get(r["Cluster_ID"], 0) < min_cluster_genes:
            statuses.append("Unconfirmed_Synteny")
        elif sym in lsr_genes:
            statuses.append("Confirmed_Cluster")
        else:
            statuses.append("Non_Lsr_Target")
    df["Synteny_Status"] = statuses
    return df


def compute_intergenic_distances(coords_df: pd.DataFrame, max_gap_bp: int = 5000, out_csv: str = "lsr_intergenic_distances.csv") -> pd.DataFrame:
    ordered = coords_df.drop_duplicates(subset=["ORF_ID"]).sort_values(["Contig", "Genomic_Start", "Genomic_End"]).reset_index(drop=True)
    rows = []
    for i in range(len(ordered) - 1):
        g1, g2 = ordered.iloc[i], ordered.iloc[i + 1]
        if g1["Contig"] != g2["Contig"] or (g2["Genomic_Start"] - g1["Genomic_End"]) > max_gap_bp:
            continue
        ig_bp = int(g2["Genomic_Start"] - g1["Genomic_End"] - 1)
        if g1["Strand"] == "fwd" and g2["Strand"] == "fwd":
            orient = "Tandem (-> ->)"
        elif g1["Strand"] == "rev" and g2["Strand"] == "rev":
            orient = "Tandem (<- <-)"
        elif g1["Strand"] == "fwd" and g2["Strand"] == "rev":
            orient = "Convergent (-> <-)"
        else:
            orient = "Divergent (<- ->)"
        rows.append({
            "Contig": g1["Contig"],
            "Gene_1": g1["ORF_ID"],
            "Gene_1_Sym": g1["Gene_Symbol"],
            "Gene_2": g2["ORF_ID"],
            "Gene_2_Sym": g2["Gene_Symbol"],
            "Intergenic_Bp": ig_bp,
            "Orientation": orient,
        })
    pairs_df = pd.DataFrame(rows)
    if not pairs_df.empty:
        pairs_df.to_csv(out_csv, index=False)
    return pairs_df


def scan_candidate_lsr_motifs(
    coords_df: pd.DataFrame,
    genome_seqs: dict,
    contig_lengths: dict,
    upstream_bp: int = 300,
    out_csv: str = "lsr_candidate_motifs.csv",
) -> pd.DataFrame:
    """
    Heuristic candidate scan (AATGAA.{6,8}TTCATT) on oriented upstream region only
    to avoid duplicate counts from palindromic symmetry.
    """
    unique_df = coords_df.drop_duplicates(subset=["ORF_ID"])
    targets = unique_df[unique_df["Gene_Symbol"] == "lsrA"]
    if targets.empty:
        targets = unique_df.sort_values(["Contig", "Genomic_Start"])

    regex = re.compile(r"AATGAA.{6,8}TTCATT")
    hits = []
    for _, row in targets.iterrows():
        contig = row["Contig"]
        seq_str = genome_seqs.get(contig)
        c_len = contig_lengths.get(contig)
        if not seq_str or not c_len:
            continue

        if row["Strand"] == "fwd":
            up_end = row["Genomic_Start"] - 1
            up_start = max(1, up_end - upstream_bp + 1)
            sub_seq = seq_str[up_start - 1 : up_end]
        else:
            up_start = row["Genomic_End"] + 1
            up_end = min(c_len, up_start + upstream_bp - 1)
            sub_seq = str(Seq(seq_str[up_start - 1 : up_end]).reverse_complement())

        for m in regex.finditer(sub_seq):
            hits.append({
                "Contig": contig,
                "Gene_ID": row["ORF_ID"],
                "Gene_Symbol": row["Gene_Symbol"],
                "Evidence": "Heuristic_Candidate",
                "Start_In_Upstream": m.start() + 1,
                "End_In_Upstream": m.end(),
                "Matched_Seq": m.group(0),
            })

    out_df = pd.DataFrame(hits)
    if not out_df.empty:
        out_df.to_csv(out_csv, index=False)
    print(f"Found {len(out_df)} heuristic candidate LsrR-like motif hits.")
    return out_df


def plot_lsr_operon_and_spacings(
    all_orfs_df: pd.DataFrame,
    coords_df: pd.DataFrame,
    pairs_df: pd.DataFrame,
    flank_bp: int = 2000,
    gap_merge_bp: int = 5000,
):
    if coords_df.empty:
        return
    ordered = coords_df.drop_duplicates(subset=["ORF_ID"]).sort_values(["Contig", "Genomic_Start", "Genomic_End"]).reset_index(drop=True)
    clusters = []
    cur_rows = [ordered.iloc[0]]
    cur_max_end = int(ordered.iloc[0]["Genomic_End"])
    cur_contig = ordered.iloc[0]["Contig"]

    for i in range(1, len(ordered)):
        r = ordered.iloc[i]
        if r["Contig"] == cur_contig and (int(r["Genomic_Start"]) - cur_max_end) <= gap_merge_bp:
            cur_rows.append(r)
            cur_max_end = max(cur_max_end, int(r["Genomic_End"]))
        else:
            clusters.append(pd.DataFrame(cur_rows))
            cur_rows = [r]
            cur_max_end = int(r["Genomic_End"])
            cur_contig = r["Contig"]
    clusters.append(pd.DataFrame(cur_rows))

    sym_map = dict(zip(ordered["ORF_ID"], ordered["Gene_Symbol"]))
    syn_map = dict(zip(ordered["ORF_ID"], ordered.get("Synteny_Status", pd.Series(index=ordered.index, dtype=str)).fillna("")))

    def role_color(sym: str) -> str:
        if sym == "luxS":
            return "#FDAE61"
        if sym in {"lsrA", "lsrB", "lsrC", "lsrD", "lsrF", "lsrG"}:
            return "#2B83BA"
        if sym in {"lsrR", "lsrK"}:
            return "#9E9AC8"
        return "#E0E0E0"

    fig, axes = plt.subplots(len(clusters), 1, figsize=(12, max(3.0, 2.2 * len(clusters))), squeeze=False)
    for idx, cl_df in enumerate(clusters):
        ax = axes[idx, 0]
        contig = cl_df["Contig"].iloc[0]
        w_start = max(1, int(cl_df["Genomic_Start"].min()) - flank_bp)
        w_end = int(cl_df["Genomic_End"].max()) + flank_bp
        span = max(1, w_end - w_start)

        sub = all_orfs_df[
            (all_orfs_df["Contig"] == contig)
            & (all_orfs_df["Genomic_End"] >= w_start)
            & (all_orfs_df["Genomic_Start"] <= w_end)
        ].drop_duplicates(subset=["ORF_ID"])

        for _, gene in sub.iterrows():
            x0, x1 = float(gene["Genomic_Start"]), float(gene["Genomic_End"])
            width = max(1.0, x1 - x0)
            sym = sym_map.get(gene["ORF_ID"], "?")
            syn = syn_map.get(gene["ORF_ID"], "")
            head = min(width * 0.25, span * 0.015)

            if gene["Forward"]:
                be = max(x0, x1 - head)
                pts = [(x0, 0.38), (be, 0.38), (x1, 0.50), (be, 0.62), (x0, 0.62)]
            else:
                bs = min(x1, x0 + head)
                pts = [(x1, 0.38), (bs, 0.38), (x0, 0.50), (bs, 0.62), (x1, 0.62)]

            ax.add_patch(
                Polygon(
                    pts,
                    closed=True,
                    facecolor=role_color(sym),
                    edgecolor="black",
                    linewidth=0.7,
                    linestyle="--" if syn == "Unconfirmed_Synteny" else "-",
                )
            )
            if sym != "?":
                suffix = "*" if syn == "Unconfirmed_Synteny" else ""
                ax.text((x0 + x1) / 2.0, 0.72, f"{sym}{suffix}", ha="center", va="bottom", fontsize=8.5, fontweight="bold")

        ax.set_xlim(w_start, w_end)
        ax.set_ylim(0.15, 1.05)
        ax.set_yticks([])
        ax.set_xlabel("Genomic Coordinate (bp)")
        ax.set_title(f"Gram-Negative lsr / AI-2 Cluster {idx + 1} ({contig})", loc="left", fontsize=10, fontweight="bold")
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{int(x):,}"))
        for sp in ["top", "right", "left"]:
            ax.spines[sp].set_visible(False)

    plt.tight_layout()
    plt.savefig("gram_neg_lsr_operon_map.png", dpi=300)
    plt.close(fig)

    if not pairs_df.empty:
        plot_pairs = pairs_df.iloc[::-1].copy()
        plot_pairs["Label"] = plot_pairs["Gene_1_Sym"] + " -> " + plot_pairs["Gene_2_Sym"]
        colors = ["#E41A1C" if bp < 0 else ("#377EB8" if bp <= 30 else "#4DAF4A") for bp in plot_pairs["Intergenic_Bp"]]

        fig2, ax2 = plt.subplots(figsize=(10, max(3.5, 0.45 * len(plot_pairs))))
        y_pos = np.arange(len(plot_pairs))
        ax2.barh(y_pos, plot_pairs["Intergenic_Bp"], color=colors, edgecolor="black", linewidth=0.5, height=0.65)
        ax2.axvline(0, color="black", linewidth=0.8)
        for y, bp in zip(y_pos, plot_pairs["Intergenic_Bp"]):
            ax2.text(bp + (-5 if bp < 0 else 5), y, f"{bp} bp", va="center", ha="right" if bp < 0 else "left", fontsize=8.5, fontweight="bold")
        ax2.set_yticks(y_pos)
        ax2.set_yticklabels(plot_pairs["Label"], fontsize=9, fontweight="bold")
        ax2.set_xlabel("Intergenic Distance (bp)")
        ax2.set_title("Adjacent lsr / AI-2 Gene Spacings", fontsize=12, fontweight="bold")
        for sp in ["top", "right"]:
            ax2.spines[sp].set_visible(False)
        plt.tight_layout()
        plt.savefig("gram_neg_lsr_spacings_bar.png", dpi=300)
        plt.close(fig2)


if __name__ == "__main__":
    genome_fa = "bacterial_genome.fasta"
    orfs_faa = "predicted_orfs_clean.faa"
    ref_fa = "gram_neg_ref.fasta"
    force_refresh = False

    contig_lengths = load_contig_lengths(genome_fa)
    genome_seqs = load_genome_seqs(genome_fa)
    ensure_orfs_faa(genome_fa, out_faa=orfs_faa, force_refresh=force_refresh)

    orf_records = list(SeqIO.parse(orfs_faa, "fasta"))
    all_orfs_df = parse_orf_coords(orf_records, contig_lengths)

    k12_taxon = '"Escherichia coli str. K-12 substr. MG1655"[Organism]'
    targets = {g: g for g in ["luxS", "lsrA", "lsrB", "lsrC", "lsrD", "lsrF", "lsrG", "lsrR", "lsrK"]}
    decoys = {
        g: "DECOY"
        for g in ["rbsA", "rbsB", "rbsC", "mglA", "mglB", "mglC", "xylF", "xylG", "xylH", "xylB", "glpK", "rhaB", "araB"]
    }
    gene_map = {**targets, **decoys}

    build_tagged_refs(jobs=[(k12_taxon, gene_map)], out_fasta=ref_fa, retmax_per_gene=100, force_download=force_refresh)
    hits_df = run_gram_neg_blastp(query_faa=orfs_faa, ref_fasta=ref_fa)

    if not hits_df.empty:
        coords_df = all_orfs_df[all_orfs_df["ORF_ID"].isin(hits_df["ORF_ID"])].merge(
            hits_df[["ORF_ID", "Gene_Symbol", "Percent_Identity", "Subject_Cov", "E_value"]],
            on="ORF_ID",
            how="left",
        )
        if not coords_df.empty:
            coords_df = annotate_lsr_synteny(coords_df, max_gap_bp=5000, min_cluster_genes=2)
            pairs_df = compute_intergenic_distances(coords_df)
            scan_candidate_lsr_motifs(coords_df, genome_seqs, contig_lengths)
            plot_lsr_operon_and_spacings(all_orfs_df, coords_df, pairs_df)
          import os
from Bio import SeqIO

from deneme1_1 import (
    annotate_synteny_status,
    build_tagged_refs,
    compute_lsr_intergenic_distances,
    load_contig_lengths,
    load_genome_seqs,
    parse_orf_coords,
    plot_lsr_operon_and_spacings,
    predict_bacterial_orfs,
    run_ai2_blastp,
    scan_lsr_upstream_motifs,
)


def run_gram_negative_workflow(
    genome_fa: str = "bacterial_genome.fasta",
    orfs_faa: str = "predicted_orfs_clean.faa",
    ai2_ref_fa: str = "gram_neg_ai2_k12_decoys_ref.fasta",
    force_refresh: bool = False,
) -> None:
    contig_lengths = load_contig_lengths(genome_fa)
    genome_seqs = load_genome_seqs(genome_fa)

    if force_refresh or not os.path.exists(orfs_faa):
        predict_bacterial_orfs(
            genome_fasta=genome_fa,
            out_faa=orfs_faa,
            min_aa_len=80,
            table_id=11,
            fmet_start=True,
        )

    orf_records = list(SeqIO.parse(orfs_faa, "fasta"))
    all_orfs_df = parse_orf_coords(orf_records, contig_lengths)
    if all_orfs_df.empty:
        print("No valid ORF coordinates parsed; check FASTA headers and contig names.")
        return

    k12_taxon = '"Escherichia coli str. K-12 substr. MG1655"[Organism]'
    ai2_gene_map = {
        g: g
        for g in ["luxS", "lsrA", "lsrB", "lsrC", "lsrD", "lsrF", "lsrG", "lsrR", "lsrK"]
    }
    ai2_gene_map.update({
        g: "DECOY"
        for g in [
            "rbsA", "rbsB", "rbsC",
            "mglA", "mglB", "mglC",
            "xylF", "xylG", "xylH",
            "xylB", "glpK", "rhaB", "araB",
        ]
    })

    build_tagged_refs(
        jobs=[(k12_taxon, ai2_gene_map)],
        out_fasta=ai2_ref_fa,
        retmax_per_gene=100,
        force_download=force_refresh,
    )

    ai2_hits = run_ai2_blastp(
        query_faa=orfs_faa,
        ref_fasta=ai2_ref_fa,
        blast_tsv="gram_neg_ai2_blast.tsv",
        out_csv="gram_neg_ai2_annotations.csv",
        extracted_faa="gram_neg_extracted_ai2_proteins.faa",
    )

    if ai2_hits.empty:
        print("No filtered AI-2 / lsr homologs found.")
        return

    ai2_coords_df = all_orfs_df[all_orfs_df["ORF_ID"].isin(ai2_hits["ORF_ID"])].merge(
        ai2_hits[["ORF_ID", "Gene_Symbol", "Annotation", "Percent_Identity", "Subject_Cov", "E_value"]],
        on="ORF_ID",
        how="left",
    )
    if ai2_coords_df.empty:
        print("Matched ORF_IDs did not align with all_orfs_df coordinates.")
        return

    ai2_coords_df = annotate_synteny_status(
        ai2_coords_df,
        max_gap_bp=5000,
        min_lsr_cluster_genes=2,
    )
    pairs_df = compute_lsr_intergenic_distances(
        ai2_coords_df,
        max_gap_bp=5000,
        out_csv="gram_neg_lsr_intergenic_distances.csv",
    )
    scan_lsr_upstream_motifs(
        ai2_coords_df,
        genome_seqs=genome_seqs,
        contig_lengths=contig_lengths,
        upstream_bp=300,
        out_promoters_csv="gram_neg_lsr_upstream_promoters.csv",
        out_motifs_csv="gram_neg_lsr_candidate_motifs.csv",
    )
    plot_lsr_operon_and_spacings(all_orfs_df, ai2_coords_df, pairs_df)


if __name__ == "__main__":
    run_gram_negative_workflow(
        genome_fa="bacterial_genome.fasta",
        orfs_faa="predicted_orfs_clean.faa",
        ai2_ref_fa="gram_neg_ai2_k12_decoys_ref.fasta",
        force_refresh=False,
    )
