"""
expression.py  --  everything about gene expression and PCA, with no web code in it.
(Like catalog.py: you can import and test it from a plain Python shell.)

THE SPEED PROBLEM AND ITS SOLUTION
  A .tpms.tsv is text: to find one gene, a program must read and parse the file from the
  top until it reaches that gene (seconds per GB, for every single query).
  So build_cache.py converts each species ONCE into a binary matrix (expr.npy):
      row i = gene i, column j = sample j, every value a 4-byte float.
  Because every row has the same size, row i starts at byte  i * n_samples * 4.
  numpy's "memory map" (mmap) opens the file without reading it and fetches only the bytes
  of the one row we ask for: a few KB instead of the whole file -> milliseconds.
  The operating system keeps recently read parts in RAM (the "buff/cache" in `free -h`),
  and that RAM is shared by all uvicorn workers, so nothing is loaded 8 times.

  PCA uses ALL genes of a species, so it can't be made fast per query. Instead
  build_cache.py computes it once per species version and stores the coordinates (pca.tsv);
  the website only draws them.

FILES PER SPECIES VERSION in CACHE_DIR/<Species>/<version>/  (written by build_cache.py)
  expr.npy         genes x samples matrix of TPM values (float32)
  genes.txt        gene IDs, one per line, in row order of expr.npy
  samples.txt      run accessions, one per line, in column order of expr.npy
  pca.tsv          sample_id, tissue, PC1, PC2, ...          (missing if PCA impossible)
  pca_variance.txt explained variance ratio of each PC, one per line
  sources.txt      which input files (and their timestamps) this cache was built from
"""

import csv
import io
import os
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import matplotlib
matplotlib.use("Agg")                 # draw into files/memory, never open a window
from matplotlib.figure import Figure  # Figure (not pyplot): safe when several requests draw at once
import colorcet as cc                 # Glasbey colour palettes for many categories
import numpy as np

# ---- Settings (override with environment variables) --------------------------------
CACHE_DIR = Path(os.environ.get("PLANTX_CACHE_DIR", "/vol/data/plantxpresso_cache"))
TISSUE_DIR = Path(os.environ.get("PLANTX_TISSUE_DIR", "/vol/data/datasets/PlantXpresso_tissue_metadata"))
TISSUE_FILE = "Specific_class_labelled_metadata.tsv"


# ======================================================================= paths ======
def cache_dir(sp):
    """/vol/data/plantxpresso_cache/Lotus_japonicus/A1E0"""
    return CACHE_DIR / sp.name / sp.version


def tissue_metadata_path(sp):
    """Tissue labels for the SAME version as the expression data."""
    return TISSUE_DIR / sp.name / sp.version / TISSUE_FILE


def read_tissue_labels(path):
    """
    {run_accession: tissue} from the metadata file. Returns {} when the file is missing
    or has only the header row -> callers then fall back to the all-samples plot.
    """
    labels = {}
    if not path.is_file():
        return labels
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            run = (row.get("run_accession") or "").strip()
            tissue = (row.get("tissue") or "").strip()
            if run and tissue and tissue not in ("NA", "NaN"):
                labels[run] = tissue
    return labels


# ============================================================ expression data =======
@lru_cache(maxsize=64)
def _open_matrix(folder, stamp):
    """
    Open one species' cache. Kept in memory for the 64 most recently used species
    (lru_cache), so repeated queries skip even this step. `stamp` is the file's
    modification time: when build_cache.py rebuilds, the stamp changes and we reopen.
    """
    folder = Path(folder)
    genes = (folder / "genes.txt").read_text().splitlines()
    samples = (folder / "samples.txt").read_text().splitlines()
    matrix = np.load(folder / "expr.npy", mmap_mode="r")      # opens, does NOT read
    row_of = {gene: i for i, gene in enumerate(genes)}         # gene ID -> row number
    keys = [(normalize(gene), gene) for gene in genes]          # for searching, see below
    return row_of, samples, matrix, keys


def open_matrix(sp):
    """(row_of, samples, matrix, keys) for a species, or None if its cache isn't built yet."""
    npy = cache_dir(sp) / "expr.npy"
    if not npy.is_file():
        return None
    return _open_matrix(str(cache_dir(sp)), npy.stat().st_mtime)


@dataclass
class GeneResult:
    species: object        # the catalog's Species object
    gene: str
    samples: list          # all run accessions with a value
    values: np.ndarray     # TPM per sample (same order as samples)
    groups: list           # [(label, np.ndarray of TPMs), ...] -> one violin each
    tissue_wise: bool      # True = one violin per tissue, False = one violin of all samples


def gene_expression(sp, gene):
    """
    Look up one gene: read its row, attach tissue labels, group values for the plot.
    Raises LookupError with a user-readable message if something is missing.
    """
    opened = open_matrix(sp)
    if opened is None:
        raise LookupError(f"Expression data for {sp.display_name} is still being prepared. "
                          "Please try again later.")
    row_of, samples, matrix, keys = opened
    typed = gene
    if gene not in row_of:
        gene = resolve_gene(keys, gene)     # e.g. "AT5G08640.1" -> "rna-TAIR12-...-AT5G08640-1"
    if gene is None:
        raise LookupError(f"'{typed}' was not found in {sp.display_name}.")

    values = np.array(matrix[row_of[gene]], dtype=np.float64)   # <- the one small read
    keep = ~np.isnan(values)                                     # drop missing values
    samples = [s for s, k in zip(samples, keep) if k]
    values = values[keep]

    # One group per tissue, using only samples that have a tissue label
    labels = read_tissue_labels(tissue_metadata_path(sp))
    by_tissue = {}
    for s, v in zip(samples, values):
        if s in labels:
            by_tissue.setdefault(labels[s], []).append(v)

    if by_tissue:
        groups = [(t, np.array(v)) for t, v in by_tissue.items()]
        # same as your reorder_by_median(): highest median tissue first
        groups.sort(key=lambda g: np.median(g[1]), reverse=True)
        tissue_wise = True
    else:
        # no labels (header-only metadata) -> one violin of all samples
        groups = [("All samples", values)]
        tissue_wise = False
    return GeneResult(sp, gene, samples, values, groups, tissue_wise)


# ============================================================== ID matching =========
# Gene IDs look very different between annotations, e.g.
#     rna-TAIR12-TAIR12-AT5G08640-1      LOC_Os03g01234.1      Solyc01g005000.3.1
# Users type part of an ID, often with different punctuation ("AT5G08640.1").
# So IDs and the typed text are both "normalized" before comparing:
#     lower case, and only letters + digits kept
#     "rna-TAIR12-TAIR12-AT5G08640-1"  ->  "rnatair12tair12at5g086401"
#     "AT5G08640.1"                    ->  "at5g086401"      (found inside the first one)

def normalize(text):
    return re.sub(r"[^a-z0-9]", "", text.lower())


def match_ids(keys, query, limit=20):
    """
    IDs whose normalized form CONTAINS the normalized query (anywhere, not only at the start).
    Best first: exact match, then IDs starting with the query, then the rest;
    within each group shorter IDs first (they are usually the closer match).
    """
    q = normalize(query)
    if len(q) < 3:                       # too short: would match thousands of IDs
        return []
    found = []
    for key, gene in keys:               # ~50k comparisons: a few milliseconds
        position = key.find(q)
        if position >= 0:
            rank = 0 if key == q else 1 if position == 0 else 2
            found.append((rank, len(gene), gene))
    found.sort()
    return [gene for _, _, gene in found[:limit]]


def resolve_gene(keys, query):
    """
    The real ID for what the user typed, if it's unambiguous:
      1. exactly one ID that is the same after normalizing ("at1g01010.1" -> "AT1G01010-1"), or
      2. exactly one ID that contains it ("AT5G08640.1" -> "rna-TAIR12-TAIR12-AT5G08640-1").
    Otherwise None (the page then shows the suggestions to choose from).
    """
    q = normalize(query)
    if not q:
        return None
    exact = [gene for key, gene in keys if key == q]
    if len(exact) == 1:
        return exact[0]
    containing = [gene for key, gene in keys if q in key]
    return containing[0] if len(containing) == 1 and len(q) >= 3 else None


def suggest_genes(sp, query, limit=20):
    """Suggestions for the autocomplete list and for the 'not found' message."""
    opened = open_matrix(sp)
    if opened is None:
        return []
    return match_ids(opened[3], query, limit)


def summary(values):
    """min / max / mean / median of raw TPM values (not log-transformed)."""
    return {"n": len(values), "min": float(np.min(values)), "max": float(np.max(values)),
            "mean": float(np.mean(values)), "median": float(np.median(values))}


# ================================================================ violin plot =======
def violin_figure(result, scale_width=False):
    """
    Your plot_violin_with_box(), unchanged in look: green violins with black outline,
    white mini box plots, log10(TPM + 1), tissues ordered by median (done in gene_expression).
    """
    labels = [label for label, _ in result.groups]
    data = [np.log2(values + 1.0) for _, values in result.groups]
    sample_counts = [len(d) for d in data]
    positions = list(range(1, len(labels) + 1))

    fig = Figure(figsize=(max(8, len(labels) * 1.2), 6))
    ax = fig.subplots()

    if scale_width:   # wider violin = more samples (square-root scaled, min 35 % width)
        max_n = max(sample_counts)
        violin_widths = [0.8 * max((n / max_n) ** 0.5, 0.35) for n in sample_counts]
        box_widths = [w * (0.15 / 0.8) for w in violin_widths]
    else:
        violin_widths = [0.8] * len(labels)
        box_widths = [0.15] * len(labels)

    # A violin needs >= 2 different values; a tissue where the gene is e.g. 0 everywhere
    # gets only its box (a flat line), instead of an error.
    can_violin = [i for i, d in enumerate(data) if len(d) > 1 and np.ptp(d) > 0]
    if can_violin:
        violin_parts = ax.violinplot(
            [data[i] for i in can_violin], positions=[positions[i] for i in can_violin],
            widths=[violin_widths[i] for i in can_violin],
            showmeans=False, showmedians=False, showextrema=False,
        )
        for body in violin_parts["bodies"]:
            body.set_facecolor("green")
            body.set_edgecolor("black")
            body.set_linewidth(1.5)
            body.set_alpha(1.0)

    box_parts = ax.boxplot(data, positions=positions, widths=box_widths,
                           patch_artist=True, showfliers=False)
    for patch in box_parts["boxes"]:
        patch.set_facecolor("white")
        patch.set_edgecolor("black")
        patch.set_linewidth(1.5)
    for element in ("whiskers", "caps", "medians"):
        for artist in box_parts[element]:
            artist.set_color("black")
            artist.set_linewidth(1.5)

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{label}\n(n={n})" for label, n in zip(labels, sample_counts)],
                       rotation=45 if result.tissue_wise else 0, ha="right" if result.tissue_wise else "center")
    ax.set_xlabel("Tissue" if result.tissue_wise else "")
    ax.set_ylabel("Gene Expression (log2(TPM + 1))")
    where = "across tissues" if result.tissue_wise else "across all samples (no tissue labels)"
    ax.set_title(f"Expression of {result.gene} {where}\n{result.species.display_name}")
    fig.tight_layout()
    return fig


# =================================================================== PCA plot =======
def load_pca(sp):
    """(rows, variance) from the precomputed files, or None if not available."""
    folder = cache_dir(sp)
    if not (folder / "pca.tsv").is_file():
        return None
    with open(folder / "pca.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    variance = [float(x) for x in (folder / "pca_variance.txt").read_text().split()]
    return rows, variance


def pca_figure(sp, rows, variance):
    """Your do_pca.py plot: PC1 vs PC2, one colour per tissue, legend on the right."""
    tissues = list(dict.fromkeys(r["tissue"] for r in rows))      # unique, in file order
    # Glasbey palette: the first 10 colours are the usual tab10 ones, then up to 256
    # maximally different colours, so no two tissues share a colour
    palette = cc.glasbey_category10
    color_map = {t: tuple(palette[i % len(palette)]) for i, t in enumerate(tissues)}

    fig = Figure(figsize=(10, 7))
    ax = fig.subplots()
    for tissue in tissues:
        pts = [r for r in rows if r["tissue"] == tissue]
        ax.scatter([float(r["PC1"]) for r in pts], [float(r["PC2"]) for r in pts],
                   label=tissue, color=color_map[tissue],
                   s=100, alpha=0.7, edgecolor="black", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({variance[0] * 100:.2f}% variance)", fontsize=13)
    ax.set_ylabel(f"PC2 ({variance[1] * 100:.2f}% variance)", fontsize=13)
    ax.set_title(f"Global PCA of Gene Expression by Tissue Type\n{sp.display_name}",
                 fontsize=15, fontweight="bold", pad=20)
    ax.legend(title="Tissue", bbox_to_anchor=(1.05, 1), loc="upper left", frameon=True, fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig


# ===================================================================== export =======
def figure_bytes(fig, fmt, dpi):
    """Render a figure to PNG or SVG bytes (in memory, nothing written to disk)."""
    if fmt == "png":
        # keep very wide plots (many tissues) below ~20000 px so 600 dpi can't explode
        dpi = min(dpi, 20000 / fig.get_figwidth())
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight")
    return buf.getvalue()
