"""
build_cache.py  --  prepares the fast expression files and the PCA for every species.

Run it after new data arrives (cron does this, see INSTRUCTIONS), or by hand:
    cd ~/plantxpresso && ./venv/bin/python build_cache.py
    ./venv/bin/python build_cache.py Lotus_japonicus      # just one species

For each species in the catalog (its CURRENT version only):
  1. skip it if the cache was built from exactly the same input files (sources.txt)
  2. read <Species>.tpms.tsv (never the .repr.tpms.tsv)  -> expr.npy, genes.txt, samples.txt
  3. read the tissue labels and compute the PCA          -> pca.tsv, pca_variance.txt
  4. build the DIAMOND protein database from <Species>.pep.fasta (never .repr) -> pep.dmnd
  5. swap the new folder in place of the old one (the website never sees a half-built cache)
Afterwards, cache folders of old versions / removed species are deleted, and
sequence-search job folders older than 7 days are deleted.

Only reads from /vol/data/datasets; writes only to CACHE_DIR (and deletes old jobs).
"""

import shutil
import subprocess
import sys
import time

import numpy as np
import polars as pl
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from app.catalog import catalog
from app.expression import CACHE_DIR, cache_dir, read_tissue_labels, tissue_metadata_path
from app.search import DATABASE_NAME, DIAMOND, cleanup_old_jobs


def fingerprint(paths):
    """
    One line per input file: name, size, modification time ('missing' if absent).
    Only the file NAME is used (not the full path), so moving or renaming the data
    folder doesn't make every species look changed.
    """
    lines = []
    for p in paths:
        if p.is_file():
            st = p.stat()
            lines.append(f"{p.name}\t{st.st_size}\t{st.st_mtime}")
        else:
            lines.append(f"{p.name}\tmissing")
    return "\n".join(lines) + "\n"


def compute_pca(matrix, samples, labels, out):
    """
    Same steps as your do_pca.py:
      keep samples with a tissue label (all samples if there are no labels),
      keep genes with mean TPM >= 1, log2(TPM + 1), z-score per gene, PCA with 10 PCs.
    """
    labelled = [j for j, s in enumerate(samples) if s in labels]
    cols = labelled if labelled else list(range(len(samples)))
    X = np.nan_to_num(matrix[:, cols])                       # genes x chosen samples, missing -> 0
    X = X[X.mean(axis=1) >= 1]                               # expressed genes only
    X = np.log2(X + 1)
    if X.shape[0] < 2 or X.shape[1] < 3:
        return "PCA skipped (fewer than 3 samples or 2 expressed genes)"

    Z = StandardScaler().fit_transform(X.T)                  # samples x genes, z-scored
    n_pcs = min(10, Z.shape[0], Z.shape[1])
    pca = PCA(n_components=n_pcs, random_state=0)            # random_state: same result every run
    pcs = pca.fit_transform(Z)

    with open(out / "pca.tsv", "w") as fh:
        fh.write("sample_id\ttissue\t" + "\t".join(f"PC{k + 1}" for k in range(n_pcs)) + "\n")
        for row, j in zip(pcs, cols):
            tissue = labels.get(samples[j], "no tissue label")
            fh.write(f"{samples[j]}\t{tissue}\t" + "\t".join(f"{v:.6g}" for v in row) + "\n")
    (out / "pca_variance.txt").write_text("\n".join(f"{v:.6g}" for v in pca.explained_variance_ratio_) + "\n")
    return f"PCA on {len(cols)} samples x {X.shape[0]} genes"


def build_diamond_db(pep, out):
    """`diamond makedb`: turns the protein FASTA into DIAMOND's fast-search format."""
    if not pep.is_file():
        return "no .pep.fasta, no sequence search"
    if shutil.which(DIAMOND) is None:
        return "DIAMOND not installed, no sequence search"
    db = out / DATABASE_NAME
    subprocess.run([DIAMOND, "makedb", "--in", str(pep), "--db", str(db), "--quiet"],
                   check=True, capture_output=True)
    return "DIAMOND db built"


def build_species(sp):
    tpms = sp.folder / f"{sp.name}.tpms.tsv"                 # the full file, not .repr
    pep = sp.folder / f"{sp.name}.pep.fasta"                 # the full file, not .repr
    meta = tissue_metadata_path(sp)
    final = cache_dir(sp)
    if not tpms.is_file():
        return "no .tpms.tsv - skipped"

    sources = fingerprint([tpms, pep, meta])
    if (final / "sources.txt").is_file() and (final / "sources.txt").read_text() == sources:
        return "up to date"

    # Build into <version>.tmp, swap in at the end
    tmp = final.with_name(final.name + ".tmp")
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)

    df = pl.read_csv(tpms, separator="\t", null_values=["NA", "NaN", ""])
    genes = df[:, 0].cast(pl.Utf8).to_list()                 # first column = gene IDs
    samples = [c.strip() for c in df.columns[1:]]
    if not genes or not samples:
        shutil.rmtree(tmp)
        return "empty or malformed .tpms.tsv - skipped"
    matrix = df.select(pl.all().exclude(df.columns[0]).cast(pl.Float32)).to_numpy()
    del df                                                   # free the text-parsed copy

    np.save(tmp / "expr.npy", np.ascontiguousarray(matrix, dtype=np.float32))
    (tmp / "genes.txt").write_text("\n".join(genes) + "\n")
    (tmp / "samples.txt").write_text("\n".join(samples) + "\n")

    labels = read_tissue_labels(meta)
    pca_note = compute_pca(matrix, samples, labels, tmp)
    db_note = build_diamond_db(pep, tmp)
    (tmp / "sources.txt").write_text(sources)                # written last = "build complete"

    old = final.with_name(final.name + ".old")
    shutil.rmtree(old, ignore_errors=True)
    if final.exists():
        final.rename(old)
    tmp.rename(final)
    shutil.rmtree(old, ignore_errors=True)
    return f"built: {len(genes)} genes x {len(samples)} samples, {len(labels)} labelled; {pca_note}; {db_note}"


def remove_stale(species):
    """Delete cache folders of versions/species that are no longer current."""
    if not CACHE_DIR.is_dir():
        return
    for sp_dir in CACHE_DIR.iterdir():
        current = species.get(sp_dir.name)
        for v_dir in sp_dir.iterdir():
            if current is None or v_dir.name != current.version:
                shutil.rmtree(v_dir, ignore_errors=True)
                print(f"removed old cache {v_dir}")
        if not any(sp_dir.iterdir()):
            sp_dir.rmdir()


def main(only):
    species = catalog.get().species
    names = only or sorted(species)
    for name in names:
        if name not in species:
            print(f"{name}: not in catalog")
            continue
        start = time.time()
        try:
            msg = build_species(species[name])
        except Exception as e:                               # one bad file must not stop the rest
            msg = f"ERROR {type(e).__name__}: {e}"
        print(f"{name} {species[name].version}: {msg} ({time.time() - start:.1f}s)", flush=True)
    if not only:
        remove_stale(species)
        print(f"deleted {cleanup_old_jobs()} search job(s) older than 7 days")


if __name__ == "__main__":
    main(sys.argv[1:])
