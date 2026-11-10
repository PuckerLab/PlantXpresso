"""
search.py  --  "I have a sequence, which gene is it?"  (no web code in here)

What happens for one search:
  1. check_query()  : is it ONE sequence in FASTA format, and a valid CDS / protein?
  2. create_job()   : make a job folder with a random ID, save the sequence as query.fasta
  3. run_job()      : wait for a free "slot" (max MAX_SEARCHES searches at once on the VM),
                      run DIAMOND  (blastx for a CDS, blastp for a protein) against the
                      species' protein database pep.dmnd, rank hits by bit score
  4. read_job()     : the website reads job.json to show "queued" / "running" / results

Why a job folder and a job.json file?
  The website runs as several worker processes. A search started by worker 1 must be
  visible when worker 3 answers the "is it done?" request. Files are the one thing all
  workers share, so the job's state lives in job.json, not in Python memory.
  The folder also keeps every search separate: two users never share a file.
"""

import fcntl
import json
import os
import shutil
import subprocess
import threading
import time
import uuid
from pathlib import Path

from .expression import cache_dir

# ---- Settings ----------------------------------------------------------------------
JOBS_DIR = Path(os.environ.get("PLANTX_JOBS_DIR", "/vol/data/plantxpresso_jobs"))
DIAMOND = os.environ.get("PLANTX_DIAMOND", "diamond")   # program name or full path
EVALUE = "1e-5"                  # hits with a higher e-value are ignored
MAX_HITS = 10                    # hits shown in the table
MAX_SEARCHES = int(os.environ.get("PLANTX_MAX_SEARCHES", "3"))   # at once, on the whole VM
THREADS = int(os.environ.get("PLANTX_SEARCH_THREADS", "4"))      # CPU cores per search
KEEP_DAYS = 7                    # job folders older than this are deleted (by build_cache.py)
STALE_MINUTES = 30               # a job "running" longer than this was interrupted

NUCLEOTIDES = set("ACGTN")
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY") | {"X"}   # 20 standard + X (unknown)
STOP_CODONS = {"TAA", "TAG", "TGA"}

DATABASE_NAME = "pep.dmnd"       # built by build_cache.py from <Species>.pep.fasta


def database(sp):
    return cache_dir(sp) / DATABASE_NAME


# ================================================================ validation ========
def check_query(text, seq_type):
    """
    Returns (header, sequence, errors, warnings).
    errors   -> the search is refused, the user sees why
    warnings -> the search runs, the user sees a note
    """
    errors, warnings = [], []
    text = (text or "").strip()

    # --- FASTA format: one '>' header line, then sequence lines ---
    if not text.startswith(">"):
        return "", "", ["The input must be in FASTA format: a first line starting with '>', "
                        "then the sequence."], []
    lines = text.splitlines()
    if any(line.startswith(">") for line in lines[1:]):
        return "", "", ["Only one sequence can be searched at a time."], []
    header = lines[0][1:].strip() or "query"
    seq = "".join(lines[1:]).replace(" ", "").replace("\t", "").upper()
    if not seq:
        return header, "", ["The FASTA entry has a header but no sequence."], []

    # --- CDS (nucleotide) checks ---
    if seq_type == "cds":
        bad = sorted(set(seq) - NUCLEOTIDES)
        if bad:
            errors.append(f"A CDS may only contain A, C, G, T or N; found: {' '.join(bad)}."
                          + (" Is this a protein sequence?" if set(seq) <= AMINO_ACIDS | {"*"} else ""))
        elif not 30 <= len(seq) <= 50000:
            errors.append(f"A CDS must be 30-50,000 nucleotides long (this one: {len(seq):,}).")
        elif len(seq) % 3:
            errors.append(f"A CDS length must be a multiple of 3 (this one: {len(seq):,} nt).")
        else:
            codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
            internal = sum(c in STOP_CODONS for c in codons[:-1])
            if internal:
                errors.append(f"{internal} stop codon(s) inside the reading frame: not a valid CDS "
                              "(does it contain UTRs or introns, or is it reverse-complemented?).")
            if codons[0] != "ATG":
                warnings.append("Does not start with ATG (partial CDS?).")
            if codons[-1] not in STOP_CODONS:
                warnings.append("Does not end with a stop codon (partial CDS?).")

    # --- protein checks ---
    elif seq_type == "protein":
        if seq.endswith("*"):
            seq = seq[:-1]                      # a final stop symbol is fine
        bad = sorted(set(seq) - AMINO_ACIDS)
        if "*" in seq:
            errors.append("The protein contains an internal stop ('*').")
        elif bad:
            errors.append(f"Not valid amino-acid letters: {' '.join(bad)}.")
        elif set(seq) <= NUCLEOTIDES:
            errors.append("This looks like a nucleotide sequence. Choose 'CDS' instead.")
        elif not 10 <= len(seq) <= 10000:
            errors.append(f"A protein must be 10-10,000 amino acids long (this one: {len(seq):,}).")
    else:
        errors.append("Choose whether the sequence is a CDS or a protein.")
    return header, seq, errors, warnings


# ===================================================================== jobs ========
def _write_job(folder, job):
    """Write job.json in one step (temp file + rename), so a reader never sees half a file."""
    tmp = folder / "job.json.tmp"
    tmp.write_text(json.dumps(job, indent=1))
    tmp.rename(folder / "job.json")


def create_job(sp, header, seq, seq_type, warnings):
    """Make the job folder, save the query, return the job ID."""
    job_id = uuid.uuid4().hex            # 32 random characters: impossible to guess
    folder = JOBS_DIR / job_id
    folder.mkdir(parents=True)
    (folder / "query.fasta").write_text(f">{header}\n{seq}\n")
    _write_job(folder, {
        "id": job_id, "status": "queued", "species": sp.name, "version": sp.version,
        "seq_type": seq_type, "program": "blastx" if seq_type == "cds" else "blastp",
        "header": header, "length": len(seq), "warnings": warnings,
        "submitted": time.time(), "started": None, "finished": None,
        "error": None, "hits": [],
    })
    return job_id


def read_job(job_id):
    """The job as a dict, or None if it doesn't exist (wrong ID or already deleted)."""
    if not (len(job_id) == 32 and job_id.isalnum()):       # only well-formed IDs
        return None
    path = JOBS_DIR / job_id / "job.json"
    if not path.is_file():
        return None
    job = json.loads(path.read_text())
    # A job can only be stuck if the server restarted while it waited or ran
    if job["status"] in ("queued", "running") and time.time() - job["submitted"] > STALE_MINUTES * 60:
        job["status"], job["error"] = "failed", "The search was interrupted. Please submit it again."
    return job


def start_job(sp, job_id):
    """Run the search in the background, so the web page can answer immediately."""
    threading.Thread(target=run_job, args=(sp, job_id), daemon=True).start()


def _acquire_slot():
    """
    Wait until one of MAX_SEARCHES "slots" is free, and take it.
    A slot is a lock on a small file (.slot0.lock, .slot1.lock, ...). Locks are shared by
    all worker processes, and the operating system releases a lock automatically if a
    process dies, so a crash can never leave a slot blocked forever.
    """
    while True:
        for i in range(MAX_SEARCHES):
            fh = open(JOBS_DIR / f".slot{i}.lock", "w")
            try:
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)   # take it, don't wait
                return fh                                        # keep it open = keep the slot
            except BlockingIOError:
                fh.close()                                       # taken, try the next one
        time.sleep(1)                                            # all taken: wait, try again


def run_job(sp, job_id):
    folder = JOBS_DIR / job_id
    job = json.loads((folder / "job.json").read_text())
    slot = _acquire_slot()
    try:
        job["status"], job["started"] = "running", time.time()
        _write_job(folder, job)

        cmd = [DIAMOND, job["program"],
               "--db", str(database(sp)),
               "--query", str(folder / "query.fasta"),
               "--out", str(folder / "hits.tsv"),
               "--outfmt", "6", "sseqid", "bitscore", "evalue", "pident", "length", "qcovhsp",
               "--evalue", EVALUE,
               "--max-target-seqs", str(MAX_HITS),
               "--threads", str(THREADS),
               "--quiet"]
        # the user's sequence is only ever inside query.fasta, never part of the command
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if proc.returncode != 0:
            raise RuntimeError(proc.stderr.strip()[-500:] or "DIAMOND failed")

        job["hits"] = rank_hits(folder / "hits.tsv")
        job["status"] = "done"
    except Exception as e:
        job["status"], job["error"] = "failed", f"Search failed: {e}"
    finally:
        job["finished"] = time.time()
        _write_job(folder, job)
        slot.close()                                             # free the slot


def rank_hits(path):
    """
    DIAMOND's tabular output -> list of hits, best first.
    One target can have several aligned pieces (HSPs); keep its best one.
    """
    best = {}
    for line in path.read_text().splitlines():
        sseqid, bitscore, evalue, pident, length, qcov = line.split("\t")
        hit = {"id": sseqid, "bitscore": float(bitscore), "evalue": float(evalue),
               "identity": float(pident), "length": int(length), "query_coverage": float(qcov)}
        if sseqid not in best or hit["bitscore"] > best[sseqid]["bitscore"]:
            best[sseqid] = hit
    return sorted(best.values(), key=lambda h: h["bitscore"], reverse=True)[:MAX_HITS]


def cleanup_old_jobs():
    """Delete job folders older than KEEP_DAYS (called nightly by build_cache.py)."""
    if not JOBS_DIR.is_dir():
        return 0
    removed = 0
    for folder in JOBS_DIR.iterdir():
        if folder.is_dir() and time.time() - folder.stat().st_mtime > KEEP_DAYS * 86400:
            shutil.rmtree(folder, ignore_errors=True)
            removed += 1
    return removed
