"""
catalog.py  --  turns the metadata TSV + the dataset folders into a browsable catalog.

Ground rules (so future updates "just work"):
  * The metadata TSV is the source of truth for WHAT exists (species, order, family, ...).
  * The folders are the source of truth for WHICH VERSION is downloadable:
	only an *untarred* version folder (e.g. Lotus_japonicus/A1E0/) counts as "current".
  * If the TSV lists several versions of a species (Acer_saccharum_A1E0, Acer_saccharum_A2E0),
	we pick the highest version whose folder exists. So if the TSV already has A2E0 but rsync
	has not delivered the folder yet, users still get A1E0 -- nothing breaks.
  * Old versions that were tar.gz-ed (e.g. A1E0.tar.gz) are listed as "archived versions".
  * The catalog reloads itself when the TSV changes, or every RELOAD_SECONDS anyway
	(to catch folder changes from rsync). No server restart needed after an update.
  * Two metadata files: the internal one (METADATA_TSV) only builds the catalog;
	the webserver one (WEBSERVER_TSV) is what users see (class MetadataTable at the bottom).
"""

import csv
import os
import re
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

# ---- Settings (override with environment variables, see README) -------------------
DATA_DIR = Path(os.environ.get("PLANTX_DATA_DIR", "/vol/data/datasets/PanXpresso"))
# Internal metadata: builds the order -> family -> species tree and picks versions (never shown).
METADATA_TSV = Path(os.environ.get("PLANTX_METADATA_TSV", "/vol/data/datasets/PlantXpresso_metadata.tsv"))
# User-facing metadata: shown on the Metadata tab and species pages, offered for download.
WEBSERVER_TSV = Path(os.environ.get("PLANTX_METADATA_WEBSERVER_TSV",
									"/vol/data/datasets/PlantXpresso_metadata_webserver.tsv"))
RELOAD_SECONDS = int(os.environ.get("PLANTX_RELOAD_SECONDS", "600"))

# Column names in the TSV that the browser depends on.
COL_SPECIES = "Species"
COL_ORDER = "Order"
COL_FAMILY = "Family"

# The files inside a version folder that the website offers, by "type".
# Filename = <Species><suffix>. Order matters: this is the order shown on the website.
# Anything not listed here (e.g. the per-version .json) is never shown or served.
FILE_TYPES = {
	"tpms":      (".tpms.tsv",       "Expression (TPM, all transcripts)"),
	"repr_tpms": (".repr.tpms.tsv",  "Expression (TPM, representative transcripts)"),
	"cds":       (".cds.fasta",      "CDS sequences (all)"),
	"repr_cds":  (".repr.cds.fasta", "CDS sequences (representative)"),
	"pep":       (".pep.fasta",      "Peptide sequences (all)"),
	"repr_pep":  (".repr.pep.fasta", "Peptide sequences (representative)"),
}

# "Acer_saccharum_A1E0" -> species "Acer_saccharum", version "A1E0", a=1, e=0
VERSION_RE = re.compile(r"^(?P<species>.+)_(?P<version>A(?P<a>\d+)E(?P<e>\d+))$")


@dataclass
class Species:
	name: str                  # Lotus_japonicus
	version: str               # A1E0
	order: str
	family: str
	meta: dict                 # the full TSV row
	files: dict = field(default_factory=dict)      # type -> Path (only files that exist)
	archives: list = field(default_factory=list)   # older versions, e.g. [Path(.../A1E0.tar.gz)]

	@property
	def display_name(self):
		return self.name.replace("_", " ")

	@property
	def folder(self):
		return DATA_DIR / self.name / self.version

	def files_of_types(self, types):
		"""Paths of the requested file types (all types if `types` is empty)."""
		return [p for t, p in self.files.items() if not types or t in types]


class Catalog:
	def __init__(self):
		self._lock = threading.Lock()
		self._loaded_mtime = None
		self._loaded_at = 0.0
		self.species = {}          # name -> Species
		self.tree = {}             # order -> family -> [Species]
		self.columns = []          # TSV header, for display / filtered TSV downloads
		self.problems = []         # human-readable warnings (shown on /status)

	# -- public ---------------------------------------------------------------------
	def get(self):
		"""Call this at the start of every request: reloads only when needed (cheap)."""
		try:
			mtime = METADATA_TSV.stat().st_mtime
		except FileNotFoundError:
			mtime = None
		stale = (time.time() - self._loaded_at) > RELOAD_SECONDS
		if mtime != self._loaded_mtime or stale:
			with self._lock:                      # only one request reloads at a time
				if mtime != self._loaded_mtime or (time.time() - self._loaded_at) > RELOAD_SECONDS:
					self._load()
					self._loaded_mtime = mtime
					self._loaded_at = time.time()
		return self

	def in_group(self, order=None, family=None):
		"""All species of one order or one family."""
		return [s for s in self.species.values()
				if (order is None or s.order == order) and (family is None or s.family == family)]

	# -- internal -------------------------------------------------------------------
	def _load(self):
		problems = []
		best = {}        # species name -> ((a, e), version, row)   highest version on disk
		columns = []
		seen_in_tsv = set()

		if not METADATA_TSV.exists():
			problems.append(f"Metadata TSV not found: {METADATA_TSV}")
		else:
			with open(METADATA_TSV, newline="", encoding="utf-8") as fh:
				reader = csv.DictReader(fh, delimiter="\t")
				columns = [c.strip() for c in (reader.fieldnames or [])]
				for row in reader:
					row = {(k or "").strip(): (v or "").strip() for k, v in row.items()}
					m = VERSION_RE.match(row.get(COL_SPECIES, ""))
					if not m:
						if row.get(COL_SPECIES):
							problems.append(f"Cannot parse version from Species '{row[COL_SPECIES]}'")
						continue
					name, version = m["species"], m["version"]
					seen_in_tsv.add(name)
					if not (DATA_DIR / name / version).is_dir():
						continue      # tarred (old) or not synced yet -> not the current version
					key = (int(m["a"]), int(m["e"]))
					if name not in best or key > best[name][0]:
						best[name] = (key, version, row)

		for name in sorted(seen_in_tsv - set(best)):
			problems.append(f"'{name}' is in the metadata TSV but has no untarred version folder (not synced yet?)")

		species = {}
		for name, (_, version, row) in best.items():
			sp = Species(
				name=name,
				version=version,
				order=row.get(COL_ORDER) or "Unassigned",
				family=row.get(COL_FAMILY) or "Unassigned",
				meta=row,
			)
			for ftype, (suffix, _label) in FILE_TYPES.items():
				path = sp.folder / f"{name}{suffix}"
				if path.is_file():
					sp.files[ftype] = path
			sp.archives = sorted((DATA_DIR / name).glob("A*E*.tar.gz"))
			species[name] = sp

		# Folders on disk that have no matching (untarred) row in the TSV
		if DATA_DIR.is_dir():
			for d in sorted(DATA_DIR.iterdir()):
				if d.is_dir() and d.name not in species:
					problems.append(f"Folder '{d.name}' has no matching current version in the metadata TSV")

		# order -> family -> [Species], everything sorted alphabetically
		tree = {}
		for sp in sorted(species.values(), key=lambda s: (s.order, s.family, s.name)):
			tree.setdefault(sp.order, {}).setdefault(sp.family, []).append(sp)

		self.species, self.tree, self.columns, self.problems = species, tree, columns, problems


class MetadataTable:
	"""
	The user-facing metadata table (PlantXpresso_metadata_webserver.tsv).
	Same idea as Catalog, much simpler: re-read the file whenever it changes on disk.
	"""

	def __init__(self, path):
		self.path = path
		self._loaded_mtime = "never"
		self.columns = []      # header names, in file order
		self.rows = []         # one dict per line: column name -> value

	def get(self):
		try:
			mtime = self.path.stat().st_mtime
		except FileNotFoundError:
			mtime = None
		if mtime != self._loaded_mtime:
			columns, rows = [], []
			if mtime is not None:
				with open(self.path, newline="", encoding="utf-8") as fh:
					reader = csv.DictReader(fh, delimiter="\t")
					columns = [c.strip() for c in (reader.fieldnames or [])]
					rows = [{(k or "").strip(): (v or "").strip() for k, v in r.items()} for r in reader]
			self.columns, self.rows, self._loaded_mtime = columns, rows, mtime
		return self

	@property
	def key_column(self):
		"""The column naming the species: 'Species' if present, else the first column."""
		if COL_SPECIES in self.columns:
			return COL_SPECIES
		return self.columns[0] if self.columns else None

	def species_name(self, row):
		"""'Lotus_japonicus_A1E0' / 'Lotus_japonicus' / 'Lotus japonicus' -> 'Lotus_japonicus'"""
		value = row.get(self.key_column, "")
		m = VERSION_RE.match(value)
		return (m["species"] if m else value).replace(" ", "_")

	def row_for(self, sp):
		"""This species' row: the exact version if listed, otherwise any row with its name."""
		exact = f"{sp.name}_{sp.version}"
		rows = [r for r in self.rows if self.species_name(r) == sp.name]
		return next((r for r in rows if r.get(self.key_column) == exact), rows[0] if rows else None)


catalog = Catalog()
metadata_table = MetadataTable(WEBSERVER_TSV)
