# PlantXpresso
Scripts for backend development of PlantXpresso webserver hosting the plant species datasets from the PanXpresso data collection.

PlantXpresso is a FastAPI web application running on a de.NBI SimpleVM (Web Services, port 3838,
access via Life Science Login).

## Tabs

| Tab | What users can do |
|---|---|
| **Download** | Browse order → family → species; choose file types; download a species, family or order as one `.tar.gz` (browser) or with ready-made `wget` commands (command line, also for the entire collection). |
| **Expression analysis** | **By gene ID:** pick a species, type part of an ID (autocomplete matches anywhere in the ID, ignoring case and punctuation) → violin plot of log2(TPM + 1) per tissue (or across all samples if no tissue labels exist) + min/max/mean/median TPM table; download PNG (600 dpi), SVG or values (TSV). **By sequence:** paste one CDS (DIAMOND blastx) or protein (DIAMOND blastp) in FASTA format → hits with e-value ≤ 1e-5 ranked by bit score; the best hit's expression is plotted. For identifying genes within the same species only, not for ortholog search. |
| **PCA** | Precomputed PCA per species (labelled samples, genes with mean TPM ≥ 1, log2(TPM + 1), z-scored); PC1 vs PC2 coloured by tissue (Glasbey palette); download PNG, SVG or PC coordinates. |
| **Species metadata** | One row per species (assembly, annotation, BUSCO, samples …); search, sort, download as TSV. |

Each species also has a page (`/species/<name>`) with its files, older archived versions and metadata.

## Code layout

```
plantxpresso/
├── app/
│   ├── catalog.py         # metadata TSV + species folders -> catalog (current version per species)
│   ├── expression.py      # gene lookup from the cache, ID matching/autocomplete, stats,
│   │                      #   violin + PCA figures (no web code)
│   ├── search.py          # sequence checks, search jobs, DIAMOND runs, hit ranking (no web code)
│   ├── main.py            # web routes: pages, downloads, plots, search, API
│   └── templates/
│       ├── base.html          # header, logo, tab bar, shared CSS
│       ├── index.html         # Download tab
│       ├── expression.html    # Expression analysis tab (gene ID form + sequence form)
│       ├── _gene_block.html   # stats table + violin + downloads (shared by expression and search results)
│       ├── search_job.html    # sequence search: queued / running / results
│       ├── pca.html           # PCA tab
│       ├── metadata.html      # Species metadata tab
│       ├── species.html       # one species
│       └── cli.html           # command-line download instructions
├── build_cache.py         # nightly: builds the expression cache, PCA and DIAMOND databases
├── plantxpresso.env       # ALL paths/settings, read by the service, cron and run.sh
├── plantxpresso.service   # systemd: the always-on public server (port 3838)
├── run.sh                 # start by hand for development/testing
└── requirements.txt
```

Every feature = a route in `main.py` (+ a template). New tab = one line in `TABS` (main.py) + a route + a template.

## Data and generated files

All paths are set in `plantxpresso.env`.

| Path | Content | Written by |
|---|---|---|
| `/vol/data/datasets/PlantXpresso/<Species>/<A#E#>/` | species data: `.tpms.tsv`, `.cds.fasta`, `.pep.fasta` (+ `.repr` variants); old versions as `A#E#.tar.gz` | rsync from the storage VM (never written by the app) |
| `/vol/data/datasets/PlantXpresso_metadata.tsv` | internal metadata: builds the browse tree, selects the current version | storage VM |
| `/vol/data/datasets/PlantXpresso_metadata_webserver.tsv` | species metadata shown to users | storage VM |
| `/vol/data/datasets/PlantXpresso_tissue_metadata/<Species>/<A#E#>/Specific_class_labelled_metadata.tsv` | tissue label per sample (`run_accession`, `tissue`) | you |
| `/vol/data/datasets/PlantXpresso_logo.png` | logo | you |
| `/vol/data/plantxpresso_cache/<Species>/<A#E#>/` | `expr.npy` + `genes.txt` + `samples.txt` (fast gene lookup), `pca.tsv` + `pca_variance.txt`, `pep.dmnd` (DIAMOND db), `sources.txt` (input fingerprint) | `build_cache.py` |
| `/vol/data/plantxpresso_jobs/<job-id>/` | one folder per sequence search; deleted after 7 days | website / `build_cache.py` |

The website only reads from `/vol/data/datasets`. The cache uses only the full `.tpms.tsv` and `.pep.fasta` files, never the `.repr` ones.

## What runs where

| What | How | When |
|---|---|---|
| Website | systemd service `plantxpresso` (uvicorn, workers, port 3838) | always; starts at boot after `/vol/data` is mounted |
| Cache build | `build_cache.py` from the `ubuntu` crontab | nightly at 02:00 (VM time, usually UTC); only changed species are rebuilt |
| Data sync | rsync + TSV transfer from the storage VM | storage VM's own schedule (not at 02:00) |

Changes to the tissue labels show up immediately in the violin plots; the PCA follows after the next nightly build.

## Setup (once, on the webserver VM)

```bash
sudo apt update && sudo apt install -y python3-venv pigz
# DIAMOND: apt install diamond-aligner, or a binary (here /vol/data/tools/diamond) + PLANTX_DIAMOND in plantxpresso.env

sudo mkdir -p /vol/data/plantxpresso_cache /vol/data/plantxpresso_jobs
sudo chown ubuntu:ubuntu /vol/data/plantxpresso_cache /vol/data/plantxpresso_jobs

cd ~/plantxpresso
python3 -m venv venv
./venv/bin/pip install --no-cache-dir -r requirements.txt

sudo cp plantxpresso.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now plantxpresso

# first full cache build (long; run in tmux)
set -a && . ./plantxpresso.env && set +a && ./venv/bin/python build_cache.py 2>&1 | tee ~/first_build.log
```

Crontab (`crontab -e`, as ubuntu, no sudo):

```
0 2 * * * cd /home/ubuntu/plantxpresso && set -a && . ./plantxpresso.env && set +a && flock -n /tmp/plantx_cache.lock nice ./venv/bin/python build_cache.py >> /home/ubuntu/plantx_cache.log 2>&1
```

## Everyday tasks

| Task | Command(s) |
|---|---|
| Health check | `systemctl status plantxpresso` · `curl -s localhost:3838/status` |
| Logs | `journalctl -u plantxpresso -n 100` (website) · `tail ~/plantx_cache.log` (nightly build) |
| Code or template changed | copy files to the VM, then `sudo systemctl restart plantxpresso` |
| New Python package in `requirements.txt` | `./venv/bin/pip install -r requirements.txt`, then restart |
| A path/setting changed | edit `plantxpresso.env`, then restart (cron picks it up automatically) |
| Service file changed | `sudo cp plantxpresso.service /etc/systemd/system/ && sudo systemctl daemon-reload && sudo systemctl restart plantxpresso` |
| Rebuild one species now | `set -a && . ./plantxpresso.env && set +a && ./venv/bin/python build_cache.py <Species_name>` |
| Test changes privately | `PLANTX_PORT=8000 ./run.sh --reload` + SSH tunnel to port 8000 |

Copy code to the VM without IDE/venv/git folders:

```bash
rsync -av --exclude '.git/' --exclude '.venv/' --exclude 'venv/' --exclude '.idea/' --exclude '__pycache__/' \
  -e "ssh -p <port> -i <key>" ./ ubuntu@<vm-ip>:/home/ubuntu/plantxpresso/
```

Full operations guide (reboots, maintenance, troubleshooting): the "PlantXpresso – server operations guide" doc.

## License

GNU General Public License v3.0, see `LICENSE`.
