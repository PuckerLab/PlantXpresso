# PlantXpresso
Scripts for backend development of PlantXpresso webserver hosting the plant species datasets from the PanXpresso data collection.

Tabs: **Download** · **Metadata** · **Expression** (violin plot of one gene) · **PCA** (per species) ·
**Sequence search** (DIAMOND blastx for a CDS / blastp for a protein, best hit's expression).
New tab = one line in `TABS` (app/main.py) + a route + a template.

```
plantxpresso/
├── app/
│   ├── catalog.py        # metadata TSV + folders -> species catalog
│   ├── expression.py     # gene lookup, tissue labels, stats, violin + PCA figures (no web code)
│   ├── search.py         # sequence checks, search jobs, DIAMOND runs (no web code)
│   ├── main.py           # web routes (pages, downloads, plots)
│   └── templates/        # base, index (Download), metadata, expression, pca, species, cli
├── build_cache.py        # .tpms.tsv -> fast cache, PCA, DIAMOND db from .pep.fasta (cron, 2 am)
├── requirements.txt
├── run.sh                # start by hand (development)
└── plantxpresso.service  # systemd: always-on public server on port 3838
```

Data read: `/vol/data/datasets/PlantXpresso/` (species folders), `PlantXpresso_metadata.tsv`,
`PlantXpresso_metadata_webserver.tsv`, `PlantXpresso_logo.png`,
`PlantXpresso_tissue_metadata/<Species>/<version>/Specific_class_labelled_metadata.tsv`.
Cache written: `/vol/data/plantxpresso_cache/<Species>/<version>/` (never inside the rsync target).
Search jobs: `/vol/data/plantxpresso_jobs/<job-id>/` (deleted after 7 days).
Needs the `diamond` program on the VM: `sudo apt install diamond-aligner`.
