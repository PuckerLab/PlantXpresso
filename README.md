
# PlantXpresso


<p align="center">
  <img width="200" height="200" alt="PlantXpresso_Botanical_RNA_Coffee_Cup-removebg-preview" src="https://github.com/user-attachments/assets/e45d4c35-2973-49ec-a867-0a41ee51d3c6" />
</p>

Scripts for backend development of PlantXpresso webserver hosting the plant species datasets from the PanXpresso data collection.

PlantXpresso is a FastAPI web application running on a de.NBI SimpleVM

## PlantXpresso access link

https://simplevm.bi.denbi.de/panxpressowebserver_100

Access requires Life Science AAI Login
 

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

#### AI Statement
The core expression plot, and PCA code blocks were manually scripted. Claude Sonnet 5.5
was used to debug and optimize the scripts. The tissue metadata annotations were done manually
and PlantScience.ai was used to obtain information about some specific plant tissue types and
anatomical structures. Claude Opus 5.5 was used to construct the web browser backend scripts 
involving FastAPI, Uvicorn, and Jinja2. ChatGPT's GPT6 was used to design the logo of the 
web browser based on a manually described concept prompt.
