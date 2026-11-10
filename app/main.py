"""
main.py  --  the PlantXpresso web application (FastAPI).

Every feature = one URL + one function below (+ optionally an HTML template).

Download levels: "species", "family", "order" (and "all" = whole collection, command line only).
For every level the user gets two options:
  * Browser:      /download/{level}/{name}   -> one compressed .tar.gz, built on the fly
  * Command line: /cli/{level}/{name}        -> a page with ready-to-copy wget commands that
                                                fetch the individual files listed by
                                                /manifest/{level}/{name} (resumable, no
                                                compression work for the server)
Both accept ?types=tpms&types=cds ... to choose file types (none given = all types).

Tabs (the bar under the header; add a tab = one line in TABS + its route + template):
  /                                       Download tab: browse order -> family -> species
  /metadata                               Metadata tab: the webserver metadata table
  /expression                             Expression tab: violin plot for one gene
  /pca                                    PCA tab: precomputed PCA of one species
  /search                                 Sequence search tab: DIAMOND blastp / blastx

Other routes:
  /species/{name}                         one species: files, metadata, older versions
  /files/{name}/{version}/{filename}      one single file
  /archive/{name}/{archive}               an older, archived version (e.g. A1E0.tar.gz)
  /metadata/download                      the webserver metadata table as a TSV file
  /expression/plot?species=&gene=&fmt=     the violin plot as png / svg, or the values as tsv
  /api/genes?species=&q=                  autocomplete: matching gene IDs (JSON list)
  /pca/plot?species=&fmt=                  the PCA plot as png / svg, or the coordinates as tsv
  /search/{job_id}                        one search: queued / running / results
  /logo.png                               the site logo
  /api/species                            the catalog as JSON
  /status                                 health check + warnings
"""

import os
import shutil
import subprocess
import threading
from pathlib import Path
from typing import List
from urllib.parse import quote, urlencode, urlparse

from fastapi import FastAPI, Form, HTTPException, Query, Request
from fastapi.responses import (FileResponse, HTMLResponse, PlainTextResponse, RedirectResponse,
                               Response, StreamingResponse)
from fastapi.templating import Jinja2Templates

from . import expression as ex
from . import search
from .catalog import DATA_DIR, FILE_TYPES, WEBSERVER_TSV, catalog, metadata_table

# ---- Settings (override with environment variables, see INSTRUCTIONS) -------------
LOGO_PATH = Path(os.environ.get("PLANTX_LOGO", "/vol/data/datasets/PlantXpresso_logo.png"))
PUBLIC_URL = os.environ.get("PLANTX_PUBLIC_URL", "").rstrip("/")   # only needed if CLI links look wrong
MAX_ARCHIVES = int(os.environ.get("PLANTX_MAX_ARCHIVES", "3"))      # parallel .tar.gz builds per worker

LEVELS = ("species", "family", "order", "all")

# The tab bar: (key, label, URL). To add a tab later: add a line here, a route, a template.
TABS = [
    ("download", "Download", "/"),
    ("expression", "Expression analysis", "/expression"),
    ("pca", "PCA", "/pca"),
    ("metadata", "Species metadata", "/metadata"),
]

app = FastAPI(title="PlantXpresso")
templates = Jinja2Templates(directory=Path(__file__).parent / "templates")


# =============================================================== helpers ============
def render(request, template, active, **context):
    """
    Fill an HTML template. `active` = which tab is highlighted ("download", "metadata", ...).
    `base` = URL prefix, in case the site runs under a sub-path.
    """
    context.update(request=request, base=request.scope.get("root_path", ""), file_types=FILE_TYPES,
                   tabs=TABS, active=active)
    return templates.TemplateResponse(request, template, context)


def public_base(request):
    """The site's public address, e.g. https://plantxpresso.example.org  (used in wget commands)."""
    return PUBLIC_URL or str(request.base_url).rstrip("/")


def members_for(level, name):
    """The Species objects belonging to one download level, e.g. ('family', 'Fabaceae')."""
    cat = catalog.get()
    if level == "species":
        members = [cat.species[name]] if name in cat.species else []
    elif level == "family":
        members = cat.in_group(family=name)
    elif level == "order":
        members = cat.in_group(order=name)
    elif level == "all":
        members = list(cat.species.values())
    else:
        raise HTTPException(404, f"Unknown level '{level}'. Use one of {LEVELS}")
    if not members:
        raise HTTPException(404, f"No {level} called '{name}'")
    return sorted(members, key=lambda s: s.name)


def check_types(types):
    bad = [t for t in types if t not in FILE_TYPES]
    if bad:
        raise HTTPException(400, f"Unknown file type(s): {bad}. Allowed: {list(FILE_TYPES)}")
    return types


def archive_name(level, name, members, types):
    """e.g. Lotus_japonicus_A1E0  /  PlantXpresso_Fabaceae_tpms_cds"""
    stem = f"{name}_{members[0].version}" if level == "species" else f"PlantXpresso_{name}"
    if types and set(types) != set(FILE_TYPES):
        stem += "_" + "_".join(types)
    return stem


def file_url_path(sp, path):
    """URL path of one file: /files/<species>/<version>/<filename>"""
    return f"/files/{sp.name}/{sp.version}/{path.name}"


# ---- streaming .tar.gz (browser option) --------------------------------------------
_active = 0                      # .tar.gz downloads currently running in this worker
_active_lock = threading.Lock()  # so two threads don't change _active at the same moment


def stream_tar_gz(files, name):
    """
    Stream a .tar.gz to the browser while it is being created: no temporary file on disk,
    and the download starts immediately. At most MAX_ARCHIVES run at once per worker;
    beyond that users get a friendly "try again" message instead of everyone slowing down.
    """
    if not files:
        raise HTTPException(404, "No files match this selection.")
    if _active >= MAX_ARCHIVES:
        raise HTTPException(503, "Many downloads are running right now - please try again in a few "
                                 "minutes, or use the command-line download.")
    # paths inside the archive look like  Lotus_japonicus/A1E0/Lotus_japonicus.tpms.tsv
    rel_paths = [str(p.relative_to(DATA_DIR)) for p in files]
    compress = ["--use-compress-program=pigz -p 2"] if shutil.which("pigz") else ["-z"]
    cmd = ["tar", "-c", *compress, "-f", "-", "-C", str(DATA_DIR), "--", *rel_paths]

    def generate():
        global _active
        with _active_lock:
            _active += 1
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        try:
            while chunk := proc.stdout.read(1024 * 1024):     # 1 MB at a time
                yield chunk
        finally:                                              # also runs if the user cancels
            if proc.poll() is None:
                proc.kill()
            proc.wait()
            with _active_lock:
                _active -= 1

    return StreamingResponse(
        generate(),
        media_type="application/gzip",
        headers={"Content-Disposition": f'attachment; filename="{name}.tar.gz"'},
    )


# ================================================================= pages ============
@app.get("/", response_class=HTMLResponse)
def index(request: Request):
    cat = catalog.get()
    return render(request, "index.html", "download", tree=cat.tree, n_species=len(cat.species))


@app.get("/species/{name}", response_class=HTMLResponse)
def species_page(request: Request, name: str):
    sp = catalog.get().species.get(name)
    if sp is None:
        raise HTTPException(404, f"Species '{name}' not found")
    table = metadata_table.get()
    return render(request, "species.html", "download", sp=sp,
                  meta_columns=table.columns, meta_row=table.row_for(sp),
                  file_url=lambda p: file_url_path(sp, p))


@app.get("/metadata", response_class=HTMLResponse)
def metadata_page(request: Request):
    """Metadata tab: the whole webserver metadata table, searchable and sortable."""
    table = metadata_table.get()
    known = catalog.get().species
    # for each row: the species name to link to, or None if it has no downloadable data
    rows = [(r, table.species_name(r) if table.species_name(r) in known else None) for r in table.rows]
    return render(request, "metadata.html", "metadata",
                  columns=table.columns, key_column=table.key_column, rows=rows)


@app.get("/cli/{level}/{name}", response_class=HTMLResponse)
def cli_page(request: Request, level: str, name: str, types: List[str] = Query(default=[])):
    """Page with copy-paste wget commands for one level."""
    members = members_for(level, name)
    check_types(types)
    files = [p for sp in members for p in sp.files_of_types(types)]
    site = public_base(request)
    query = ("?" + urlencode([("types", t) for t in types])) if types else ""
    stem = "PlantXpresso_all" if level == "all" else archive_name(level, name, members, types)
    return render(
        request, "cli.html", "download",
        level=level, name=name, stem=stem, n_species=len(members), n_files=len(files),
        size_gb=sum(p.stat().st_size for p in files) / 1e9,
        types=[FILE_TYPES[t][1] for t in types] or ["all file types"],
        manifest_url=f"{site}/manifest/{level}/{quote(name)}{query}",
        metadata_url=f"{site}/metadata/download",
        # wget --cut-dirs: drop the leading "files/" (and any sub-path) from saved file paths
        cut_dirs=1 + len([p for p in urlparse(site).path.split("/") if p]),
    )


# ===================================================== expression + PCA tabs ========
def species_or_404(name):
    sp = catalog.get().species.get(name)
    if sp is None:
        raise HTTPException(404, f"Species '{name}' not found")
    return sp


def all_species():
    """Species for the drop-down menus, sorted by name."""
    return sorted(catalog.get().species.values(), key=lambda s: s.name)


def send(content, media_type, filename, download):
    """Return bytes to the browser; download=True makes it a "save as" file."""
    headers = {"Content-Disposition": f'attachment; filename="{filename}"'} if download else {}
    return Response(content, media_type=media_type, headers=headers)


MEDIA = {"png": "image/png", "svg": "image/svg+xml", "tsv": "text/tab-separated-values"}


def gene_block(sp, gene):
    """
    Everything the shared gene block (_gene_block.html) needs: the expression result
    and its statistics table. Used by the Expression tab and by search results.
    """
    result = ex.gene_expression(sp, gene)            # raises LookupError if not found
    stats = [("All samples", ex.summary(result.values))]
    if result.tissue_wise:
        stats += [(label, ex.summary(values)) for label, values in result.groups]
    return {"result": result, "stats": stats}


@app.get("/expression", response_class=HTMLResponse)
def expression_page(request: Request, species: str = "", gene: str = "", scale: bool = False):
    """Expression tab. Without a gene: just the form. With a gene: stats + plot + downloads."""
    context = dict(species_list=all_species(), species=species, gene=gene.strip(), scale=scale,
                   result=None, stats=[], error=None, suggestions=[], **search_form())
    if species and gene.strip():
        sp = species_or_404(species)
        try:
            context.update(gene_block(sp, gene.strip()))
            context["gene"] = context["result"].gene     # the real ID, if the typed one was resolved
        except LookupError as e:
            context["error"] = str(e)
            context["suggestions"] = ex.suggest_genes(sp, gene.strip())
    return render(request, "expression.html", "expression", **context)


@app.get("/api/genes")
def api_genes(species: str, q: str = ""):
    """
    Autocomplete for the Expression tab: up to 20 IDs of this species matching what was typed
    (anywhere in the ID, ignoring case and punctuation). The browser calls this while the user types.
    """
    sp = catalog.get().species.get(species)
    return ex.suggest_genes(sp, q) if sp else []


@app.get("/expression/plot")
def expression_plot(species: str, gene: str, fmt: str = "png", dpi: int = 110,
                    scale: bool = False, download: bool = False):
    sp = species_or_404(species)
    if fmt not in MEDIA:
        raise HTTPException(400, "fmt must be png, svg or tsv")
    try:
        result = ex.gene_expression(sp, gene)
    except LookupError as e:
        raise HTTPException(404, str(e))
    stem = f"{sp.name}_{sp.version}_{result.gene}_expression"
    if fmt == "tsv":   # the numbers behind the plot: one line per sample
        labels = ex.read_tissue_labels(ex.tissue_metadata_path(sp))
        lines = ["run_accession\ttissue\tTPM"] + [
            f"{s}\t{labels.get(s, '')}\t{v:g}" for s, v in zip(result.samples, result.values)]
        content = "\n".join(lines) + "\n"
    else:
        content = ex.figure_bytes(ex.violin_figure(result, scale), fmt, min(dpi, 600))
    return send(content, MEDIA[fmt], f"{stem}.{fmt}", download)


@app.get("/pca", response_class=HTMLResponse)
def pca_page(request: Request, species: str = ""):
    context = dict(species_list=all_species(), species=species, pca=None, error=None)
    if species:
        sp = species_or_404(species)
        pca = ex.load_pca(sp)
        if pca is None:
            context["error"] = (f"No PCA available for {sp.display_name} yet "
                                "(still being prepared, or too few samples).")
        else:
            rows, variance = pca
            context["pca"] = dict(n_samples=len(rows),
                                  n_tissues=len({r["tissue"] for r in rows}),
                                  variance=[v * 100 for v in variance])
    return render(request, "pca.html", "pca", **context)


@app.get("/pca/plot")
def pca_plot(species: str, fmt: str = "png", dpi: int = 110, download: bool = False):
    sp = species_or_404(species)
    if fmt not in MEDIA:
        raise HTTPException(400, "fmt must be png, svg or tsv")
    pca = ex.load_pca(sp)
    if pca is None:
        raise HTTPException(404, "No PCA available for this species")
    stem = f"{sp.name}_{sp.version}_PCA"
    if fmt == "tsv":   # the precomputed coordinates, exactly as stored
        return send((ex.cache_dir(sp) / "pca.tsv").read_bytes(), MEDIA["tsv"], f"{stem}.tsv", download)
    rows, variance = pca
    content = ex.figure_bytes(ex.pca_figure(sp, rows, variance), fmt, min(dpi, 600))
    return send(content, MEDIA[fmt], f"{stem}.{fmt}", download)


# =========================================================== sequence search ========
SEARCH_NOTE = ("This search finds the most similar gene in the selected species. Use it for a "
               "sequence from that same species (e.g. a gene known under another ID). It is "
               "based on local alignment and is not suitable for finding orthologs in other species.")


def search_species():
    """Species that have a DIAMOND database (others can't be searched yet)."""
    return [sp for sp in all_species() if search.database(sp).is_file()]


def search_form(species="", seq_type="cds", sequence="", errors=None):
    """What the sequence-search form on the Expression analysis tab needs."""
    return dict(search_species_list=search_species(), search_species=species,
                seq_type=seq_type, sequence=sequence, search_errors=errors or [], note=SEARCH_NOTE)


@app.get("/search")
def search_page(request: Request):
    """Old address of the Sequence search tab, which is now part of Expression analysis."""
    return RedirectResponse(f"{request.scope.get('root_path', '')}/expression#sequence", status_code=303)


@app.post("/search", response_class=HTMLResponse)
def search_submit(request: Request, species: str = Form(""), seq_type: str = Form(""),
                  sequence: str = Form("", max_length=200_000)):
    """Check the input; if valid, create a job, start it, and send the user to its page."""
    sp = catalog.get().species.get(species)
    header, seq, errors, warnings = search.check_query(sequence, seq_type)
    if sp is None or not search.database(sp).is_file():
        errors.insert(0, "Choose a species from the list.")
    if errors:  # show the Expression analysis tab again, with the user's input and what's wrong
        return render(request, "expression.html", "expression",
                      species_list=all_species(), species="", gene="", scale=False,
                      result=None, stats=[], error=None, suggestions=[],
                      **search_form(species, seq_type, sequence, errors))
    job_id = search.create_job(sp, header, seq, seq_type, warnings)
    search.start_job(sp, job_id)
    # 303 = "now go to this page with GET"; refreshing it never resubmits the search
    return RedirectResponse(f"{request.scope.get('root_path', '')}/search/{job_id}", status_code=303)


@app.get("/search/{job_id}", response_class=HTMLResponse)
def search_job_page(request: Request, job_id: str, scale: bool = False):
    job = search.read_job(job_id)
    if job is None:
        raise HTTPException(404, "Search not found (results are kept for 7 days).")
    sp = catalog.get().species.get(job["species"])
    context = dict(job=job, sp=sp, note=SEARCH_NOTE, block=None, block_error=None, scale=scale)
    if job["status"] == "done" and job["hits"] and sp is not None:
        top = job["hits"][0]["id"]               # highest bit score
        try:
            block = gene_block(sp, top)
            block.update(species=sp.name, gene=top)
            context["block"] = block
        except LookupError as e:
            context["block_error"] = f"Top hit {top}: {e}"
    return render(request, "search_job.html", "expression", **context)

# ============================================================= downloads ============
@app.get("/download/{level}/{name}")
def download_archive_of_level(level: str, name: str, types: List[str] = Query(default=[])):
    """Browser option: one .tar.gz for a species, family or order."""
    if level == "all":
        raise HTTPException(400, "The whole collection is too large for one browser download - "
                                 "please use the command-line download.")
    members = members_for(level, name)
    check_types(types)
    files = [p for sp in members for p in sp.files_of_types(types)]
    return stream_tar_gz(files, archive_name(level, name, members, types))


@app.get("/manifest/{level}/{name}", response_class=PlainTextResponse)
def manifest(request: Request, level: str, name: str, types: List[str] = Query(default=[])):
    """Command-line option: plain list of file URLs, one per line (input for wget -i)."""
    members = members_for(level, name)
    check_types(types)
    site = public_base(request)
    return "\n".join(site + file_url_path(sp, p) for sp in members for p in sp.files_of_types(types)) + "\n"


@app.get("/files/{name}/{version}/{filename}")
def download_file(name: str, version: str, filename: str):
    sp = catalog.get().species.get(name)
    # Only serve files the catalog knows about -> users can never reach other paths.
    path = None
    if sp and sp.version == version:
        path = next((p for p in sp.files.values() if p.name == filename), None)
    if path is None:
        raise HTTPException(404, "File not found (it may have been replaced by a newer version)")
    return FileResponse(path, filename=filename)


@app.get("/archive/{name}/{archive}")
def download_old_version(name: str, archive: str):
    sp = catalog.get().species.get(name)
    path = next((p for p in sp.archives if p.name == archive), None) if sp else None
    if path is None:
        raise HTTPException(404, "Archived version not found")
    return FileResponse(path, filename=f"{name}_{archive}")


@app.get("/metadata/download")
def download_metadata():
    """The webserver metadata table, exactly as the file on disk."""
    if not WEBSERVER_TSV.is_file():
        raise HTTPException(404, "Metadata table not found")
    return FileResponse(WEBSERVER_TSV, media_type="text/tab-separated-values",
                        filename="PlantXpresso_species_metadata.tsv")


@app.get("/logo.png")
def logo():
    if not LOGO_PATH.is_file():
        raise HTTPException(404, "Logo not found")
    return FileResponse(LOGO_PATH, headers={"Cache-Control": "public, max-age=86400"})


# ============================================================ API/status ============
@app.get("/api/species")
def api_species():
    return [
        {"species": sp.name, "version": sp.version, "order": sp.order, "family": sp.family,
         "files": {t: p.name for t, p in sp.files.items()},
         "archived_versions": [a.name for a in sp.archives]}
        for sp in sorted(catalog.get().species.values(), key=lambda s: s.name)
    ]


@app.get("/status")
def status():
    cat = catalog.get()
    return {"species": len(cat.species), "orders": len(cat.tree),
            "families": sum(len(f) for f in cat.tree.values()),
            "metadata_tsv_found": catalog.columns != [],
            "metadata_webserver_tsv_found": WEBSERVER_TSV.is_file(),
            "logo_found": LOGO_PATH.is_file(),
            "expression_cache_missing": [s for s in sorted(cat.species)
                                         if not (ex.cache_dir(cat.species[s]) / "expr.npy").is_file()],
            "warnings": cat.problems}
