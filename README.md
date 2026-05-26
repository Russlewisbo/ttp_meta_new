# TTP Meta-Analysis — v4 Quarto site

Updated Bayesian meta-analysis of time-to-positivity (TTP) of blood cultures as a predictor of mortality and clinical failure in bacteremia. Built as a Quarto website intended for GitHub Pages.

**Primary pool:** 100 studies (57 baseline + 43 new strict-pool extractions)
**Sensitivity pool:** 140 studies (+ 40 inclusive-pool studies)

## Pages

- `index.qmd` — main mortality meta-analysis (forest plot, subgroups, TTP cutpoint dose-response, sensitivity analyses)
- `clinical_failure.qmd` — persistent bacteremia / clearance / relapse outcomes
- `methods.qmd` — search strategy, PRISMA flow, eligibility, effect-size calculations, prior specification
- `rob.qmd` — QUIPS risk-of-bias traffic-light and summary plots

## Folder layout

```
quarto_site/
├── _quarto.yml                  ← website config + theme
├── index.qmd                    ← main mortality analysis
├── clinical_failure.qmd
├── methods.qmd
├── rob.qmd
├── README.md                    ← you are here
├── .gitignore
├── R/
│   ├── 01_prepare_data.R        ← reads v4 workbook → effect-size table
│   └── 02_fit_models.R          ← 9 brms models, cached
├── data/
│   ├── TTP_MetaAnalysis_Extraction_Complete_v4.xlsx  ← copy from project root
│   └── PRISMA_v4_update.svg                          ← copy from project root
└── .github/
    └── workflows/
        └── quarto-publish.yml   ← GitHub Actions deploy to gh-pages
```

## Local rendering

Prerequisites:
- R ≥ 4.3
- Quarto ≥ 1.5
- C++ toolchain for Stan (Rtools on Windows, Xcode CLT on macOS, build-essential on Linux)

```bash
# 1. Copy the v4 workbook and PRISMA SVG into data/
cp ../TTP_MetaAnalysis_Extraction_Complete_v4.xlsx data/
cp ../PRISMA_v4_update.svg data/

# 2. Install R packages (one-time)
Rscript -e 'install.packages(c("tidyverse","brms","metafor","posterior","bayesplot","ggdist","patchwork","gt","readxl","knitr","scales"))'

# 3. Render the site
quarto render

# Output appears in _site/. Open _site/index.html in a browser.
```

First render compiles all Stan models — **expect 10–20 minutes**. Subsequent renders use the `_brms_cache/` and the `_freeze/` directories, so they finish in seconds.

## Deploying to GitHub Pages

1. Push this folder to a GitHub repo (e.g. `your-username/ttp-meta-v4`).
2. In repo Settings → Pages → Build and deployment source: **GitHub Actions** (or pick the `gh-pages` branch if you prefer the simpler peaceiris-action path).
3. The `.github/workflows/quarto-publish.yml` workflow will:
   - Install R 4.4.1 and Quarto 1.5.57
   - Restore renv library (if you commit a `renv.lock`)
   - Cache Stan compiled models between runs (saves ~15 min/build)
   - Cache Quarto freeze (saves more on re-renders)
   - Render the site and push to `gh-pages`
4. Your site appears at `https://your-username.github.io/ttp-meta-v4/`.

Update the `navbar.right.href` URL in `_quarto.yml` to point at your repo before pushing.

### Cold-start time on GitHub Actions

The first CI run with no cached Stan binaries takes ~30–45 minutes (the matrix of 9 brms models compiles from scratch). After that, the cache keys at the bottom of the workflow keep CI runs under 5 minutes unless `R/02_fit_models.R` or `data/*.xlsx` change.

## Re-running with the inclusive pool as primary

In `index.qmd`, swap the `setup` chunk:

```r
es_all  <- prepare_effect_sizes("data/TTP_MetaAnalysis_Extraction_Complete_v4.xlsx",
                                include_inclusive_pool = TRUE,    # <- changed
                                drop_contaminated_ttp  = TRUE)
```

The cache directory paths in `fit_all_models(..., cache_dir = "_brms_cache/primary")` should be changed to a different folder name to avoid model-collision with the cached strict-pool fits.

## Files referenced from the project root

These need to be copied into `data/` before rendering (the workflow assumes they live alongside the .qmd files for portability):

| Source | Destination |
|---|---|
| `TTP_MetaAnalysis_Extraction_Complete_v4.xlsx` | `data/TTP_MetaAnalysis_Extraction_Complete_v4.xlsx` |
| `PRISMA_v4_update.svg` | `data/PRISMA_v4_update.svg` |

Or symlink them if you prefer:

```bash
ln -s ../TTP_MetaAnalysis_Extraction_Complete_v4.xlsx data/
ln -s ../PRISMA_v4_update.svg data/
```

## Citation

If you use this analysis, please cite both the original v1 report (Feb 2026) and the v4 update. The PROSPERO protocol number is in `PROSPERO_TTP_protocol.md`.
