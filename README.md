# hiGN expression browser

A small Shiny app that lets anyone in the lab look up whether a gene is expressed
in hiPSC-derived glutamatergic neurons at DIV0 / DIV14 / DIV28 / DIV50, using the
same normalisation and detection rule as `hiGN_gene_check_from_raw.R`.

```
hign-expression-browser/
├── app.R                 the app
├── prepare_data.R        run once, builds data/app_data.rds from raw counts
├── raw_counts_table.txt  input (you supply)
├── data/app_data.rds     built artefact, ~1 MB, commit this
├── manifest.json         built by rsconnect, needed by Connect Cloud
└── .github/workflows/    optional: publish to GitHub Pages
```

---

## Step 1 — build the data file

Put `raw_counts_table.txt` in the project folder and run once:

```r
source("prepare_data.R")
```

This strips Ensembl version suffixes, collapses the 45 duplicate IDs, computes
library sizes on the full matrix, maps every Ensembl ID to a symbol with
`org.Hs.eg.db`, and saves everything to `data/app_data.rds` (well under 1 MB).

All Bioconductor work happens here, once. The app then needs only CRAN packages
— `shiny`, `ggplot2`, `dplyr`, `tidyr`, `DT` — which is what makes it deployable
anywhere, and fast.

Re-run this script whenever the counts change.

## Step 2 — run it locally

```r
shiny::runApp()
```

Type `TLR3` in the search box. You should get a green card, a per-timepoint
table, and a trajectory plot with ACTB/GAPDH and CD19/INS alongside for scale.

## Step 3 — put it on GitHub

```bash
git init
git add app.R prepare_data.R data/app_data.rds README.md .gitignore
git commit -m "hiGN expression browser"
git remote add origin git@github.com:YOURNAME/hign-expression-browser.git
git push -u origin main
```

## Step 4 — deploy to Posit Connect Cloud

Generate the dependency file Connect Cloud reads, then push it:

```r
install.packages("rsconnect")
rsconnect::writeManifest()
```

```bash
git add manifest.json && git commit -m "Add manifest" && git push
```

Then in Connect Cloud: **Publish → Shiny → select this repository → confirm the
branch → set `app.R` as the primary file → Publish.** You get a URL to send to
the lab. Every later `git push` can redeploy automatically.

---

## Optional — GitHub Pages instead

`.github/workflows/deploy-shinylive.yml` compiles the app to WebAssembly with
`shinylive` and publishes it to Pages, so it runs in the visitor's browser with
no server. Enable Pages for the repo with source **GitHub Actions**, then push.

Trade-offs: first load takes 10–30 s while the browser downloads R, and the app
must stay on CRAN packages that have WebAssembly builds (the current five all do
— check <https://repo.r-wasm.org/> before adding any others).

---

## Updating the app

| You changed | Do this |
|---|---|
| The counts | Re-run `prepare_data.R`, commit `data/app_data.rds`, push |
| The detection thresholds | Edit the constants at the top of `prepare_data.R`, re-run, push |
| The app UI or plots | Edit `app.R`, push |
| Added an R package | Re-run `rsconnect::writeManifest()`, push |
