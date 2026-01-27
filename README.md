# Phenotypic Plasticity: Expression Variability × Methylation Conservation

This repository contains analysis code for the **phenotypic plasticity project**, which links **gene expression variability across epithelial cell types (single‑cell RNA‑seq)** to **DNA methylation divergence (WGBS PWD)** across regulatory contexts.

The core scientific question is:

> Do genes that are more *expression‑plastic* across epithelial cell types show **more conserved** or **more divergent** DNA methylation in promoters and enhancers?

The repository is organized around a **single main analysis notebook**, supported by several **mapping, preprocessing, and quality‑control Rmd files**. This README explains the role of each file and how they fit together.

---

## 🧠 Conceptual Overview

### Expression (Plasticity Proxy)

* Data sources: **SCP259** (normal vs inflamed colon) and **SCP1162** (normal vs tumor; MMRd/MMRp)
* Epithelial cells only
* Pipeline:

  1. Collapse transcript‑level rows → **gene‑level counts**
  2. Per‑cell normalization: **CP10k**
  3. Log transform: **log2(CP10k + 1)**
  4. Average within epithelial **cell types**
  5. For each gene:

     * `mean_expr`: mean across epithelial cell types
     * `var_expr`: variance across epithelial cell types (plasticity proxy)

### Methylation (PWD)

* Data source: Whole‑genome bisulfite sequencing (WGBS)

* Metric: **Pairwise methylation divergence (PWD)**

* Interval‑level PWD is precomputed elsewhere and imported here

* Regulatory contexts:

  * **pELS** (proximal enhancers)
  * **dELS** (distal enhancers)
  * **TSS500** (±500 bp around 5′‑most TSS per gene)

* Gene‑level PWD is computed as:

```
wmean_pwd = weighted.mean(pwd_mean, w = n_cpg)
```

* Enhancers are assigned to genes via **nearest 5′‑most TSS**

---

## ⭐ Main Analysis File (Start Here)

### **`Expression_Phenotypic.Rmd`**

🚨 **This is the main analysis notebook for the project.**

This file integrates **gene‑level expression plasticity** and **gene‑level methylation divergence** and produces all primary results and figures.

**What this file does:**

1. **Processes single‑cell RNA‑seq data**

   * SCP259 and SCP1162 (epithelial only)
   * Gene‑level expression matrices per epithelial cell type
   * Computes `mean_expr` and `var_expr`

2. **Imports gene‑level PWD tables**

   * From pELS, dELS, and TSS500
   * Uses weighted mean PWD per gene per sample pair

3. **Merges expression and methylation at the gene level**

   * Ensures consistent gene identifiers
   * Keeps joins explicit and reproducible

4. **Generates figures**

   * `mean_expr` vs `wmean_pwd`
   * `var_expr` vs `wmean_pwd`
   * Stratified by:

     * Regulatory feature (pELS / dELS / TSS500)
     * Sample‑pair category (e.g. Normal pairs)

5. **Creates summary tables**

   * Top genes by expression
   * Top genes by expression variance
   * Genes with highest / lowest methylation divergence

👉 **If you want the biological results, this is the file to read and run.**

---

## 🧬 Gene & Interval Mapping Files

These files define how transcripts, genes, TSSs, and intervals are mapped. They are critical for correctness but are not meant to be run repeatedly once validated.

### `TSS_gene_mapping.Rmd`

Defines the **gene‑level TSS representation** used throughout the project.

* Uses RefSeq annotations (`refseq_all`)
* Selects the **5′‑most TSS per gene** (strand‑aware)
* Builds **TSS500 windows**: `[TSS − 500, TSS + 500]`
* Output is used for:

  * Promoter‑level PWD (TSS500)
  * Enhancer → gene assignment

This file encodes a **key methodological assumption** of the project.

---

### `Refseq_PWD.Rmd`

Handles **gene‑level methylation divergence calculation** from interval‑level PWD.

* Inputs:

  * Interval‑level PWD tables (from `run_all_pairs_one_interval_dt`)
  * Interval annotations (pELS / dELS / TSS500)
* Steps:

  * Assign intervals to genes (directly for TSS500; nearest TSS for enhancers)
  * Collapse interval‑level PWD → **gene‑level weighted mean PWD**

This file produces the **gene‑level PWD tables** consumed by `Expression_Phenotypic.Rmd`.

---

## 🔁 Transcript / Gene ID Mapping (Support Files)

These files were created to **resolve transcript‑to‑gene mapping issues** and to validate consistency across annotation systems. They are mainly for transparency and reproducibility.

### `ENST_mapping_Ensembl.Rmd`

* Maps Ensembl transcript IDs (ENST) → gene identifiers
* Used to understand transcript duplication and gene collapsing
* Informs how expression matrices are summed at the gene level

---

### `ENST_mapping_UCSC.Rmd`

* Alternative transcript‑to‑gene mapping using UCSC / RefSeq‑style annotations
* Cross‑checks consistency with Ensembl mappings

---

### `ENST_data_check.Rmd`

**Quality‑control notebook**

* Diagnoses:

  * Duplicate transcripts
  * One‑to‑many transcript–gene mappings
  * Potential gene inflation or loss

This file does **not** generate final analysis outputs, but documents decisions made upstream.

---

## 🧪 Which Files Matter for What?

| Goal                        | File(s)                     |
| --------------------------- | --------------------------- |
| Main biological results     | `Expression_Phenotypic.Rmd` |
| Gene‑level PWD construction | `Refseq_PWD.Rmd`            |
| TSS definition              | `TSS_gene_mapping.Rmd`      |
| Transcript ↔ gene mapping   | `ENST_mapping_*.Rmd`        |
| Data QC / validation        | `ENST_data_check.Rmd`       |

---

## 📌 Reproducibility Notes

* All joins are performed at the **gene level** using consistent identifiers
* Enhancer–gene assignment uses **nearest 5′‑most TSS**, not overlap
* Gene‑level PWD uses **CpG‑weighted means**
* Expression plasticity is defined as **variance across epithelial cell types**, not across cells

---

## 🧬 Final Note

If you are new to this repository:

➡️ **Start with `Expression_Phenotypic.Rmd`**

All other files exist to support, justify, or validate the assumptions used there.
