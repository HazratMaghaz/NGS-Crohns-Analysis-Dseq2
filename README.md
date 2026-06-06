# RNA-seq Differential Expression Analysis of Crohn's Disease using DESeq2 and GSEA

This repository contains a reproducible RNA-seq differential expression analysis workflow for Crohn's Disease using **R**, **DESeq2**, and **gene set enrichment analysis (GSEA)**.

The aim of this project is to identify differentially expressed genes between Crohn's Disease and control samples, visualize expression patterns, and interpret pathway-level biological signals using Hallmark gene sets.

---

## Project Overview

- **Disease context:** Crohn's Disease
- **Sample type:** Pediatric ileal biopsy RNA-seq count data
- **Total samples:** 304
  - 254 Crohn's Disease samples
  - 50 non-IBD control samples
- **Main analysis:** Differential expression analysis using DESeq2
- **Downstream analysis:** PCA, heatmap, volcano plot, and Hallmark GSEA

---

## Repository Structure

```text
.
├── data/
│   ├── Counts.tsv
│   └── Metadata.tsv
├── scripts/
│   └── analysis_pipeline.R
├── results/
│   ├── figures/
│   └── tables/
├── docs/
│   └── project_notes.md
├── README.md
├── LICENSE
└── .gitignore
```

---

## Analysis Workflow

1. Load RNA-seq count matrix and sample metadata.
2. Prepare the DESeq2 dataset using sample group information.
3. Normalize count data and apply variance-stabilizing transformation.
4. Explore sample-level patterns using PCA and heatmap visualization.
5. Run differential expression analysis using DESeq2.
6. Export DEG result tables.
7. Generate volcano plot for significant gene expression changes.
8. Perform Hallmark gene set enrichment analysis for pathway-level interpretation.

---

## Reported Results

Using an adjusted p-value / FDR threshold of 0.05:

- **Total genes analyzed:** 65,218
- **Significant DEGs:** 11,613
- **Upregulated genes:** 5,352
- **Downregulated genes:** 6,261

---

## Key Outputs

| Output | Description |
|---|---|
| `results/figures/Heatmap.png` | Sample-to-sample distance heatmap |
| `results/figures/PCA_plot.png` | PCA plot for disease/control separation |
| `results/figures/Volcano_plot.png` | Volcano plot of DESeq2 differential expression results |
| `results/figures/Hallmark_plot.png` | Hallmark GSEA visualization |
| `results/tables/` | Generated DEG and enrichment result tables |

---

## Tools and Packages

- R
- DESeq2
- tidyverse
- ggplot2
- pheatmap
- EnhancedVolcano
- clusterProfiler
- msigdbr
- RColorBrewer

---

## How to Run

Clone the repository:

```bash
git clone https://github.com/HazratMaghaz/NGS-Crohns-Analysis-Dseq2.git
cd NGS-Crohns-Analysis-Dseq2
```

Run the analysis script from the repository root:

```r
source("scripts/analysis_pipeline.R")
```

The script uses relative paths, so it should be executed from the main repository directory.

---

## Skills Demonstrated

- RNA-seq count matrix handling
- Metadata-based experimental design
- Differential expression analysis using DESeq2
- PCA and heatmap-based sample exploration
- Volcano plot visualization
- Hallmark pathway enrichment analysis
- Reproducible project organization in R

---

## Author

**Hazrat Maghaz**  
Bioinformatician | Computational Biologist  

- Website: https://hazratmaghaz.tech
- GitHub: https://github.com/HazratMaghaz
- LinkedIn: https://www.linkedin.com/in/hazrat-maghaz-6967b9374/

---

## License

This repository is available under the MIT License. See the `LICENSE` file for details.
