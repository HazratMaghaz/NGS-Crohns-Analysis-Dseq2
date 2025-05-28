# NGS-Crohns-Analysis-Dseq2
A comprehensive analysis of differential gene expression in Crohn’s Disease using RNA-seq data, DESeq2, and GSEA in R.
# Differential Gene Expression Analysis of Crohn’s Disease using RNA-seq

This repository contains the full workflow and analysis of a bioinformatics project conducted as part of the *Next Generation Sequencing* course at **NUST – School of Interdisciplinary Engineering and Sciences (SINES)**. The goal was to identify differentially expressed genes (DEGs) in **Crohn’s Disease (CD)** using RNA-seq count data from pediatric ileal biopsy samples and perform downstream analysis, including visualization and gene set enrichment analysis (GSEA).

## 🧬 Project Overview

- **Course Title:** Next Generation Sequencing  
- **Instructor:** Dr. Masood Ur Rehman Kayani  
- **Student:** Hazrat Maghaz  
- **Submitted on:** December 8, 2024

## 📁 Project Structure

📦ngs-crohns-analysis/
├── 📄 NGS_Project_Report_HazratMaghaz.pdf # Final project report
├── 📂 code/ # R scripts used for analysis
├── 📂 figures/ # All plots and visualizations
│ ├── heatmap.png
│ ├── pca_plot.png
│ ├── deg_plot.png
│ ├── volcano_plot.png
│ ├── bubble_plot.png
│ └── gsea_plots/
│ ├── pathway1.png
│ ├── pathway2.png
│ └── ...


---

## 🔬 Analysis Pipeline

### 1. **Dataset Overview**
- **Total samples:** 304  
  - 254 CD (Crohn’s Disease) samples  
  - 50 control (non-IBD) samples  
- **Source:** Pediatric ileal biopsies  
- **Data format:** RNA-seq raw count matrix (genes × samples)

### 2. **Correlation Heatmap**
- Created using Manhattan distance to visualize pairwise similarity.
- Darker blocks = more similar samples; lighter blocks = more different.

### 3. **Principal Component Analysis (PCA)**
- Showed clear separation between CD and control groups.
- PC1 and PC2 together explained 58% of data variance.

### 4. **DESeq2 Differential Expression Analysis**
- **Tool:** DESeq2 package in R.
- **Filtering threshold:** FDR < 0.05
- **Total genes analyzed:** 65,218  
  - Significant DEGs: 11,613  
    - Upregulated: 5,352 (46.08%)  
    - Downregulated: 6,261 (53.91%)

### 5. **Visualization of DEGs**
- **Volcano Plot:** Highlights logFC vs adjusted p-values.
- **Bubble Plot:** Top significant gene with p-value = 3.53e-71 and padj = 1.12e-66.

### 6. **Gene Set Enrichment Analysis (GSEA)**
- **Gene sets used:** HALLMARK from MSigDB
- **Result:**  
  - 19 pathways significantly positively enriched  
  - 3 positively enriched but not significant  
  - 1 pathway negatively enriched (e.g., oxidative phosphorylation)

---

## 📊 Key Figures

- Heatmap of all 304 samples.
- PCA plot showing case/control clustering.
- DEG volcano plot and top gene bubble plot.
- GSEA enrichment plots (top positively and negatively enriched pathways).

---

## 💻 Technologies Used

- **R**  
- **RStudio**  
- **DESeq2**  
- **ggplot2**  
- **clusterProfiler (for GSEA)**

---

## 🧠 Learning Outcomes

- Hands-on experience with RNA-seq analysis pipeline.
- Mastery of DESeq2 for identifying DEGs.
- GSEA interpretation using real disease-related data.
- Bioinformatics visualization and critical result interpretation.

---

## 📜 License

This project is academic coursework and is not licensed for commercial use.

---

## 📩 Contact

For questions or collaboration, contact: **Hazrat Maghaz** – maaz28608@gmail.com.


