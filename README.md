# Huntington's Disease RNA-seq Explorer

An interactive R Shiny application for exploring RNA-seq data from post-mortem human prefrontal cortex samples comparing individuals with Huntington's Disease (HD) and neurologically healthy controls.

The application provides tools for exploring sample metadata, normalized gene expression, differential expression results, and expression patterns of individual genes through an intuitive graphical interface.

---

## Features

### Sample Information

* Browse sample metadata in an interactive table
* View summary statistics for disease conditions
* Visualize sample distribution across experimental groups

### Counts Explorer

Explore normalized RNA-seq counts with adjustable filtering options.

Features include:

* Gene variance filtering
* Minimum detectable expression filtering
* Mean vs. variance visualization
* Mean vs. non-zero expression plots
* Principal Component Analysis (PCA)
* Heatmap of the most variable genes

### Differential Expression

Investigate differential gene expression between Huntington's Disease and control samples.

Features include:

* Searchable differential expression table
* Interactive volcano plot
* Significance thresholds using adjusted p-values and log2 fold change

### Gene Expression

Visualize expression of individual genes across experimental conditions.

Features include:

* Search by gene symbol or Ensembl Gene ID
* Violin plots
* Boxplots
* Individual sample points

---

## Technologies Used

* R
* Shiny
* tidyverse
* ggplot2
* DT
* pheatmap
* bslib

---

## Data

The application uses:

* Sample metadata
* Normalized RNA-seq count matrix
* Differential expression results generated with DESeq2

Input files are stored in the `clean_data/` directory.

---

## Project Structure

```text
├── app.R
├── clean_data/
│   ├── sample_info.csv
│   ├── normalized_counts.csv
│   └── deseq_results.csv
└── README.md
```

---

## Running the App

1. Clone this repository.

```bash
git clone https://github.com/yourusername/hd-rnaseq-explorer.git
```

2. Open `app.R` in RStudio.

3. Install the required packages if necessary.

```r
install.packages(c(
  "shiny",
  "tidyverse",
  "DT",
  "ggplot2",
  "pheatmap",
  "bslib"
))
```

4. Run the application.

```r
shiny::runApp()
```

---

## Learning Goals

This project demonstrates:

* Interactive data visualization with Shiny
* RNA-seq exploratory data analysis
* Differential expression result exploration
* PCA and heatmap generation
* Reactive programming in R
* Scientific data communication through interactive dashboards

---

## Future Improvements

Potential future additions include:

* Downloadable figures and tables
* Additional sample metadata filters
* Gene set enrichment analysis (GSEA)
* Interactive heatmaps
* Pathway visualization
* Custom statistical comparisons

---

## License

This project is intended for educational and research purposes.
