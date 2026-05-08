library(shiny)
library(tidyverse)
library(DT)
library(ggplot2)
library(pheatmap)
library(bslib)

# Load cleaned data
sample_info <- read.csv("clean_data/sample_info.csv")


counts <- read.csv(
  "clean_data/normalized_counts.csv",
  check.names = FALSE,
  row.names = 1
)
counts <- tibble::rownames_to_column(counts, var = "gene_id")
counts[, -1] <- lapply(counts[, -1], as.numeric)

de <- read.csv("clean_data/deseq_results.csv")
gene_lookup <- de %>%
  select(gene_id, symbol) %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct(gene_id, symbol)

gene_choices <- setNames(
  gene_lookup$gene_id,
  paste0(gene_lookup$symbol, " | ", gene_lookup$gene_id)
)

condition_colors <- c(
  "Control" = "#7B8FA1",
  "HD" = "#8E5EA2"
)

ui <- fluidPage(
  theme = bs_theme(
    bootswatch = "flatly",
    primary = "#2C3E50"
  ),
  
  titlePanel("Huntington's Disease RNA-seq Explorer"),
  
  p("This Shiny app explores RNA-seq data from post-mortem human prefrontal cortex samples comparing Huntington's Disease cases with neurologically healthy controls."),
  
  p("Use the tabs below to explore sample metadata, normalized counts, differential expression results, and individual gene expression patterns."),
  
  tabsetPanel(
    tabPanel("Sample Info",
             h3("Sample Metadata"),
             DTOutput("sample_table"),
             h3("Condition Summary"),
             tableOutput("condition_summary"),
             plotOutput("condition_plot", height = "500px")
    ),
    
    tabPanel("Counts Explorer",
             h3("Counts Matrix Explorer"),
             p("Explore normalized gene counts and filter genes based on variability and detection across samples."),
             p(em("Higher variance genes are more informative for PCA and heatmaps. The non-zero filter removes genes with little or no expression across samples.")),
             
             sidebarLayout(
               sidebarPanel(
                 sliderInput("variance_percentile", "Minimum gene variance percentile:",
                             min = 0, max = 100, value = 50),
                 sliderInput("nonzero_min", "Minimum samples with detectable expression:",
                             min = 0, max = ncol(counts) - 1, value = 5),
                 h4("Filtering Summary"),
                 verbatimTextOutput("counts_summary")
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Mean vs Variance",
                            plotOutput("mean_variance_plot", height = "500px")),
                   tabPanel("Mean vs Non-Zero",
                            plotOutput("mean_zero_plot", height = "500px")),
                   tabPanel("PCA",
                            plotOutput("pca_plot", height = "500px")),
                   tabPanel("Heatmap",
                            plotOutput("heatmap_plot", height = "500px"))
                 )
               )
             )
    ),
    
    tabPanel("Differential Expression",
             h3("Differential Expression Results"),
             p("Explore genes that differ between Huntington's Disease and Control samples."),
             p(em("log2 fold change shows the direction and size of expression change. Adjusted p-value controls for multiple testing.")),
             
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_search", "Search gene:", "")
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Table", DTOutput("de_table")),
                   tabPanel("Volcano Plot",
                            p("Each point represents one gene. The vertical dashed lines mark log2 fold change cutoffs of -1 and 1, and the horizontal dashed line marks adjusted p-value = 0.05."),
                            plotOutput("volcano_plot", height = "500px"))
                 )
               )
             )
    ),
    
    tabPanel("Gene Expression",
             h3("Individual Gene Expression"),
             p("Search for a gene symbol or Ensembl gene ID to compare normalized expression between HD and Control samples."),
             p(em("Example genes: HTT, BDNF, GFAP, MBP, SNAP25")),
             
             selectizeInput("selected_gene", "Choose gene:",
                            choices = NULL,
                            selected = NULL,
                            options = list(
                              placeholder = "Type a gene symbol or ENSG ID",
                              maxOptions = 20
                            )
             ),
             
             p(strong("Grouped by:"), "Condition (HD vs Control)"),
             
             plotOutput("gene_boxplot", height = "500px")
    )
  )
)

server <- function(input, output, session) {
  updateSelectizeInput(
    session,
    "selected_gene",
    choices = gene_choices,
    selected = character(0),
    server = TRUE
  )
  
  # Sample info tab
  output$sample_table <- renderDT({
    datatable(sample_info)
  })
  
  output$condition_summary <- renderTable({
    table(sample_info$condition)
  })
  
  output$condition_plot <- renderPlot({
    ggplot(sample_info, aes(x = condition, fill = condition)) +
      geom_bar() +
      scale_fill_manual(values = condition_colors) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "none"
      ) +
      labs(
        title = "Sample Distribution by Condition",
        x = "Condition",
        y = "Number of Samples"
      )
  })
  
  # Counts explorer tab
  gene_stats <- reactive({
    count_values <- counts[, -1]
    
    data.frame(
      gene_id = counts$gene_id,
      mean_count = rowMeans(count_values),
      variance = apply(count_values, 1, var),
      nonzero_samples = rowSums(count_values > 0)
    )
  })
  
  filtered_gene_stats <- reactive({
    stats <- gene_stats()
    var_cutoff <- quantile(stats$variance, input$variance_percentile / 100)
    
    stats %>%
      mutate(pass_filter = variance >= var_cutoff &
               nonzero_samples >= input$nonzero_min)
  })
  
  output$counts_summary <- renderPrint({
    stats <- filtered_gene_stats()
    
    total_genes <- nrow(stats)
    passing <- sum(stats$pass_filter)
    
    cat("Number of samples:", ncol(counts) - 1, "\n")
    cat("Total genes:", total_genes, "\n")
    cat("Genes passing filter:", passing, "\n")
    cat("Genes filtered out:", total_genes - passing, "\n")
    cat("Percent passing:", round(100 * passing / total_genes, 2), "%\n")
  })
  
  output$mean_variance_plot <- renderPlot({
    ggplot(filtered_gene_stats(),
           aes(x = mean_count, y = variance, color = pass_filter)) +
      geom_point(alpha = 0.5) +
      scale_x_log10() +
      scale_y_log10() +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold")
      ) +
      labs(
        title = "Mean Count vs Variance",
        x = "Mean normalized count",
        y = "Variance",
        color = "Included after filter"
      )
  })
  
  output$mean_zero_plot <- renderPlot({
    ggplot(filtered_gene_stats(),
           aes(x = mean_count,
               y = nonzero_samples,
               color = pass_filter)) +
      
      geom_point(alpha = 0.5) +
      scale_x_log10() +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold")
      ) +
      labs(
        title = "Mean Count vs Non-Zero Samples",
        x = "Mean normalized count",
        y = "Number of non-zero samples",
        color = "Included after filter"
      )
  })
  
  output$pca_plot <- renderPlot({
    filtered <- filtered_gene_stats() %>%
      filter(pass_filter)
    
    validate(
      need(nrow(filtered) >= 2, "Move the filters lower to include at least 2 genes for PCA.")
    )
    
    filtered_counts <- counts %>%
      filter(gene_id %in% filtered$gene_id)
    
    mat <- as.matrix(filtered_counts[, -1])
    
    pca <- prcomp(t(mat), scale. = TRUE)
    
    pc1_var <- round(100 * summary(pca)$importance[2,1], 1)
    pc2_var <- round(100 * summary(pca)$importance[2,2], 1)
    
    pca_df <- data.frame(
      PC1 = pca$x[,1],
      PC2 = pca$x[,2],
      condition = sample_info$condition
    )
    
    ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = 3) +
      scale_color_manual(values = condition_colors) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold")
      ) +
      labs(
        title = "Principal Component Analysis of Filtered Counts",
        x = paste0("PC1 (", pc1_var, "% variance)"),
        y = paste0("PC2 (", pc2_var, "% variance)")
      )
  })
  
  output$heatmap_plot <- renderPlot({
    filtered <- filtered_gene_stats() %>%
      filter(pass_filter) %>%
      arrange(desc(variance)) %>%
      slice_head(n = 50)
    
    validate(
      need(nrow(filtered) >= 2, "Move the filters lower to include at least 2 genes for the heatmap.")
    )
    
    filtered_counts <- counts %>%
      filter(gene_id %in% filtered$gene_id)
    
    mat <- as.matrix(filtered_counts[, -1])
    rownames(mat) <- filtered_counts$gene_id
    
    pheatmap(
      log2(mat + 1),
      scale = "row",
      show_rownames = FALSE,
      show_colnames = FALSE,
      main = "Top Variable Genes Heatmap"
    )
  })
  
  # Differential expression tab
  filtered_de <- reactive({
    if (input$gene_search == "") {
      de
    } else {
      de %>%
        filter(
          grepl(input$gene_search, gene_id, ignore.case = TRUE) |
            grepl(input$gene_search, symbol, ignore.case = TRUE)
        )
    }
  })
  
  output$de_table <- renderDT({
    results <- filtered_de()
    
    validate(
      need(nrow(results) > 0, "No genes matched your search. Try a different gene symbol or Ensembl ID.")
    )
    
    datatable(results)
  })
  
  output$volcano_plot <- renderPlot({
    de$significant <- de$padj < 0.05 & abs(de$log2FoldChange) >= 1
    
    ggplot(de,
           aes(x = log2FoldChange,
               y = -log10(padj),
               color = significant)) +
      geom_point(alpha = 0.6) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold")
      ) +
      labs(
        title = "Volcano Plot of Differential Expression",
        subtitle = "Dashed lines show |log2 fold change| ≥ 1 and adjusted p-value < 0.05 thresholds",
        x = "log2 Fold Change",
        y = "-log10 adjusted p-value",
        color = "Significant"
      )
  })
  
  # Gene expression tab
  output$gene_boxplot <- renderPlot({
    validate(
      need(input$selected_gene != "", "Please select a gene to view its expression.")
    )
    
    selected_counts <- counts %>%
      filter(gene_id == input$selected_gene) %>%
      pivot_longer(
        cols = -gene_id,
        names_to = "sample_id",
        values_to = "expression"
      ) %>%
      left_join(sample_info, by = "sample_id")
    
    validate(
      need(nrow(selected_counts) > 0, "Gene not found in the normalized counts matrix."),
      need(any(!is.na(selected_counts$expression)), "No expression values were found for this gene.")
    )
    
    gene_label <- gene_lookup %>%
      filter(gene_id == input$selected_gene) %>%
      pull(symbol)
    
    if (length(gene_label) == 0 || is.na(gene_label) || gene_label == "") {
      gene_label <- input$selected_gene
    }
    
    ggplot(selected_counts, aes(x = condition, y = expression, fill = condition)) +
      geom_violin(alpha = 0.4, trim = TRUE) +
      geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
      scale_fill_manual(values = condition_colors) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "none"
      ) +
      labs(
        title = paste0("Normalized Expression of ", gene_label),
        subtitle = input$selected_gene,
        x = "Condition",
        y = "Normalized expression"
      )
  })
}

shinyApp(ui = ui, server = server)