
# FUNCTION: volcanoPlot() ####
# This function generates volcano plots and from the output of differential expression analysis.
# input = limma/TopTable output, title_1/2 = string for the groups in order that they appear in the contrast (i.e. lvi vs stroma), batch = string
volcanoPlot <- function(input, title_1, title_2, batch, gene = 10, logFC_threshold = 1, pval_threshold = 1.3, upcolor = "firebrick3", downcolor = "steelblue3", directory) {
  # Create output directory if it does not exist
  if (!file.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  # Preprocess the data
  data <- input %>%
    mutate(neg.log.padj = -1 * log10(adj.P.Val)) %>% 
    filter(!is.na(neg.log.padj)) %>%
    mutate(topstat = neg.log.padj * abs(logFC * 2))
  
  # Filter data into Uppers, Mids, and Lowers based on adjustable thresholds
  Uppers <- data %>%
    filter((logFC >= logFC_threshold) & (neg.log.padj >= pval_threshold)) %>%
    mutate(Group = 'Uppers')
  
  Mids <- data %>%
    filter(((logFC < logFC_threshold) & (logFC > -logFC_threshold)) | (neg.log.padj < pval_threshold)) %>%
    mutate(Group = 'Mids')
  
  Lowers <- data %>%
    filter((logFC <= -logFC_threshold) & (neg.log.padj >= pval_threshold)) %>%
    mutate(Group = 'Lowers')
  if (nrow(Uppers) == 0 & nrow(Lowers) == 0) {
    message("No significant genes found for the given thresholds.")
    
    # Write empty DEGs and GSEA tables
    write.table(data, file = paste0(directory, "/", batch, "_", title_1, "vs", title_2, "_all_genes.csv"),
                row.names = FALSE, col.names = TRUE, sep = ',')
    write.table(data.frame(Genes = character(), logFC = numeric(), AveExpr = numeric(),
                           P.Value = numeric(), adj.P.Val = numeric(), neg.log.padj = numeric(),
                           Group = character()),
                file = paste0(directory, "/", batch, "_", title_1, "vs", title_2, "_DEGS.csv"),
                row.names = FALSE, col.names = TRUE, sep = ',')
    write.table(data.frame(Genes = character(), statistic = numeric()),
                file = paste0(directory, "/", batch, "_", title_1, "vs", title_2, "_gseadegs.csv"),
                row.names = FALSE, col.names = TRUE, sep = ',')
    
    # Create an empty volcano plot
    empty_plot <- ggplot() +
      ggtitle(paste0(title_2, " vs. ", title_1)) +
      theme_void() +
      annotate("text", x = 0, y = 0, label = "No significant genes found", size = 6, color = "red")
    
    # Save the empty plot
    ggsave(file = paste0(directory, "/", batch, "_", title_1, "vs", title_2, "_Batch_Corrected.png"),
           plot = empty_plot,
           width = 11,
           height = 8,
           units = "in",
           dpi = 300)
    
    return(empty_plot)  # Return the empty plot to avoid further processing
  }
  # Create the topstat labels for volcano plot
  inUpperstopstat <- character()
  inLowerstopstat <- character()
  inMidstopstat <- character()
  
  # Determine top genes for Uppers and Lowers
  if (nrow(Uppers) > 0) {
    for (x in 1:nrow(Uppers)) {
      if (round(Uppers$topstat[x]) %in% round(sort(Uppers$topstat, decreasing = TRUE)[1:gene])) {
        inUpperstopstat <- append(inUpperstopstat, TRUE)
      } else {
        inUpperstopstat <- append(inUpperstopstat, FALSE)
      }
    }
  } else {
    inUpperstopstat <- logical(0)  # Empty logical vector if no rows in Uppers
  }
  
  for(x in 1:nrow(Mids)) {
    inMidstopstat <- append(inMidstopstat, FALSE)
  }
  
  if (nrow(Lowers) > 0 && gene <= nrow(Lowers)) {
    for(x in 1:nrow(Lowers)) {
      if(round(Lowers$topstat[x]) %in% round(sort(Lowers$topstat, decreasing = TRUE)[1:gene])) {
        inLowerstopstat <- append(inLowerstopstat, TRUE)
      } else {
        inLowerstopstat <- append(inLowerstopstat, FALSE)
      }
    }
  }
  
  # Combine groups and topstat
  midupp <- append(inUpperstopstat, inMidstopstat)
  tops <- append(midupp, inLowerstopstat)
  Volcano_groups <- rbind(Uppers, Mids)
  Volcano_groups2 <- rbind(Volcano_groups, Lowers)
  Volcano_groups3 <- cbind(Volcano_groups2, tops)
  setDT(Volcano_groups3, keep.rownames = TRUE)[] # Keep row names
  
  # Write output files for all genes and DEGs
  write.table(Volcano_groups3, file = paste0(directory, "/", batch, "_", title_1, "vs", title_2, "_all_genes.csv"), row.names = FALSE, col.names = TRUE, sep = ',')
  
  # DEGs processing and saving
  degs <- rbind(Uppers, Lowers) %>%
    mutate(statistic = (logFC / abs(logFC)) * neg.log.padj)
  
  degs_ordered <- degs[order(degs$statistic, decreasing = TRUE),]
  setDT(degs_ordered, keep.rownames = TRUE)[]
  
  degsforgsea <- degs_ordered %>% dplyr::select(c('rn', 'statistic'))
  degslist <- degs_ordered %>% dplyr::select(c('rn', 'logFC', 'AveExpr', 'P.Value', 'adj.P.Val', 'neg.log.padj', 'Group'))
  colnames(degslist)[1] <- "Genes"
  
  # Write DEGs and GSEA files
  write.table(degsforgsea, file = paste0(directory, "/", batch, "_", title_1, "vs", title_2, "_gseadegs.csv"), row.names = FALSE, col.names = TRUE, sep = ',')
  write.table(degslist, file = paste0(directory, "/", batch, "_", title_1, "vs", title_2, "_DEGS.csv"), row.names = FALSE, col.names = TRUE, sep = ',')
  
  # Generate volcano plot
  volcano_plot <- ggplot(Volcano_groups3, mapping = aes(x = logFC, y = neg.log.padj, label = rn)) +
    geom_point(mapping = aes(color = Group), size = 1) +
    xlab(expression(paste(log[2], " fold change"))) +
    ylab(expression(paste(-log[10], "(q)"))) +
    geom_text_repel(data = . %>%
                      mutate(label = ifelse(Group %in% c("Uppers", "Lowers") & tops == TRUE, rn, "")),
                    aes(label = label), 
                    show.legend = FALSE,
                    max.overlaps = Inf,
                    force = 2,
                    box.padding = 0.2) +
    ggtitle(paste0(title_2, " vs. ", title_1)) +
    scale_color_manual(name = "Groups",
                       labels = c(if(nrow(Lowers) > 0) paste0("Increased in ", title_2, ": ", nrow(Lowers)) else paste0("Unchanged: ", nrow(Mids)),
                                  if(nrow(Lowers) > 0) paste0("Unchanged: ", nrow(Mids)) else paste0("Increased in ", title_1, ": ", nrow(Uppers)),
                                  paste0("Increased in ", title_1, ": ", nrow(Uppers))),
                       values = c(if(nrow(Lowers) > 0) downcolor else "#CBCBCB",
                                  if(nrow(Lowers) > 0) "#CBCBCB" else upcolor,
                                  upcolor)) +
    theme_bw() +
    theme(text = element_text(family = "Arial", color = "#2C2C2E"),
          plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
          legend.title = element_blank(),
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = c(1.02, 0.5),
          legend.justification = c(0, 1),
          legend.key.width = unit(1, "lines"),
          legend.key.height = unit(1, "lines"),
          plot.margin = unit(c(1, 13, 0.5, 0.5), "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = pval_threshold, color = "azure4", linetype = "dashed") +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), color = "azure4", linetype = "dashed")
  
  # Save volcano plot
  ggsave(file = paste0(directory, "/", batch, "_", title_1, "vs", title_2, "_Batch_Corrected.png"),
         width = 11,
         height = 8,
         units = "in",
         dpi = 300)
  
  return(volcano_plot)
}



# simple volcano plot
plot_DE <- function(data, title_1, title_2) {
  plot <- data %>%
    mutate(DE = case_when(
      logFC > 0 & adj.P.Val < 0.05 & logFC > 0.5 ~ "UP",
      logFC < 0 & adj.P.Val < 0.05 & logFC < -0.5 ~ "DOWN",
      TRUE ~ "NOT DE"
    )) %>%
    ggplot(aes(logFC, adj.P.Val, col = DE)) + 
    geom_point(shape = 19, size = 1) + 
    geom_text_repel(data = data %>% 
                      mutate(DE = ifelse(logFC > 0 & adj.P.Val < 0.05, "UP", 
                                         ifelse(logFC < 0 & adj.P.Val < 0.05, "DOWN", "NOT DE"))) %>%
                      rownames_to_column() %>%
                      filter(adj.P.Val < 0.05) %>%
                      top_n(-300, adj.P.Val), aes(label = rowname), show.legend = FALSE,
                    max.overlaps = Inf) +
    theme_bw() +
    geom_vline(xintercept = c(-0.5, 0.5), color = "azure4", linetype = "dashed") +
    geom_hline(yintercept = 0.05, color = "azure4", linetype = "dashed") +
    xlab("Average log-expression") +
    ylab("P-value") +
    ggtitle(paste0(title_2, " vs ", title_1)) +
    scale_color_manual(values = c("steelblue3", "gray", "firebrick3"),
                       labels = c(paste0("Increased in ", title_2), "Unchanged", paste0("Increased in ", title_1)),
                       name = "Expression") +
    scale_y_continuous(trans = trans_reverser('log10'), breaks = c(1e0, 1e-5, 1e-10, 1e-15, 1e-20)) +
    theme(text = element_text(size = 15))
  
  return(plot)
}

plot_DE_list_v1 <- function(data, title_1, title_2, gene_list) {
  plot <- data %>%
    mutate(DE = case_when(
      logFC > 0 & adj.P.Val < 0.05 & logFC > 0.5 ~ "UP",
      logFC < 0 & adj.P.Val < 0.05 & logFC < -0.5 ~ "DOWN",
      TRUE ~ "NOT DE"
    )) %>%
    ggplot(aes(logFC, adj.P.Val, col = DE)) + 
    geom_point(shape = 19, size = 1) + 
    
    # Label only genes in gene_list
    geom_text_repel(data = data %>%
                      rownames_to_column() %>%
                      filter(rowname %in% gene_list), 
                    aes(label = rowname), 
                    show.legend = FALSE,
                    max.overlaps = Inf) +
    
    theme_bw() +
    geom_vline(xintercept = c(-0.5, 0.5), color = "azure4", linetype = "dashed") +
    geom_hline(yintercept = 0.05, color = "azure4", linetype = "dashed") +
    xlab("Average log-expression") +
    ylab("P-value") +
    ggtitle(paste0(title_2, " vs ", title_1)) +
    scale_color_manual(values = c("steelblue3", "gray", "firebrick3"),
                       labels = c(paste0("Increased in ", title_2), "Unchanged", paste0("Increased in ", title_1)),
                       name = "Expression") +
    scale_y_continuous(trans = trans_reverser('log10'), breaks = c(1e0, 1e-5, 1e-10, 1e-15, 1e-20)) +
    theme(text = element_text(size = 15))
  
  return(plot)
}

plot_DE_list <- function(data, title_1, title_2, gene_list) {
  # Create a column to mark genes in gene_list
  data <- data %>%
    rownames_to_column() %>%
    mutate(label_flag = ifelse(rowname %in% gene_list, "Labeled", "Not Labeled")) %>%
    mutate(DE = case_when(
      logFC > 0 & adj.P.Val < 0.05 & logFC > 0.5 ~ "UP",
      logFC < 0 & adj.P.Val < 0.05 & logFC < -0.5 ~ "DOWN",
      TRUE ~ "NOT DE"
    ))
  
  # Create the plot
  plot <- data %>%
    ggplot(aes(logFC, adj.P.Val, col = DE)) + 
    geom_point(shape = 19, size = 1) + 
    
    # Overlay points for labeled genes with a new color (e.g., purple)
    geom_point(data = data %>% filter(label_flag == "Labeled"),
               aes(logFC, adj.P.Val), color = "#800080", size = 1.2, shape = 19) +
    
    # Label only genes in gene_list
    geom_text_repel(data = data %>%
                      filter(rowname %in% gene_list), 
                    aes(label = rowname), 
                    show.legend = FALSE,
                    color = "black",
                    box.padding   = 0.2,
                    point.padding = 0.5,
                    max.overlaps = Inf) +
    
    theme_bw() +
    geom_vline(xintercept = c(-0.5, 0.5), color = "azure4", linetype = "dashed") +
    geom_hline(yintercept = 0.05, color = "azure4", linetype = "dashed") +
    xlab("Average log-expression") +
    ylab("P-value") +
    #    ggtitle(paste0(title_2, " vs ", title_1)) +
    ggtitle(paste0(title_1, "/", title_2)) +
    
    # Original DE colors
    scale_color_manual(values = c("UP" = "firebrick3", "NOT DE" = "gray", "DOWN" = "steelblue3"), 
                       labels = c(paste0("Increased in ", title_2), "Unchanged", paste0("Increased in ", title_1)), 
                       name = "Expression") +
    
    scale_y_continuous(trans = trans_reverser('log10'), breaks = c(0.05, 0.01)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),   
      axis.text.x = element_text(size = 12),   
      axis.text.y = element_text(size = 12),   
      axis.title = element_text(size = 14),
      legend.position = "none")
  
  return(plot)
}