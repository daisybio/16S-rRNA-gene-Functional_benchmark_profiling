process_file <- function(file_path, map_data, combined_matrices) {
  print(file_path)
  KO_table <- read.table(file_path, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
  KO_table <- KO_table %>% dplyr::select(one_of(dput(as.character(map_data$Sample_ID))))
 KO_table <- OTUtable::filter_taxa(KO_table, abundance = 0, persistence = 5)
 combined_matrices[[length(combined_matrices) + 1]] <- KO_table
 names(combined_matrices)[length(combined_matrices)] <- basename(file_path)
 return(combined_matrices)
}

process_element <- function(j) {
  KO_table <- as.data.frame(t(combined_matrices[[j]]))
  KO_table$Sample_ID <- rownames(KO_table)
  KO_table <- dplyr::left_join(KO_table, map_data[, c(1, coln)])
  rownames(KO_table) <- KO_table$Sample_ID
  KO_table$Sample_ID <- NULL
  KO_table <- reshape2::melt(KO_table)
  colnames(KO_table) <- c("CASE", "KO", "Abundance")
  Wilcox_p1 <- ggpubr::compare_means(Abundance ~ CASE, data = KO_table, group.by = 'KO',
                                     paired = FALSE, p.adjust.method = "BH", method = "wilcox.test")
  Wilcox_p1_sig <- dplyr::filter(Wilcox_p1, p < 0.05)
  colnames(Wilcox_p1) <- paste(colnames(Wilcox_p1), "p1", sep = "_")
  colnames(Wilcox_p1)[1] <- "KO"
  colnames(Wilcox_wgs)[1] <- "KO"
  common_genes <- intersect(Wilcox_wgs$KO, Wilcox_p1$KO)
  Wilcox_wgs_1 <- subset(Wilcox_wgs, KO %in% common_genes)
  Wilcox_p1_1 <- subset(Wilcox_p1, KO %in% common_genes)
  threshold <- 0.05
  tp <- sum(Wilcox_wgs_1$p <= threshold & Wilcox_p1_1$p.format_p1 <= threshold)
  tn <- sum(Wilcox_wgs_1$p > threshold & Wilcox_p1_1$p.format_p1 > threshold)
  fp <- sum(Wilcox_wgs_1$p > threshold & Wilcox_p1_1$p.format_p1 <= threshold)
  fn <- sum(Wilcox_wgs_1$p <= threshold & Wilcox_p1_1$p.format_p1 > threshold)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * precision * recall / (precision + recall)
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  result <- data.frame(
    Data = sub("\\..*", "", names(combined_matrices)[j]),
    precision = precision,
    recall = recall,
    f1 = f1,
    accuracy = accuracy
  )
  return(result)
}
