#' Organize ligand-receptor interaction database
#'
#' @export
InteractionDB <- function(
    pairdb,
    rawcount,
    sparsity_threshold) {

  ## filter out ligands and receptors that are extremely sparse in count data
  gene_names <- rownames(rawcount)
  n <- dim(rawcount)[2]
  gene_prop <- rowSums(rawcount > 0) / n
  gene_names_qc <- names(gene_prop)[gene_prop > sparsity_threshold] # genes that pass QC sparsity threshold
  gene_names <- names(gene_prop)[gene_prop > 0.01] # genes that at least have a very small proportion of expressing cells or spots

  target_list <- unique(pairdb$ligand_organized)
  message(paste0("Total number of ligands = ", length(target_list)))
  ligand_qc_original <- c() # record the original combination
  ligand_qc <- c() # some of the genes may not be included in the analyzed data, so the ligand combination names could be changed
  for (i in 1:length(target_list)) {
    t_name <- target_list[i]
    t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))
    n_t_name <- length(t_name)

    if ((n_t_name == 1) && (t_name %in% gene_names_qc)) {
      ligand_qc_original <- c(ligand_qc_original, target_list[i])
      ligand_qc <- c(ligand_qc, target_list[i])
    }

    if (n_t_name > 1) {
      t_name <- t_name[t_name %in% gene_names]
      t_values <- rawcount[t_name, , drop = FALSE]
      t_values <- as.vector(apply(t_values, 2, sum))
      t_prop <- sum(t_values > 0) / n

      if (t_prop > sparsity_threshold) {
        if (length(t_name) == n_t_name) {
          stopifnot(str_c(t_name, collapse = "+") == target_list[i])
          ligand_qc_original <- c(ligand_qc_original, target_list[i])
          ligand_qc <- c(ligand_qc, target_list[i])
        } else if (!(str_c(t_name, collapse = "+") %in% target_list)) {
          stopifnot(length(t_name) < n_t_name)
          ligand_qc_original <- c(ligand_qc_original, target_list[i])
          ligand_qc <- c(ligand_qc, str_c(t_name, collapse = "+"))
        }
      }
    }
  }
  message(paste0("Total number of ligands after sparsity QC = ", length(ligand_qc)))

  target_list <- unique(pairdb$receptor_organized)
  message(paste0("Total number of receptors = ", length(target_list)))
  receptor_qc_original <- c() # record the original combination
  receptor_qc <- c() # some of the genes may not be included in the analyzed data, so the receptor combination names could be changed
  for (i in 1:length(target_list)) {
    t_name <- target_list[i]
    t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))
    n_t_name <- length(t_name)

    if ((n_t_name == 1) && (t_name %in% gene_names_qc)) {
      receptor_qc_original <- c(receptor_qc_original, target_list[i])
      receptor_qc <- c(receptor_qc, target_list[i])
    }

    if (n_t_name > 1) {
      t_name <- t_name[t_name %in% gene_names]
      t_values <- rawcount[t_name, , drop = FALSE]
      t_values <- as.vector(apply(t_values, 2, sum))
      t_prop <- sum(t_values > 0) / n

      if (t_prop > sparsity_threshold) {
        if (length(t_name) == n_t_name) {
          stopifnot(str_c(t_name, collapse = "+") == target_list[i])
          receptor_qc_original <- c(receptor_qc_original, target_list[i])
          receptor_qc <- c(receptor_qc, target_list[i])
        } else if (!(str_c(t_name, collapse = "+") %in% target_list)) {
          stopifnot(length(t_name) < n_t_name)
          receptor_qc_original <- c(receptor_qc_original, target_list[i])
          receptor_qc <- c(receptor_qc, str_c(t_name, collapse = "+"))
        }
      }
    }
  }
  message(paste0("Total number of receptors after sparsity QC = ", length(receptor_qc)))

  pairdb <- pairdb[(pairdb$ligand_organized %in% ligand_qc_original) & (pairdb$receptor_organized %in% receptor_qc_original), ]
  names(ligand_qc) <- ligand_qc_original
  names(receptor_qc) <- receptor_qc_original
  pairdb$ligand_organized_qc <- ligand_qc[pairdb$ligand_organized]
  pairdb$receptor_organized_qc <- receptor_qc[pairdb$receptor_organized]

  return(pairdb)
}
