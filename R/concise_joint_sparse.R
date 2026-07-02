#' CONCISE ligand-receptor interaction analysis: Joint ligand-receptor fitting with sparse matrix for efficient computation
#'
#' @export
CONCISE_joint_sparse <- function(
    result_path,
    cell_type_anno,
    ligand_cell_type,
    receptor_cell_type,
    pairdb,
    rawcount,
    loc,
    Wx_distance_threshold,
    remove_invalid_genes = TRUE,
    one_side_test = TRUE,
    lib_size = NULL) {

  ### load and summarize marginal fitting results
  cat("Loading and summarizing marginal fitting results for ligands and receptors...\n")
  ## load ligand results
  ligand_res <- read.csv(paste0(result_path, "/MoM_marginal_fit_ligand_res_Kchoice_1.csv"), row.names = 1)
  colnames(ligand_res)[2:length(colnames(ligand_res))] <- paste(colnames(ligand_res)[2:length(colnames(ligand_res))], "1", sep = ".")
  ligand_res.1 <- ligand_res
  ligand_res <- read.csv(paste0(result_path, "/MoM_marginal_fit_ligand_res_Kchoice_2.csv"), row.names = 1)
  colnames(ligand_res)[2:length(colnames(ligand_res))] <- paste(colnames(ligand_res)[2:length(colnames(ligand_res))], "2", sep = ".")
  ligand_res.2 <- ligand_res
  ligand_res <- read.csv(paste0(result_path, "/MoM_marginal_fit_ligand_res_Kchoice_3.csv"), row.names = 1)
  colnames(ligand_res)[2:length(colnames(ligand_res))] <- paste(colnames(ligand_res)[2:length(colnames(ligand_res))], "3", sep = ".")
  ligand_res.3 <- ligand_res
  ligand_res <- merge(ligand_res.1, ligand_res.2, by = "ligand_organized_qc")
  ligand_res <- merge(ligand_res, ligand_res.3, by = "ligand_organized_qc")
  rm(ligand_res.1, ligand_res.2, ligand_res.3)

  ## load receptor results
  receptor_res <- read.csv(paste0(result_path, "/MoM_marginal_fit_receptor_res_Kchoice_1.csv"), row.names = 1)
  colnames(receptor_res)[2:length(colnames(receptor_res))] <- paste(colnames(receptor_res)[2:length(colnames(receptor_res))], "1", sep = ".")
  receptor_res.1 <- receptor_res
  receptor_res <- read.csv(paste0(result_path, "/MoM_marginal_fit_receptor_res_Kchoice_2.csv"), row.names = 1)
  colnames(receptor_res)[2:length(colnames(receptor_res))] <- paste(colnames(receptor_res)[2:length(colnames(receptor_res))], "2", sep = ".")
  receptor_res.2 <- receptor_res
  receptor_res <- read.csv(paste0(result_path, "/MoM_marginal_fit_receptor_res_Kchoice_3.csv"), row.names = 1)
  colnames(receptor_res)[2:length(colnames(receptor_res))] <- paste(colnames(receptor_res)[2:length(colnames(receptor_res))], "3", sep = ".")
  receptor_res.3 <- receptor_res
  receptor_res <- merge(receptor_res.1, receptor_res.2, by = "receptor_organized_qc")
  receptor_res <- merge(receptor_res, receptor_res.3, by = "receptor_organized_qc")
  rm(receptor_res.1, receptor_res.2, receptor_res.3)

  ## summarize the single variant MoM result and organize the cellchatdb accordingly
  n_K <- 3
  # select kernels
  ligand_res$select.kernel <- apply(ligand_res[paste("F_score.", 1:n_K, sep = "")], 1, which.min)
  ligand_res$a_hat <- as.vector(ligand_res$a_hat.1)
  ligand_res$sqsigma_hat <- 0
  ligand_res$sqsigma_e_hat <- 0
  for (r in 1:dim(ligand_res)[1]) {
    idx_kernel <- ligand_res$select.kernel[r]
    ligand_res$sqsigma_hat[r] <- ligand_res[r, paste0("sqsigma_hat.", idx_kernel)]
    ligand_res$sqsigma_e_hat[r] <- ligand_res[r, paste0("sqsigma_e_hat.", idx_kernel)]
  }

  receptor_res$select.kernel <- apply(receptor_res[paste("F_score.", 1:n_K, sep = "")], 1, which.min)
  receptor_res$a_hat <- as.vector(receptor_res$a_hat.1)
  receptor_res$sqsigma_hat <- 0
  receptor_res$sqsigma_e_hat <- 0
  for (r in 1:dim(receptor_res)[1]) {
    idx_kernel <- receptor_res$select.kernel[r]
    receptor_res$sqsigma_hat[r] <- receptor_res[r, paste0("sqsigma_hat.", idx_kernel)]
    receptor_res$sqsigma_e_hat[r] <- receptor_res[r, paste0("sqsigma_e_hat.", idx_kernel)]
  }

  # remove invalid genes
  if (remove_invalid_genes) {
    for (i_k in 1:n_K) {
      ligand_res <- ligand_res[ligand_res[paste0("valid.", i_k)] != "False", ]
      receptor_res <- receptor_res[receptor_res[paste0("valid.", i_k)] != "False", ]
    }
  }
  rownames(ligand_res) <- as.vector(ligand_res$ligand_organized_qc)
  rownames(receptor_res) <- as.vector(receptor_res$receptor_organized_qc)

  pairdb <- pairdb[(pairdb$ligand_organized_qc %in% ligand_res$ligand_organized_qc) & (pairdb$receptor_organized_qc %in% receptor_res$receptor_organized_qc), ]

  ligand_res_tmp <- ligand_res[pairdb$ligand_organized_qc, ]
  stopifnot(sum(ligand_res_tmp$ligand_organized_qc != pairdb$ligand_organized_qc) == 0)
  pairdb$a1_hat <- ligand_res_tmp$a_hat
  pairdb$sqsigma1_hat <- ligand_res_tmp$sqsigma_hat
  pairdb$sqsigma_e1_hat <- ligand_res_tmp$sqsigma_e_hat
  pairdb$select_kernel1 <- ligand_res_tmp$select.kernel

  receptor_res_tmp <- receptor_res[pairdb$receptor_organized_qc, ]
  stopifnot(sum(receptor_res_tmp$receptor_organized_qc != pairdb$receptor_organized_qc) == 0)
  pairdb$a2_hat <- receptor_res_tmp$a_hat
  pairdb$sqsigma2_hat <- receptor_res_tmp$sqsigma_hat
  pairdb$sqsigma_e2_hat <- receptor_res_tmp$sqsigma_e_hat
  pairdb$select_kernel2 <- receptor_res_tmp$select.kernel

  pairdb$kernel_type <- paste(pairdb$select_kernel1, pairdb$select_kernel2, sep = "-")
  cat("Done!\n")



  ### MoM parameter estimation for delta
  celltype_anno <- cell_type_anno
  L_celltype <- ligand_cell_type
  R_celltype <- receptor_cell_type
  index_cell <- which(celltype_anno %in% c(L_celltype, R_celltype))

  ### pre-screening expression level
  cat("Pre-screening expression level\n")
  if (is.null(lib_size)) {
    lib_size <- as.vector(apply(rawcount, 2, sum))
  }
  aver_sd <- ligand_res$aver_sd.1[1]
  gt_s <- lib_size * aver_sd

  #Lgene
  cat("\tLigands\n")
  pairdb$exp_val_Lgene_L <- 0.
  pairdb$exp_val_Lgene_R <- 0.
  for (i in 1:dim(pairdb)[1]) {
    t_name <- pairdb$ligand_organized_qc[i]
    t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))
    exp_val <- rawcount[t_name, , drop = FALSE]
    exp_val <- as.vector(apply(exp_val, 2, sum))
    exp_val <- exp_val - pairdb$a1_hat[i] * gt_s

    stopifnot(length(celltype_anno) == dim(rawcount)[2])
    exp_val_L <- mean(exp_val[celltype_anno == L_celltype])
    exp_val_R <- mean(exp_val[celltype_anno == R_celltype])
    pairdb$exp_val_Lgene_L[i] <- exp_val_L
    pairdb$exp_val_Lgene_R[i] <- exp_val_R
  }

  #Rgene
  cat("\tReceptors\n")
  pairdb$exp_val_Rgene_L <- 0.
  pairdb$exp_val_Rgene_R <- 0.
  for (i in 1:dim(pairdb)[1]) {
    t_name <- pairdb$receptor_organized_qc[i]
    t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))
    exp_val <- rawcount[t_name, , drop = FALSE]
    exp_val <- as.vector(apply(exp_val, 2, sum))
    exp_val <- exp_val - pairdb$a2_hat[i] * gt_s

    stopifnot(length(celltype_anno) == dim(rawcount)[2])
    exp_val_L <- mean(exp_val[celltype_anno == L_celltype])
    exp_val_R <- mean(exp_val[celltype_anno == R_celltype])
    pairdb$exp_val_Rgene_L[i] <- exp_val_L
    pairdb$exp_val_Rgene_R[i] <- exp_val_R
  }

  ## only include cells of specific cell types
  n <- dim(rawcount)[2]
  sl <- sqrt(n)
  rawcount <- rawcount[, index_cell]
  loc <- loc[index_cell, ]
  aver_sd <- ligand_res$aver_sd.1[1]
  stopifnot(aver_sd == receptor_res$aver_sd.1[1])
  gt_s <- lib_size * aver_sd
  gt_s <- gt_s[index_cell]
  celltype_anno <- as.vector(celltype_anno)[index_cell]
  stopifnot(length(unique(celltype_anno)) == 2)
  cat("Done!\n")

  ## calculate interaction kernel
  cat("Creating interaction kernel...\n")
  lvalx <- Wx_distance_threshold
  D <- dist(loc, method = "euclidean")
  D <- as.matrix(D)
  Wx <- D
  Wx[Wx <= lvalx] <- 1
  Wx[Wx > lvalx] <- 0.
  diag(Wx) <- 0.
  Wx[as.vector(celltype_anno) != L_celltype, ] <- 0
  Wx[, as.vector(celltype_anno) != R_celltype] <- 0
  cat("Done!\n")

  ## pre-calculation for parameter estimation
  cat("Parameter estimation...\n")
  cat("\tPre-calculation for parameter estimation\n")
  n <- length(index_cell)
  stopifnot(length(gt_s) == n)
  stopifnot(dim(rawcount)[2] == n)
  stopifnot(dim(loc)[1] == n)
  I_s_dense <- diag((gt_s^2) / sl)
  I_s <- Matrix(I_s_dense, sparse = TRUE)
  tcrgt_s <- tcrossprod(gt_s)

  Kx_s <- Wx / sl * tcrgt_s # all.equal(Kx * tcrossprod(gt_s), t(t(Kx * gt_s) * gt_s))
  Kx_s <- Matrix(Kx_s, sparse = TRUE)
  crKx_s <- crossprod(Kx_s, Kx_s) #equal to t(Kx_s) %*% Kx_s
  trcrKx_s <- sum(diag(crKx_s))
  rm(Wx)

  ## scan each ligand-receptor pair
  cat("\tScanning each ligand-receptor pair...\n")
  ligand_list_qc <- pairdb$ligand_organized_qc
  gen_y_L <- matrix(0, nrow = n, ncol = length(ligand_list_qc))
  for (r in 1:length(ligand_list_qc)) {
    t_name <- ligand_list_qc[r]
    t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))

    gen_y_tmp <- rawcount[t_name, , drop = FALSE]
    gen_y_tmp <- as.vector(apply(gen_y_tmp, 2, sum))
    gen_y_L[, r] <- gen_y_tmp
  }
  receptor_list_qc <- pairdb$receptor_organized_qc
  gen_y_R <- matrix(0, nrow = n, ncol = length(receptor_list_qc))
  for (r in 1:length(receptor_list_qc)) {
    t_name <- receptor_list_qc[r]
    t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))

    gen_y_tmp <- rawcount[t_name, , drop = FALSE]
    gen_y_tmp <- as.vector(apply(gen_y_tmp, 2, sum))
    gen_y_R[, r] <- gen_y_tmp
  }
  y_v_L <- gen_y_L - gt_s %*% t(pairdb$a1_hat) # cell by gene
  y_v_R <- gen_y_R - gt_s %*% t(pairdb$a2_hat) # cell by gene
  delta_hat_mat <- crossprod(y_v_L, Kx_s %*% y_v_R) / sl / trcrKx_s
  pairdb$delta_hat <- diag(delta_hat_mat)
  cat("Done!\n")



  ### conduct statistical inference for delta
  ## construct kernel list
  cat("Creating spatial kernels...\n")
  lval.list <- c(ligand_res$lval.1[1], ligand_res$lval.2[1], ligand_res$lval.3[1])
  n_K <- length(lval.list)
  K_s_list <- list()
  for (i_k in 1:n_K) {
    cat(paste0("\tSpatial kernel ", i_k, "\n"))
    lval <- lval.list[i_k]
    K_s_list[[i_k]] <- exp(-D^2/(2*(lval^2))) / sl * tcrgt_s ##lval has to be matched
  }
  rm(D)

  ## pre-calculation for statistical inference
  cat("Statistical inference...\n")
  cat("\tPre-calculation for statistical inference\n")
  term_se_names <- c("term_se_1_1", "term_se_1_2",
                     "term_se_2_1", "term_se_2_2", "term_se_2_3", "term_se_2_4", "term_se_2_5",
                     "term_se_3_1", "term_se_3_2", "term_se_3_3", "term_se_3_4")
  stopifnot(length(K_s_list) == n_K)
  term_se_precalculate <- matrix(0, nrow = n_K^2, ncol = length(term_se_names))
  colnames(term_se_precalculate) <- term_se_names
  kernel_type_names <- paste(rep(1:n_K, rep(n_K, n_K)), rep(1:n_K, n_K), sep = "-")
  rownames(term_se_precalculate) <- kernel_type_names

  ## terms not related to marginal kernels
  cat("\t\tTerms not related to marginal kernels\n")
  Kx_s_dot_Kx_s <- Kx_s^2
  term_se_1_1 <- sum(as.matrix(Kx_s_dot_Kx_s) * tcrgt_s)
  term_se_1_2 <- sum(Kx_s_dot_Kx_s * Kx_s)
  term_se_precalculate[, "term_se_1_1"] <- term_se_1_1
  term_se_precalculate[, "term_se_1_2"] <- term_se_1_2
  rm(Kx_s_dot_Kx_s)

  term_se_2_1 <- sum(crKx_s * crKx_s) # crKx_s symmetric
  term_se_precalculate[, "term_se_2_1"] <- term_se_2_1

  Kx_sTI_s <- crossprod(Kx_s, I_s) #check: all.equal(crossprod(Kx_s, I_s), t(Kx_s) %*% I_s)
  Kx_sI_s <- Kx_s %*% I_s
  term_se_2_5 <- sum(diag(Kx_sTI_s %*% Kx_sI_s))
  term_se_precalculate[, "term_se_2_5"] <- term_se_2_5

  s_Kx_s_Kx_s_outer <- crossprod(Kx_s, Kx_s*gt_s) #check: all.equal(crossprod(Kx_s, Kx_s*gt_s), t(Kx_s) %*% (Kx_s * gt_s))
  Kx_sT <- t(Kx_s)
  s_Kx_sT_Kx_sT_outer <- crossprod(Kx_sT, Kx_sT*gt_s) #check: all.equal(crossprod(Kx_sT, Kx_sT*gt_s), Kx_s %*% (t(Kx_s) * gt_s))
  term_se_3_2 <- sum(s_Kx_s_Kx_s_outer * I_s)
  term_se_3_4 <- sum(s_Kx_sT_Kx_sT_outer * I_s)
  term_se_precalculate[, "term_se_3_2"] <- term_se_3_2
  term_se_precalculate[, "term_se_3_4"] <- term_se_3_4
  rm(Kx_sT)

  ## terms related to only one marginal kernel
  cat("\t\tTerms related to only one marginal kernel\n")
  Kx_sTK_s_list <- list()
  Kx_sK_s_list <- list()
  for (i_k in 1:n_K) {
    K_s <- K_s_list[[i_k]]

    Kx_sTK_s_list[[i_k]] <- crossprod(Kx_s, K_s) #check: all.equal(crossprod(Kx_s, K_s), t(Kx_s) %*% K_s)
    Kx_sK_s_list[[i_k]] <- Kx_s %*% K_s
  }

  s_Kx_s_Kx_s_outer <- as.matrix(s_Kx_s_Kx_s_outer)
  s_Kx_sT_Kx_sT_outer <- as.matrix(s_Kx_sT_Kx_sT_outer)
  names_1 <- rep(1:n_K, rep(n_K, n_K))
  names_2 <- rep(1:n_K, n_K)
  for (i_k in 1:n_K) {
    K_s <- K_s_list[[i_k]]
    Kx_sTK_s <- Kx_sTK_s_list[[i_k]]
    Kx_sK_s <- Kx_sK_s_list[[i_k]]

    term_se_2_3 <- sum(t(Kx_sI_s) * Kx_sTK_s)
    term_se_precalculate[which(names_1 == i_k), "term_se_2_3"] <- term_se_2_3

    term_se_2_4 <- sum(t(Kx_sTI_s) * Kx_sK_s)
    term_se_precalculate[which(names_2 == i_k), "term_se_2_4"] <- term_se_2_4

    term_se_3_1 <- sum(s_Kx_s_Kx_s_outer * K_s)
    term_se_precalculate[which(names_2 == i_k), "term_se_3_1"] <- term_se_3_1

    term_se_3_3 <- sum(s_Kx_sT_Kx_sT_outer * K_s)
    term_se_precalculate[which(names_1 == i_k), "term_se_3_3"] <- term_se_3_3
  }
  rm(K_s_list)
  rm(s_Kx_s_Kx_s_outer)
  rm(s_Kx_sT_Kx_sT_outer)

  ## terms related to both two marginal kernels
  cat("\t\tTerms related to both two marginal kernels\n")
  for (i_k1 in 1:n_K) {
    # kernel 1
    Kx_sTK_s <- Kx_sTK_s_list[[i_k1]]

    for (i_k2 in 1:n_K) {
      # kernel 2
      Kx_sK_s <- Kx_sK_s_list[[i_k2]]

      term_se_2_2 <- sum(t(Kx_sTK_s) * Kx_sK_s)
      term_se_precalculate[paste(i_k1, i_k2, sep = "-"), "term_se_2_2"] <- term_se_2_2
    }
  }
  rm(Kx_sTK_s_list)
  rm(Kx_sK_s_list)

  term_se_precalculate[, "term_se_1_1"] <- term_se_precalculate[, "term_se_1_1"] / (sl^2)
  term_se_precalculate[, "term_se_1_2"] <- term_se_precalculate[, "term_se_1_2"] / sl
  term_se_precalculate[, "term_se_3_1"] <- term_se_precalculate[, "term_se_3_1"] / sl
  term_se_precalculate[, "term_se_3_2"] <- term_se_precalculate[, "term_se_3_2"] / sl
  term_se_precalculate[, "term_se_3_3"] <- term_se_precalculate[, "term_se_3_3"] / sl
  term_se_precalculate[, "term_se_3_4"] <- term_se_precalculate[, "term_se_3_4"] / sl

  ## scan each ligand-receptor pair in pairdb
  cat("\tScanning each ligand-receptor pair in pairdb\n")
  n_pair <- dim(pairdb)[1]

  pairdb$se_delta_hat <- 0

  for (r in 1:n_pair) {
    # get parameter estimate
    pairdb_tmp <- pairdb[r, ]

    a1_hat <- pairdb_tmp$a1_hat
    sqsigma1_hat <- pairdb_tmp$sqsigma1_hat
    sqsigma_e1_hat <- pairdb_tmp$sqsigma_e1_hat

    a2_hat <- pairdb_tmp$a2_hat
    sqsigma2_hat <- pairdb_tmp$sqsigma2_hat
    sqsigma_e2_hat <- pairdb_tmp$sqsigma_e2_hat

    delta_hat <- pairdb_tmp$delta_hat
    kernel_type <- pairdb_tmp$kernel_type

    if (delta_hat/sqrt(sqsigma1_hat+sqsigma_e1_hat)/sqrt(sqsigma2_hat+sqsigma_e2_hat) < -1) {
      delta_hat <- sqrt(sqsigma1_hat+sqsigma_e1_hat)*sqrt(sqsigma2_hat+sqsigma_e2_hat)*(-0.9)
    }

    # calculate se
    coef_se_vec <- c(a1_hat*a2_hat, delta_hat,
                     delta_hat^2, sqsigma1_hat*sqsigma2_hat, sqsigma1_hat*sqsigma_e2_hat, sqsigma2_hat*sqsigma_e1_hat, sqsigma_e1_hat*sqsigma_e2_hat,
                     a1_hat*sqsigma2_hat, a1_hat*sqsigma_e2_hat, a2_hat*sqsigma1_hat, a2_hat*sqsigma_e1_hat)
    term_se_vec <- as.vector(term_se_precalculate[kernel_type, ])
    se_delta_hat <- sqrt(sum(coef_se_vec * term_se_vec)) / trcrKx_s

    # store se
    pairdb$se_delta_hat[r] <- se_delta_hat
  }

  if (one_side_test) {
    pairdb$z_score <- pairdb$delta_hat / pairdb$se_delta_hat
    P_value <- pnorm(pairdb$z_score, lower.tail=F)
    pairdb$p_value <- P_value
  } else {
    pairdb$z_score <- pairdb$delta_hat / pairdb$se_delta_hat
    W <- pairdb$z_score ^ 2
    P_value <- pchisq(W, 1, lower.tail=F)
    pairdb$p_value <- P_value
  }
  cat("Done!\n")

  ### save results
  write.csv(pairdb, paste0(result_path, "/MoM_joint_fit_LRI_res_", L_celltype, "_to_", R_celltype, ".csv", sep = ""))

  return(pairdb)
}


