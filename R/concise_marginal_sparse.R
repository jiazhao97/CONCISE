#' CONCISE ligand-receptor interaction analysis: Marginal fitting with sparse matrix for efficient computation
#'
#' @export
CONCISE_marginal_sparse <- function(
    kernel_choice,
    result_path,
    pairdb,
    rawcount,
    loc,
    lib_size = NULL) {

  K_choice <- kernel_choice

  ## calculate libarary size
  if (is.null(lib_size)) {
    lib_size <- as.vector(apply(rawcount, 2, sum))
  }
  n <- length(lib_size)

  ## spatial kernels
  cat("Creating spatial kernels...\n")
  cat(paste0("\tmarginal kernel ", K_choice, "\n"))
  D <- dist(loc, method = "euclidean")
  D <- as.matrix(D)
  lmin <- min(D[D > 0])
  lmax <- max(D[D > 0])
  lval.list <- 10^(seq(log10(lmin),log10(lmax),length.out=10))[5:7]
  lval <- lval.list[K_choice]
  W <- exp(-D^2/(2*(lval^2)))
  rm(D)
  cat("Done!\n")

  ## get ground truth library size to make the scale of normalized gene expressions roughly having sd = 1
  ligand_list <- unique(pairdb$ligand_organized_qc)
  ligand_list <- unlist(lapply(ligand_list, function(t){unlist(strsplit(t, "+", fixed = TRUE))}))
  receptor_list <- unique(pairdb$receptor_organized_qc)
  receptor_list <- unlist(lapply(receptor_list, function(t){unlist(strsplit(t, "+", fixed = TRUE))}))
  genes_use <- union(unique(ligand_list), unique(receptor_list))

  rawcount_adj_lbs <- t(t(rawcount) / lib_size)
  sd_exp_level <- apply(rawcount_adj_lbs[genes_use, ], 1, sd)
  aver_sd <- mean(sd_exp_level)
  gt_s <- lib_size * aver_sd


  ### pre-calculation for parameter estimation
  cat("Pre-calculation for parameter estimation...\n")
  stopifnot(dim(W)[1] == n)
  sl <- sqrt(n)
  I_s_dense <- diag((gt_s^2) / sl)
  I_s <- Matrix(I_s_dense, sparse = TRUE)
  trsqI_s <- sum(I_s * I_s)
  tcrgt_s <- tcrossprod(gt_s)
  K_s <- W / sl * tcrgt_s
  rm(W)
  trsqK_s <- sum(K_s * K_s)
  trK_sI_s <- sum(K_s * I_s_dense)

  S <- matrix(c(trsqK_s, trK_sI_s, trK_sI_s, trsqI_s), nrow = 2, ncol = 2)
  cat("Done!\n")


  ### parameter estimation for ligands
  cat("Parameter estimation for ligands...\n")
  ligand_list <- unique(pairdb$ligand_organized_qc)
  ligand_res <- data.frame(ligand_organized_qc = ligand_list)
  ligand_res$lval <- lval
  ligand_res$aver_sd <- aver_sd
  ligand_res$a_hat <- 0
  ligand_res$sqsigma_hat <- 0
  ligand_res$sqsigma_e_hat <- 0
  ligand_res$valid <- "True"
  gen_y <- matrix(0, nrow = n, ncol = length(ligand_list))
  for (r in 1:length(ligand_list)) {
    t_name <- ligand_list[r]
    t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))

    gen_y_tmp <- rawcount[t_name, , drop = FALSE]
    gen_y_tmp <- as.vector(apply(gen_y_tmp, 2, sum))
    gen_y[, r] <- gen_y_tmp
  }
  # first order
  a_hat <- as.vector(colSums(gen_y * gt_s)) / sum(gt_s^2)
  ligand_res$a_hat <- a_hat
  # second order
  mu_s <- gt_s %*% t(a_hat) #cell by gene
  y_v <- gen_y - mu_s
  q <- matrix(0, nrow = 2, ncol = length(ligand_list))
  q[1, ] <- colSums(y_v * (K_s %*% y_v)) - colSums(diag(K_s) * mu_s)
  q[2, ] <- colSums(y_v * (I_s %*% y_v)) - colSums(diag(I_s) * mu_s)
  q <- q / sl
  # solve theta
  mean_q <- mean(q)
  S_adj <- S / mean_q
  q_adj <- q / mean_q
  theta_hat <- solve(S_adj) %*% q_adj
  sqsigma_hat <- theta_hat[1, ]
  sqsigma_e_hat <- theta_hat[2, ]
  ligand_res$sqsigma_hat <- sqsigma_hat
  ligand_res$sqsigma_e_hat <- sqsigma_e_hat
  # correction
  ligand_res[(sqsigma_hat + sqsigma_e_hat) < 0, "valid"] <- "False"
  idx_sqsigma_hat <- (sqsigma_hat < 0)
  ligand_res[idx_sqsigma_hat, "sqsigma_hat"] <- 0.
  ligand_res[idx_sqsigma_hat, "sqsigma_e_hat"] <- q_adj[2, idx_sqsigma_hat] / S_adj[2,2]
  ligand_res[idx_sqsigma_hat, "valid"] <- "corrected"
  idx_sqsigma_e_hat <- (sqsigma_e_hat < 0)
  ligand_res[idx_sqsigma_e_hat, "sqsigma_e_hat"] <- 0.
  ligand_res[idx_sqsigma_e_hat, "sqsigma_hat"] <- q_adj[1, idx_sqsigma_e_hat] / S_adj[1,1]
  ligand_res[idx_sqsigma_e_hat, "valid"] <- "corrected"
  # fitting score
  sqsigma_hat <- ligand_res$sqsigma_hat
  sqsigma_e_hat <- ligand_res$sqsigma_e_hat
  F_score <- (sqsigma_hat^2)*trsqK_s + (sqsigma_e_hat^2)*trsqI_s + 2*sqsigma_hat*sqsigma_e_hat*trK_sI_s - 2*sqsigma_hat*q[1, ] - 2*sqsigma_e_hat*q[2, ]
  ligand_res$F_score <- F_score
  cat("Done!\n")


  ### parameter estimation for receptors
  cat("Parameter estimation for receptors...\n")
  receptor_list <- unique(pairdb$receptor_organized_qc)
  receptor_res <- data.frame(receptor_organized_qc = receptor_list)
  receptor_res$lval <- lval
  receptor_res$aver_sd <- aver_sd
  receptor_res$a_hat <- 0
  receptor_res$sqsigma_hat <- 0
  receptor_res$sqsigma_e_hat <- 0
  receptor_res$valid <- "True"
  gen_y <- matrix(0, nrow = n, ncol = length(receptor_list))
  for (r in 1:length(receptor_list)) {
    t_name <- receptor_list[r]
    t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))

    gen_y_tmp <- rawcount[t_name, , drop = FALSE]
    gen_y_tmp <- as.vector(apply(gen_y_tmp, 2, sum))
    gen_y[, r] <- gen_y_tmp
  }
  # first order
  a_hat <- as.vector(colSums(gen_y * gt_s)) / sum(gt_s^2)
  receptor_res$a_hat <- a_hat
  # second order
  mu_s <- gt_s %*% t(a_hat) #cell by gene
  y_v <- gen_y - mu_s
  q <- matrix(0, nrow = 2, ncol = length(receptor_list))
  q[1, ] <- colSums(y_v * (K_s %*% y_v)) - colSums(diag(K_s) * mu_s)
  q[2, ] <- colSums(y_v * (I_s %*% y_v)) - colSums(diag(I_s) * mu_s)
  q <- q / sl
  # solve theta
  mean_q <- mean(q)
  S_adj <- S / mean_q
  q_adj <- q / mean_q
  theta_hat <- solve(S_adj) %*% q_adj
  sqsigma_hat <- theta_hat[1, ]
  sqsigma_e_hat <- theta_hat[2, ]
  receptor_res$sqsigma_hat <- sqsigma_hat
  receptor_res$sqsigma_e_hat <- sqsigma_e_hat
  # correction
  receptor_res[(sqsigma_hat + sqsigma_e_hat) < 0, "valid"] <- "False"
  idx_sqsigma_hat <- (sqsigma_hat < 0)
  receptor_res[idx_sqsigma_hat, "sqsigma_hat"] <- 0.
  receptor_res[idx_sqsigma_hat, "sqsigma_e_hat"] <- q_adj[2, idx_sqsigma_hat] / S_adj[2,2]
  receptor_res[idx_sqsigma_hat, "valid"] <- "corrected"
  idx_sqsigma_e_hat <- (sqsigma_e_hat < 0)
  receptor_res[idx_sqsigma_e_hat, "sqsigma_e_hat"] <- 0.
  receptor_res[idx_sqsigma_e_hat, "sqsigma_hat"] <- q_adj[1, idx_sqsigma_e_hat] / S_adj[1,1]
  receptor_res[idx_sqsigma_e_hat, "valid"] <- "corrected"
  # fitting score
  sqsigma_hat <- receptor_res$sqsigma_hat
  sqsigma_e_hat <- receptor_res$sqsigma_e_hat
  F_score <- (sqsigma_hat^2)*trsqK_s + (sqsigma_e_hat^2)*trsqI_s + 2*sqsigma_hat*sqsigma_e_hat*trK_sI_s - 2*sqsigma_hat*q[1, ] - 2*sqsigma_e_hat*q[2, ]
  receptor_res$F_score <- F_score
  cat("Done!\n")


  ### save results
  write.csv(ligand_res, paste0(result_path, "/MoM_marginal_fit_ligand_res_Kchoice_", K_choice, ".csv", sep = ""))
  write.csv(receptor_res, paste0(result_path, "/MoM_marginal_fit_receptor_res_Kchoice_", K_choice, ".csv", sep = ""))

  return(list(ligand_res, receptor_res))
}

