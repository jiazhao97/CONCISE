#' CONCISE ligand-receptor interaction analysis
#'
#' @export
CONCISE <- function(
    pairdb,
    rawcount,
    loc,
    Wx_distance_threshold,
    remove_invalid_genes = TRUE,
    lib_size = NULL) {

  ## spatial kernels
  message("Creating spatial kernels...")
  D <- pdist(loc, metric = "euclidean", p = 2)
  # create K_p1, K_p2 kernel lists
  lmin <- min(D[D > 0])
  lmax <- max(D[D > 0])
  lval.list <- 10^(seq(log10(lmin),log10(lmax),length.out=10))[4:6]
  n_K <- length(lval.list)
  W_list <- list()
  for (i_k in 1:n_K) {
    message(paste0("\tmarginal kernel ", i_k))
    lval <- lval.list[i_k]
    W_list[[i_k]] <- exp(-D^2/(2*(lval^2)))
    W_list[[i_k]] <- 0.5 * (W_list[[i_k]] + t(W_list[[i_k]]))
  }
  # create Kx
  message("\tinteraction kernel")
  Wx <- D
  Wx[D <= Wx_distance_threshold] <- 1
  Wx[D > Wx_distance_threshold] <- 0
  diag(Wx) <- 0
  message("Done!")

  ## calculate libarary size
  if (is.null(lib_size)) {
    lib_size <- as.vector(apply(rawcount, 2, sum))
  }
  n <- length(lib_size)

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


  ### iterate over choices of kernels and finish all the estimate based on the current kernel
  message("Iterating over choices of kernels and estimating parameters of marginal distributions...")
  stopifnot(length(W_list) == n_K)
  stopifnot(dim(W_list[[1]])[1] == n)
  sl <- sqrt(n)

  I <- diag(n) / sl
  I_s <- t(I * gt_s) * gt_s
  I_s <- 0.5 * (I_s + t(I_s))
  sqI_s <- I_s %*% I_s
  trsqI_s <- sum(diag(sqI_s))

  K_s_list <- list()
  for (i_k in 1:n_K) {
    K <- W_list[[i_k]] / sl
    K_s <- t(K * gt_s) * gt_s
    K_s <- 0.5 * (K_s + t(K_s))
    K_s_list[[i_k]] <- K_s
  }
  rm(W_list)

  ligand_list <- unique(pairdb$ligand_organized_qc)
  ligand_res <- data.frame(ligand_organized_qc = ligand_list)
  ligand_res$a_hat <- 0
  for (i_k in 1:n_K) {
    ligand_res[paste0("sqsigma_hat.", i_k)] <- 0
    ligand_res[paste0("sqsigma_e_hat.", i_k)] <- 0
    ligand_res[paste0("valid.", i_k)] <- "True"
    ligand_res[paste0("F_score.", i_k)] <- 0
  }

  receptor_list <- unique(pairdb$receptor_organized_qc)
  receptor_res <- data.frame(receptor_organized_qc = receptor_list)
  receptor_res$a_hat <- 0
  for (i_k in 1:n_K) {
    receptor_res[paste0("sqsigma_hat.", i_k)] <- 0
    receptor_res[paste0("sqsigma_e_hat.", i_k)] <- 0
    receptor_res[paste0("valid.", i_k)] <- "True"
    receptor_res[paste0("F_score.", i_k)] <- 0
  }

  for (i_k in 1:n_K) {

    message(paste0("\tmarginal kernel ", i_k))

    ## pre-calculation for parameter estimation
    K_s <- K_s_list[[i_k]]

    sqK_s <- K_s %*% K_s
    K_sI_s <- K_s %*% I_s

    trsqK_s <- sum(diag(sqK_s))
    trK_sI_s <- sum(diag(K_sI_s))

    S <- matrix(c(trsqK_s, trK_sI_s, trK_sI_s, trsqI_s), nrow = 2, ncol = 2)


    ## parameter estimate for ligands
    for (r in 1:length(ligand_list)) {
      t_name <- ligand_list[r]
      t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))

      gen_y <- rawcount[t_name, , drop = FALSE]
      gen_y <- as.vector(apply(gen_y, 2, sum))

      # first order
      a_hat <- sum(gt_s * gen_y) / sum(gt_s^2)

      # second order
      mu_s <- gt_s * a_hat
      y_v <- gen_y - mu_s

      q <- rep(0, 2)
      q[1] <- t(y_v) %*% K_s %*% y_v - sum(diag(K_s) * mu_s)
      q[2] <- t(y_v) %*% I_s %*% y_v - sum(diag(I_s) * mu_s)
      q <- q / sl

      # solve theta
      mean_q <- mean(q)
      S_adj <- S / mean_q
      q_adj <- q / mean_q
      theta_hat <- solve(S_adj) %*% q_adj
      sqsigma_hat <- theta_hat[1,1]
      sqsigma_e_hat <- theta_hat[2,1]

      # correction
      if ((sqsigma_hat + sqsigma_e_hat) < 0) {
        ligand_res[r, paste0("valid.", i_k)] <- "False"
      } else if (sqsigma_hat < 0) {
        sqsigma_hat <- 0.
        sqsigma_e_hat <- q_adj[2]/S_adj[2,2]
        ligand_res[r, paste0("valid.", i_k)] <- "corrected"
      } else if (sqsigma_e_hat < 0) {
        sqsigma_e_hat <- 0
        sqsigma_hat <- q_adj[1]/S_adj[1,1]
        ligand_res[r, paste0("valid.", i_k)] <- "corrected"
      }

      # fitting score
      F_score <- (sqsigma_hat^2)*trsqK_s + (sqsigma_e_hat^2)*trsqI_s + 2*sqsigma_hat*sqsigma_e_hat*trK_sI_s - 2*sqsigma_hat*q[1] - 2*sqsigma_e_hat*q[2]

      # store result
      ligand_res$a_hat[r] <- a_hat
      ligand_res[r, paste0("sqsigma_hat.", i_k)] <- sqsigma_hat
      ligand_res[r, paste0("sqsigma_e_hat.", i_k)] <- sqsigma_e_hat
      ligand_res[r, paste0("F_score.", i_k)] <- F_score
    }


    ## parameter estimate for receptors
    for (r in 1:length(receptor_list)) {
      t_name <- receptor_list[r]
      t_name <- unlist(strsplit(t_name, "+", fixed = TRUE))

      gen_y <- rawcount[t_name, , drop = FALSE]
      gen_y <- as.vector(apply(gen_y, 2, sum))

      # first order
      a_hat <- sum(gt_s * gen_y) / sum(gt_s^2)

      # second order
      mu_s <- gt_s * a_hat
      y_v <- gen_y - mu_s

      q <- rep(0, 2)
      q[1] <- t(y_v) %*% K_s %*% y_v - sum(diag(K_s) * mu_s)
      q[2] <- t(y_v) %*% I_s %*% y_v - sum(diag(I_s) * mu_s)
      q <- q / sl

      # solve theta
      mean_q <- mean(q)
      S_adj <- S / mean_q
      q_adj <- q / mean_q
      theta_hat <- solve(S_adj) %*% q_adj
      sqsigma_hat <- theta_hat[1,1]
      sqsigma_e_hat <- theta_hat[2,1]

      # correction
      if ((sqsigma_hat + sqsigma_e_hat) < 0) {
        receptor_res[r, paste0("valid.", i_k)] <- "False"
      } else if (sqsigma_hat < 0) {
        sqsigma_hat <- 0.
        sqsigma_e_hat <- q_adj[2]/S_adj[2,2]
        receptor_res[r, paste0("valid.", i_k)] <- "corrected"
      } else if (sqsigma_e_hat < 0) {
        sqsigma_e_hat <- 0
        sqsigma_hat <- q_adj[1]/S_adj[1,1]
        receptor_res[r, paste0("valid.", i_k)] <- "corrected"
      }

      # fitting score
      F_score <- (sqsigma_hat^2)*trsqK_s + (sqsigma_e_hat^2)*trsqI_s + 2*sqsigma_hat*sqsigma_e_hat*trK_sI_s - 2*sqsigma_hat*q[1] - 2*sqsigma_e_hat*q[2]

      # store result
      receptor_res$a_hat[r] <- a_hat
      receptor_res[r, paste0("sqsigma_hat.", i_k)] <- sqsigma_hat
      receptor_res[r, paste0("sqsigma_e_hat.", i_k)] <- sqsigma_e_hat
      receptor_res[r, paste0("F_score.", i_k)] <- F_score
    }
  }
  message("Done!")

  ## summarize the single variant MoM result and organize the cellchatdb accordingly
  # select kernels
  ligand_res$select.kernel <- apply(ligand_res[paste("F_score.", 1:n_K, sep = "")], 1, which.min)
  ligand_res$sqsigma_hat <- 0
  ligand_res$sqsigma_e_hat <- 0
  for (r in 1:length(ligand_list)) {
    idx_kernel <- ligand_res$select.kernel[r]
    ligand_res$sqsigma_hat[r] <- ligand_res[r, paste0("sqsigma_hat.", idx_kernel)]
    ligand_res$sqsigma_e_hat[r] <- ligand_res[r, paste0("sqsigma_e_hat.", idx_kernel)]
  }

  receptor_res$select.kernel <- apply(receptor_res[paste("F_score.", 1:n_K, sep = "")], 1, which.min)
  receptor_res$sqsigma_hat <- 0
  receptor_res$sqsigma_e_hat <- 0
  for (r in 1:length(receptor_list)) {
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


  ### MoM parameter estimation for delta
  message("Estimating parameters of the joint spatial process and performing statistical inference...")
  ## pre-calculation for parameter estimation
  message("\tpre-calculation for parameter estimation")
  stopifnot(dim(K_s_list[[1]])[1] == n)
  Kx <- Wx / sl
  Kx_s <- t(t(Kx * gt_s) * gt_s) #note: Kx can be non-symmetric
  crKx_s <- crossprod(Kx_s, Kx_s) #equal to t(Kx_s) %*% Kx_s
  trcrKx_s <- sum(diag(crKx_s))

  ## scan each ligand-receptor pair in pairdb
  message("\tparameter estimation for each pair in pairdb")
  n_pair <- dim(pairdb)[1]

  pairdb$delta_hat <- 0

  for (r in 1:n_pair) {
    pairdb_tmp <- pairdb[r, ]

    l_name <- pairdb_tmp$ligand_organized_qc
    r_name <- pairdb_tmp$receptor_organized_qc
    l_name <- unlist(strsplit(l_name, "+", fixed = TRUE))
    r_name <- unlist(strsplit(r_name, "+", fixed = TRUE))

    gen_y1 <- rawcount[l_name, , drop = FALSE]
    gen_y1 <- as.vector(apply(gen_y1, 2, sum))
    gen_y2 <- rawcount[r_name, , drop = FALSE]
    gen_y2 <- as.vector(apply(gen_y2, 2, sum))

    a1_hat <- pairdb_tmp$a1_hat
    a2_hat <- pairdb_tmp$a2_hat

    mu1_s <- gt_s * a1_hat
    mu2_s <- gt_s * a2_hat
    y1_v <- gen_y1 - mu1_s
    y2_v <- gen_y2 - mu2_s

    # estimate delta
    delta_hat <- t(y1_v) %*% Kx_s %*% y2_v / sl / trcrKx_s
    pairdb$delta_hat[r] <- delta_hat
  }


  ### MoM statistical inference for delta
  ## pre-calculation for statistical inference
  message("\tpre-calculation for statistical inference")
  term_se_names <- c("term_se_1_1", "term_se_1_2",
                     "term_se_2_1", "term_se_2_2", "term_se_2_3", "term_se_2_4", "term_se_2_5",
                     "term_se_3_1", "term_se_3_2", "term_se_3_3", "term_se_3_4")
  stopifnot(length(K_s_list) == n_K)
  term_se_precalculate <- matrix(0, nrow = n_K^2, ncol = length(term_se_names))
  colnames(term_se_precalculate) <- term_se_names
  kernel_type_names <- paste(rep(1:n_K, rep(n_K, n_K)), rep(1:n_K, n_K), sep = "-")
  rownames(term_se_precalculate) <- kernel_type_names

  ## terms not related to marginal kernels
  message("\tterms not related to marginal kernels")
  Kx_s_dot_Kx_s <- Kx_s^2
  term_se_1_1 <- sum(t(Kx_s_dot_Kx_s * gt_s) * gt_s)
  term_se_1_2 <- sum(Kx_s_dot_Kx_s * Kx_s)
  term_se_precalculate[, "term_se_1_1"] <- term_se_1_1
  term_se_precalculate[, "term_se_1_2"] <- term_se_1_2

  term_se_2_1 <- sum(diag(crKx_s %*% crKx_s))
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

  ## terms related to only one marginal kernel
  message("\tterms related to only one marginal kernel")
  for (i_k in 1:n_K) {
    K_s <- K_s_list[[i_k]]

    Kx_sTK_s <- crossprod(Kx_s, K_s) #check: all.equal(crossprod(Kx_s, K_s), t(Kx_s) %*% K_s)
    Kx_sK_s <- Kx_s %*% K_s

    names_1 <- rep(1:n_K, rep(n_K, n_K))
    names_2 <- rep(1:n_K, n_K)

    term_se_2_3 <- sum(diag(Kx_sTK_s %*% Kx_sI_s))
    term_se_precalculate[which(names_1 == i_k), "term_se_2_3"] <- term_se_2_3

    term_se_2_4 <- sum(diag(Kx_sK_s %*% Kx_sTI_s))
    term_se_precalculate[which(names_2 == i_k), "term_se_2_4"] <- term_se_2_4

    term_se_3_1 <- sum(s_Kx_s_Kx_s_outer * K_s)
    term_se_precalculate[which(names_2 == i_k), "term_se_3_1"] <- term_se_3_1

    term_se_3_3 <- sum(s_Kx_sT_Kx_sT_outer * K_s)
    term_se_precalculate[which(names_1 == i_k), "term_se_3_3"] <- term_se_3_3
  }

  ## terms related to both two marginal kernels
  message("\tterms related to both two marginal kernels")
  for (i_k1 in 1:n_K) {
    # kernel 1
    K1_s <- K_s_list[[i_k1]]

    Kx_sTK1_s <- crossprod(Kx_s, K1_s) #check: all.equal(crossprod(Kx_s, K1_s), t(Kx_s) %*% K1_s)
    Kx_sTK1_s_Kx_s <- Kx_sTK1_s %*% Kx_s

    for (i_k2 in 1:n_K) {
      # kernel 2
      K2_s <- K_s_list[[i_k2]]

      term_se_2_2 <- sum(diag(Kx_sTK1_s_Kx_s %*% K2_s))
      term_se_precalculate[paste(i_k1, i_k2, sep = "-"), "term_se_2_2"] <- term_se_2_2
    }
  }

  term_se_precalculate[, "term_se_1_1"] <- term_se_precalculate[, "term_se_1_1"] / (sl^2)
  term_se_precalculate[, "term_se_1_2"] <- term_se_precalculate[, "term_se_1_2"] / sl
  term_se_precalculate[, "term_se_3_1"] <- term_se_precalculate[, "term_se_3_1"] / sl
  term_se_precalculate[, "term_se_3_2"] <- term_se_precalculate[, "term_se_3_2"] / sl
  term_se_precalculate[, "term_se_3_3"] <- term_se_precalculate[, "term_se_3_3"] / sl
  term_se_precalculate[, "term_se_3_4"] <- term_se_precalculate[, "term_se_3_4"] / sl

  ## scan each ligand-receptor pair in pairdb
  message("\tstatistical inference for each pair in pairdb")
  stopifnot(dim(pairdb)[1] == n_pair)

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

    # calculate se
    coef_se_vec <- c(a1_hat*a2_hat, delta_hat,
                     delta_hat^2, sqsigma1_hat*sqsigma2_hat, sqsigma1_hat*sqsigma_e2_hat, sqsigma2_hat*sqsigma_e1_hat, sqsigma_e1_hat*sqsigma_e2_hat,
                     a1_hat*sqsigma2_hat, a1_hat*sqsigma_e2_hat, a2_hat*sqsigma1_hat, a2_hat*sqsigma_e1_hat)
    term_se_vec <- as.vector(term_se_precalculate[kernel_type, ])
    se_delta_hat <- sqrt(sum(coef_se_vec * term_se_vec)) / trcrKx_s

    # store se
    pairdb$se_delta_hat[r] <- se_delta_hat
  }

  ## calculate p-values
  pairdb$z_score <- pairdb$delta_hat / pairdb$se_delta_hat
  W <- pairdb$z_score ^ 2
  P_value <- pchisq(W, 1, lower.tail=F)
  pairdb$p_value <- P_value
  message("Done!")

  return(pairdb)
}

