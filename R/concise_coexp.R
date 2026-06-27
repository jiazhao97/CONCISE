#' CONCISE spatial gene-gene co-expression analysis
#'
#' @export
CONCISE_coexp <- function(
    gene_list,
    rawcount,
    loc,
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
  message("Done!")

  ## calculate libarary size
  if (is.null(lib_size)) {
    lib_size <- as.vector(apply(rawcount, 2, sum))
  }
  n <- length(lib_size)

  ## get ground truth library size to make the scale of normalized gene expressions roughly having sd = 1
  genes_use <- gene_list
  rawcount_adj_lbs <- t(t(rawcount) / lib_size)
  sd_exp_level <- apply(rawcount_adj_lbs[genes_use, ], 1, sd)
  aver_sd <- mean(sd_exp_level)
  gt_s <- lib_size * aver_sd


  ### iterate over choices of kernels and finish all the estimate based on the current kernel
  message("Iterating over choices of kernels and estimating parameters of marginal distributions...")
  stopifnot(length(W_list) == n_K)
  stopifnot(dim(W_list[[1]])[1] == n)
  sl <- sqrt(n)

  I_s <- diag((gt_s^2) / sl)
  I_s <- Matrix(I_s, sparse = TRUE)
  trsqI_s <- sum(I_s * I_s)

  K_s_list <- list()
  for (i_k in 1:n_K) {
    K <- W_list[[i_k]] / sl
    K_s <- K * tcrossprod(gt_s)
    K_s <- 0.5 * (K_s + t(K_s))
    K_s_list[[i_k]] <- K_s
  }
  rm(W_list)

  hvg_list <- gene_list
  hvg_res <- data.frame(hvg = hvg_list)
  hvg_res$a_hat <- 0
  for (i_k in 1:n_K) {
    hvg_res[paste0("sqsigma_hat.", i_k)] <- 0
    hvg_res[paste0("sqsigma_e_hat.", i_k)] <- 0
    hvg_res[paste0("valid.", i_k)] <- "True"
    hvg_res[paste0("F_score.", i_k)] <- 0
  }

  for (i_k in 1:n_K) {

    message(paste0("\tmarginal kernel ", i_k))

    ## pre-calculation for parameter estimation
    K_s <- K_s_list[[i_k]]

    trsqK_s <- sum(K_s * K_s)
    trK_sI_s <- sum(K_s * I_s)

    S <- matrix(c(trsqK_s, trK_sI_s, trK_sI_s, trsqI_s), nrow = 2, ncol = 2)

    ## parameter estimate for hvg
    gen_y <- t(rawcount[hvg_list, ]) #cell by gene

    # first order
    a_hat <- as.vector(colSums(gen_y * gt_s)) / sum(gt_s^2)
    hvg_res$a_hat <- a_hat

    # second order
    mu_s <- gt_s %*% t(a_hat) #cell by gene
    y_v <- gen_y - mu_s

    q <- matrix(0, nrow = 2, ncol = length(hvg_list))
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
    hvg_res[paste0("sqsigma_hat.", i_k)] <- sqsigma_hat
    hvg_res[paste0("sqsigma_e_hat.", i_k)] <- sqsigma_e_hat

    # correction
    hvg_res[(sqsigma_hat + sqsigma_e_hat) < 0, paste0("valid.", i_k)] <- "False"
    idx_valid_false <- ((sqsigma_hat + sqsigma_e_hat) < 0)

    idx_sqsigma_hat <- (sqsigma_hat < 0)
    hvg_res[idx_sqsigma_hat, paste0("sqsigma_hat.", i_k)] <- 0.
    hvg_res[idx_sqsigma_hat, paste0("sqsigma_e_hat.", i_k)] <- q_adj[2, idx_sqsigma_hat] / S_adj[2,2]
    hvg_res[idx_sqsigma_hat, paste0("valid.", i_k)] <- "corrected"

    idx_sqsigma_e_hat <- (sqsigma_e_hat < 0)
    hvg_res[idx_sqsigma_e_hat, paste0("sqsigma_e_hat.", i_k)] <- 0.
    hvg_res[idx_sqsigma_e_hat, paste0("sqsigma_hat.", i_k)] <- q_adj[1, idx_sqsigma_e_hat] / S_adj[1,1]
    hvg_res[idx_sqsigma_e_hat, paste0("valid.", i_k)] <- "corrected"

    hvg_res[idx_valid_false, paste0("valid.", i_k)] <- "False"

    # fitting score
    F_score <- (sqsigma_hat^2)*trsqK_s + (sqsigma_e_hat^2)*trsqI_s + 2*sqsigma_hat*sqsigma_e_hat*trK_sI_s - 2*sqsigma_hat*q[1, ] - 2*sqsigma_e_hat*q[2, ]
    hvg_res[paste0("F_score.", i_k)] <- F_score
  }
  message("Done!")

  ## summarize the single variant MoM result and organize the corrdb accordingly
  # select kernels
  if (remove_invalid_genes) {
    for (i_k in 1:n_K) {
      hvg_res <- hvg_res[hvg_res[paste0("valid.", i_k)] != "False", ]
    }
  }
  hvg_list <- as.vector(hvg_res$hvg)
  hvg_res$select.kernel <- apply(hvg_res[paste("F_score.", 1:n_K, sep = "")], 1, which.min)
  hvg_res$sqsigma_hat <- 0
  hvg_res$sqsigma_e_hat <- 0
  for (r in 1:length(hvg_list)) {
    idx_kernel <- hvg_res$select.kernel[r]
    hvg_res$sqsigma_hat[r] <- hvg_res[r, paste0("sqsigma_hat.", idx_kernel)]
    hvg_res$sqsigma_e_hat[r] <- hvg_res[r, paste0("sqsigma_e_hat.", idx_kernel)]
  }
  rownames(hvg_res) <- as.vector(hvg_res$hvg)

  corrdb <- data.frame(gene1 = rep(hvg_list, rep(length(hvg_list), length(hvg_list))),
                       gene2 = rep(hvg_list, length(hvg_list)))
  stopifnot(dim(corrdb)[1] == (dim(hvg_list)[1])^2)

  hvg_res_tmp <- hvg_res[corrdb$gene1, ]
  stopifnot(sum(hvg_res_tmp$hvg != corrdb$gene1) == 0)
  corrdb$a1_hat <- hvg_res_tmp$a_hat
  corrdb$sqsigma1_hat <- hvg_res_tmp$sqsigma_hat
  corrdb$sqsigma_e1_hat <- hvg_res_tmp$sqsigma_e_hat
  corrdb$select_kernel1 <- hvg_res_tmp$select.kernel

  hvg_res_tmp <- hvg_res[corrdb$gene2, ]
  stopifnot(sum(hvg_res_tmp$hvg != corrdb$gene2) == 0)
  corrdb$a2_hat <- hvg_res_tmp$a_hat
  corrdb$sqsigma2_hat <- hvg_res_tmp$sqsigma_hat
  corrdb$sqsigma_e2_hat <- hvg_res_tmp$sqsigma_e_hat
  corrdb$select_kernel2 <- hvg_res_tmp$select.kernel

  corrdb$kernel_type <- paste(corrdb$select_kernel1, corrdb$select_kernel2, sep = "-")


  ### MoM parameter estimation for delta
  message("Estimating parameters of the joint spatial process and performing statistical inference...")
  ## pre-calculation for parameter estimation
  stopifnot(dim(K_s_list[[1]])[1] == n)
  Kx <- diag(n) / sl
  Kx_s <- Kx * tcrossprod(gt_s) # all.equal(Kx * tcrossprod(gt_s), t(t(Kx * gt_s) * gt_s))
  Kx_s <- Matrix(Kx_s, sparse = TRUE)
  crKx_s <- crossprod(Kx_s, Kx_s) #equal to t(Kx_s) %*% Kx_s
  trcrKx_s <- sum(diag(crKx_s))

  ## scan each ligand-receptor pair in pairdb
  message("\tparameter estimation for each pair in pairdb")
  gen_y <- t(rawcount[hvg_list, ]) # cell by gene
  stopifnot(sum(corrdb$gene2[1:length(hvg_list)] != hvg_list) == 0)
  base_y <- gt_s %*% t(hvg_res$a_hat) # cell by gene
  stopifnot(dim(gen_y)[1] == dim(base_y)[1])
  stopifnot(dim(gen_y)[2] == dim(base_y)[2])
  y_v <- gen_y - base_y # cell by gene
  delta_hat_mat <- crossprod(y_v, Kx_s %*% y_v) / sl / trcrKx_s # all.equal(crossprod(y_v, Kx_s %*% y_v), t(y_v) %*% Kx_s %*% y_v)
  stopifnot(dim(corrdb)[1] == (dim(delta_hat_mat)[1] * dim(delta_hat_mat)[2]))
  corrdb$delta_hat <- as.vector(delta_hat_mat)


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
  term_se_1_1 <- sum(Kx_s_dot_Kx_s * tcrossprod(gt_s))
  term_se_1_2 <- sum(Kx_s_dot_Kx_s * Kx_s)
  term_se_precalculate[, "term_se_1_1"] <- term_se_1_1
  term_se_precalculate[, "term_se_1_2"] <- term_se_1_2

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

  ## terms related to only one marginal kernel
  message("\tterms related to only one marginal kernel")
  Kx_sTK_s_list <- list()
  Kx_sK_s_list <- list()
  for (i_k in 1:n_K) {
    K_s <- K_s_list[[i_k]]
    Kx_sTK_s_list[[i_k]] <- crossprod(Kx_s, K_s) #check: all.equal(crossprod(Kx_s, K_s), t(Kx_s) %*% K_s)
    Kx_sK_s_list[[i_k]] <- Kx_s %*% K_s
  }

  s_Kx_s_Kx_s_outer <- as.matrix(s_Kx_s_Kx_s_outer)
  s_Kx_sT_Kx_sT_outer <- as.matrix(s_Kx_sT_Kx_sT_outer)
  for (i_k in 1:n_K) {
    Kx_sTK_s <- Kx_sTK_s_list[[i_k]]
    Kx_sK_s <- Kx_sK_s_list[[i_k]]

    names_1 <- rep(1:n_K, rep(n_K, n_K))
    names_2 <- rep(1:n_K, n_K)

    term_se_2_3 <- sum(t(Kx_sI_s) * Kx_sTK_s)
    term_se_precalculate[which(names_1 == i_k), "term_se_2_3"] <- term_se_2_3

    term_se_2_4 <- sum(t(Kx_sTI_s) * Kx_sK_s)
    term_se_precalculate[which(names_2 == i_k), "term_se_2_4"] <- term_se_2_4

    term_se_3_1 <- sum(s_Kx_s_Kx_s_outer * K_s)
    term_se_precalculate[which(names_2 == i_k), "term_se_3_1"] <- term_se_3_1

    term_se_3_3 <- sum(s_Kx_sT_Kx_sT_outer * K_s)
    term_se_precalculate[which(names_1 == i_k), "term_se_3_3"] <- term_se_3_3
  }
  rm(s_Kx_s_Kx_s_outer)
  rm(s_Kx_sT_Kx_sT_outer)

  ## terms related to both two marginal kernels
  message("\tterms related to both two marginal kernels")
  for (i_k1 in 1:n_K) {
    # kernel 1
    Kx_sTK1_s <- Kx_sTK_s_list[[i_k1]]

    for (i_k2 in 1:n_K) {
      # kernel 2
      Kx_sK2_s <- Kx_sK_s_list[[i_k2]]

      # term_se_2_2 <- sum(diag(Kx_sTK1_s %*% Kx_sK2_s))
      term_se_2_2 <- sum(t(Kx_sTK1_s) * Kx_sK2_s)
      term_se_precalculate[paste(i_k1, i_k2, sep = "-"), "term_se_2_2"] <- term_se_2_2
    }
  }

  term_se_precalculate[, "term_se_1_1"] <- term_se_precalculate[, "term_se_1_1"] / (sl^2)
  term_se_precalculate[, "term_se_1_2"] <- term_se_precalculate[, "term_se_1_2"] / sl
  term_se_precalculate[, "term_se_3_1"] <- term_se_precalculate[, "term_se_3_1"] / sl
  term_se_precalculate[, "term_se_3_2"] <- term_se_precalculate[, "term_se_3_2"] / sl
  term_se_precalculate[, "term_se_3_3"] <- term_se_precalculate[, "term_se_3_3"] / sl
  term_se_precalculate[, "term_se_3_4"] <- term_se_precalculate[, "term_se_3_4"] / sl

  kernel_type_list <- rownames(term_se_precalculate)
  rownames(corrdb) <- paste(corrdb$gene1, corrdb$gene2, sep="-")
  corrdb$se_delta_hat <- 0
  for (kt in kernel_type_list) {

    pair_names <- rownames(corrdb)[corrdb$kernel_type == kt]
    n_pair_names <- length(pair_names)

    corrdb_pairs <- corrdb[pair_names, ]
    stopifnot(dim(corrdb_pairs)[1] == n_pair_names)
    coef_se_vec <- matrix(0, nrow = n_pair_names, ncol = 11)

    coef_se_vec[, 1] <- corrdb_pairs$a1_hat*corrdb_pairs$a2_hat
    coef_se_vec[, 2] <- corrdb_pairs$delta_hat

    coef_se_vec[, 3] <- corrdb_pairs$delta_hat^2
    coef_se_vec[, 4] <- corrdb_pairs$sqsigma1_hat*corrdb_pairs$sqsigma2_hat
    coef_se_vec[, 5] <- corrdb_pairs$sqsigma1_hat*corrdb_pairs$sqsigma_e2_hat
    coef_se_vec[, 6] <- corrdb_pairs$sqsigma2_hat*corrdb_pairs$sqsigma_e1_hat
    coef_se_vec[, 7] <- corrdb_pairs$sqsigma_e1_hat*corrdb_pairs$sqsigma_e2_hat

    coef_se_vec[, 8] <- corrdb_pairs$a1_hat*corrdb_pairs$sqsigma2_hat
    coef_se_vec[, 9] <- corrdb_pairs$a1_hat*corrdb_pairs$sqsigma_e2_hat
    coef_se_vec[, 10] <- corrdb_pairs$a2_hat*corrdb_pairs$sqsigma1_hat
    coef_se_vec[, 11] <- corrdb_pairs$a2_hat*corrdb_pairs$sqsigma_e1_hat

    term_se_vec <- as.vector(term_se_precalculate[kt, ])

    se_delta_hat <- t(t(coef_se_vec) * term_se_vec)
    se_delta_hat <- rowSums(se_delta_hat)
    se_delta_hat <- sqrt(se_delta_hat) / trcrKx_s
    stopifnot(length(se_delta_hat) == n_pair_names)

    corrdb[pair_names, ]$se_delta_hat <- se_delta_hat
  }

  ## calculate p-values
  corrdb$z_score <- corrdb$delta_hat / corrdb$se_delta_hat
  corrdb$corr_hat <- corrdb$delta_hat / sqrt(corrdb$sqsigma1_hat+corrdb$sqsigma_e1_hat) / sqrt(corrdb$sqsigma2_hat+corrdb$sqsigma_e2_hat)
  corrdb$corr_hat[corrdb$corr_hat > 1] <- 1
  corrdb$corr_hat[corrdb$corr_hat < -1] <- -1
  W <- corrdb$z_score ^ 2
  P_value <- pchisq(W, 1, lower.tail=F)
  corrdb$p_value <- P_value

  return(corrdb)
}

