generateCovarianceMatrixSynth = function(J, s_loading, g_loading){
  t.b1 = c(rep(sqrt(s_loading), J/2), rep(0, J/2))
  t.b2 = rev(t.b1)
  t.b3 = rep(sqrt(g_loading), J)
  lamb = cbind(t.b1, t.b2, t.b3)
  psi  = diag(3)
  t.rho = lamb %*% psi %*% t(lamb)
  diag(t.rho) = 1
  t.rho
}

simulateDataSynth = function(n, J, sigma, var_names = paste0("y", 1:J)){
  Y = mvtnorm::rmvnorm(n = n, mean = rep(0, J), sigma = sigma)
  colnames(Y) = var_names
  Y
}


fitSemModelSynth = function(data){
  model =   '
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  
  f1 ~~ f2
  '
  sem(model, 
      data = as.data.frame(data),
      mimic   = "lavaan",
      missing = "ML",
      std.lv  = TRUE,
      std.ov  = FALSE,
      se      = "standard",
      bootstrap = 1000)
}

fitSemModelSynthCov = function(Sigma, N){
  colnames(Sigma) = rownames(Sigma) = paste0("y", 1:ncol(Sigma))
  model =   '
  f1 =~ y1 + y2 + y3
  f2 =~ y4 + y5 + y6
  
  f1 ~~ f2
  '
  sem(model,
      sample.cov = Sigma,
      sample.nobs = N,
      mimic   = "lavaan",
      missing = "ML",
      std.lv  = TRUE,
      std.ov  = FALSE,
      se      = "standard",
      bootstrap = 1000)
}




# Extract parameter estimates
extractEstimatesSynth = function(fit){
  values_first_layer = subset(standardizedSolution(fit), op == "=~")
  vals1 = rev(values_first_layer$est.std)
  values_first_layer_ind = .map2idx(vals1)

  values_second_layer = subset(standardizedSolution(fit), 
                                op == "~~" & lhs != rhs & 
                                lhs %in% c("f1", "f2"))
  vals2 = rev(values_second_layer$est.std)
  values_second_layer_ind = .map2idx(vals2)

  list(
    vals1 = vals1,
    vals2 = vals2,
    vals1_ind = values_first_layer_ind,
    vals2_ind = values_second_layer_ind
  )
}

plotSemDiagramSynth = function(fit, shades, cols, J = 6){
  estimates = extractEstimatesSynth(fit)
  
  factors = rev(expression(F[1], F[2]))
  select_cols = paste("Task", 1:J)
  
  plot(NA, NA, xlim = c(0, 1.2), ylim = c(0, 2), axes = FALSE, 
       xlab = "", ylab = "")

  task_pos.y = seq(0.2, 1.8, length.out = J + 1)
  task_pos.y = task_pos.y[-4]
  task_pos.x = rep(0.2, J)
  task.length = c(0.1, 0.1)
  cex_val = 0.7
  
  # Draw observed variable boxes
  rect(task_pos.x + task.length[1], task_pos.y + task.length[2],
       task_pos.x - task.length[1], task_pos.y - task.length[2], lwd = 2)
  text(task_pos.x, task_pos.y, rev(select_cols), cex = .8)
  
  # Draw latent factors
  n_factors = 2
  factor_pos.y = c(task_pos.y[2], task_pos.y[5])
  factor_pos.x = rep(0.7, n_factors)
  radius_x = 0.1
  
  for (i in 1:n_factors) {
    .drawEllipse(factor_pos.x[i], factor_pos.y[i], 
                 radius_y = 0.1, radius_x = radius_x)
  }
  text(factor_pos.x, factor_pos.y, factors, cex = 1.15, col = rev(unique(cols)))
  
  # Draw factor loadings
  first_con = cbind(task_pos.x + task.length[1], task_pos.y)
  second_con = cbind(rep(factor_pos.x[1] - radius_x, J),
                     c(rep(factor_pos.y, each = 3)))
  
  cols_line_first = sapply(1:J, function(x) shades[estimates$vals1_ind[x]])
  
  for (i in 1:J) {
    lines(c(first_con[i, 1], second_con[i, 1]),
          c(first_con[i, 2], second_con[i, 2]),
          lwd = 4, col = cols_line_first[i])
    
    mid_x = (first_con[i, 1] + second_con[i, 1]) / 2 + 0.05
    mid_y = (first_con[i, 2] + second_con[i, 2]) / 2 
    rect(mid_x - 0.019, mid_y + 0.05,
         mid_x - 0.082, mid_y - 0.05, col = "white", border = "white")
    text(mid_x - 0.05, mid_y, label = round(estimates$vals1[i], 3),
         cex = .8*cex_val, col = "black")
  }
  
  # Draw factor correlations
  cols_line_second = sapply(1:length(estimates$vals2), function(x) shades[estimates$vals2_ind[x]])
  dists = 0.5
  count = 1
  
  for (i in 1:n_factors) {
    for (j in 1:n_factors) {
      if (i > j) {
        .curveFunction(
          c(factor_pos.x[i] + radius_x, factor_pos.y[i]),
          c(factor_pos.x[j] + radius_x, factor_pos.y[j]),
          dists,
          col = cols_line_second[count]
        )
        
        mid_x = factor_pos.x[i] + dists/2 + 0.1
        mid_y = (factor_pos.y[i] + factor_pos.y[j]) / 2 - 0.08
        rect(mid_x - 0.15, mid_y + 0.15,
             mid_x + 0.15, mid_y - 0, col = "white", border = "white")
        text(mid_x, mid_y + 0.08,
             labels = sprintf("%.3f", estimates$vals2[count]),
             cex = 1*cex_val, col = "black")
        count = count + 1
      }
    }
  }
}

plotCorrelationMatrixSynth = function(R_hat, shades, cols){
  h = 7
  corrplot(unname(R_hat), 
           method = "color", 
           cl.pos = "n", 
           col = shades,
           mar = c(0, 4, 2, 0),
           tl.pos = "n",
           cl.cex = 2,
           addCoef.col = "black")
  
  for (i in 1:ncol(R_hat)) {
    mtext(i, side = 2, line = .5, at = h - i,
          cex = 1, adj = 0.5, las = 1, col = cols[i])
    mtext(7 - i, side = 3, line = .5, at = h - i,
          cex = 1, adj = 0.5, las = 1, col = cols[7 - i])
  }
}


makeCorInt1 = function(lambda_A, lambda_B, rho_star, J = 6) {
    if (length(lambda_A) == 1) lambda_A = rep(lambda_A, J/2)
    if (length(lambda_B) == 1) lambda_B = rep(lambda_B, J/2)
    Lambda = matrix(c(lambda_A, rep(0, J), lambda_B), nrow = J, ncol = 2)
    Phi = matrix(c(1, rho_star, rho_star, 1), nrow = 2, ncol = 2)
    Sigma = Lambda %*% Phi %*% t(Lambda)
    diag(Sigma) = 1
    return(Sigma)
}

makeCorInt2 = function(lambda_A, lambda_B, lambda_G, J = 6) {
    if (length(lambda_A) == 1) lambda_A = rep(lambda_A, J/2)
    if (length(lambda_B) == 1) lambda_B = rep(lambda_B, J/2)
    if (length(lambda_G) == 1) lambda_G = rep(lambda_G, J)
    Lambda = matrix(c(lambda_A, rep(0, J), lambda_B, lambda_G), nrow = J, ncol = 3)
    Sigma = tcrossprod(Lambda)
    diag(Sigma) = 1
    return(Sigma)
}


getARs = function(Sigma, indA, indB){
    A = Sigma[indA, indA]
    B = Sigma[indB, indB]
    AB = Sigma[indA, indB]

    eps = 1e-10

    pairsA = combn(seq_along(indA), 2)  
    pairsB = combn(seq_along(indB), 2)  

    tetrads = c()
    for (pa in 1:ncol(pairsA)) {
        i1 = pairsA[1, pa]
        i2 = pairsA[2, pa]

        for (pb in 1:ncol(pairsB)) {
            j1 = pairsB[1, pb]
            j2 = pairsB[2, pb]

            denom = A[i1, i2] * B[j1, j2]
            num = AB[i1, j2] * AB[i2, j1]

            if (is.na(num) || is.na(denom)) next
            if (abs(denom) < eps) next

            tetrads = c(tetrads, num / denom)
        }
    }
    return(tetrads)
}



generateSigmaBifactorTemp = function(q0_mean, q0_sd, resid_var, common_var = 1, J, indA, indB){
  
  makeQ0 = function(q0_mean, q0_sd, m){
    w = rlnorm(m, meanlog = 0, sdlog = q0_sd)
    
    kEst = function(k){
      q0 = (k * w) / (1 + k * w)   
      mean(q0) - q0_mean
    }
    
    lo = 1e-12
    hi = 1
    while (kEst(hi) < 0) hi = hi * 2
    k_hat = uniroot(kEst, c(lo, hi))$root
    (k_hat * w) / (1 + k_hat * w)
  }
  
  mA = length(indA)
  mB = length(indB)
  
  q0A = makeQ0(q0_mean, q0_sd, mA)
  q0B = makeQ0(q0_mean, q0_sd, mB)
  q0  = c(q0A, q0B)
  
  r = sqrt((1 - q0) / q0)        
  lambda_g = sqrt(common_var / (1 + r^2))
  lambda_s = r * lambda_g
  
  Lambda = matrix(0, nrow = J, ncol = 3)
  colnames(Lambda) = c("g", "sA", "sB")
  Lambda[, "g"] = lambda_g
  Lambda[indA, "sA"] = lambda_s[indA]
  Lambda[indB, "sB"] = lambda_s[indB]
  
  Sig = tcrossprod(Lambda) + diag(resid_var, J)
  out = list(
    Sigma = Sig,
    Lambda = Lambda
  )
  out
}