library(lavaan)

makeSigma = function(B) {
  Sigma = tcrossprod(B)
  diag(Sigma) = 1
  colnames(Sigma) = rownames(Sigma) = paste0("V", 1:J)
  Sigma
}

getStdLoadings = function(fit, A_idx, B_idx) {
  ss = standardizedSolution(fit)
  ss = ss[ss$op == "=~", c("rhs", "est.std")]
  lam = ss$est.std
  names(lam) = ss$rhs
  list(
    lamA = lam[paste0("V", A_idx)],
    lamB = lam[paste0("V", B_idx)]
  )
}

contribStats = function(lamb_A, lamb_B, R_AB) {
  DA = drop(crossprod(lamb_A))
  DB = drop(crossprod(lamb_B))
  Wlin = (lamb_A %o% lamb_B) / (DA * DB)
  C = Wlin * R_AB
  tot = sum(abs(C))
  list(
    cum = cumsum(sort(as.vector(abs(C)), decreasing = TRUE)) / tot
  )
}

plotContributionAuc = function(lamb_A, lamb_B, R_AB, rho_star, blend_ratio) {
  out = contribStats(lamb_A, lamb_B, R_AB)

  y = c(0, out$cum)
  n = length(y)
  x = seq(0, 1, length.out = n)

  par(pty = "s")
  plot(x, y,
       type = "l", lwd = 2,
       xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1),
       axes = FALSE, asp = 1)

  poly_col = grDevices::adjustcolor("gray70", alpha.f = 0.4)
  polygon(c(x, rev(x)), c(rep(0, n), rev(y)), border = NA, col = poly_col)

  box()
  abline(0, 1, lty = 2, col = "indianred3", lwd = 3)

  axis(1, at = c(0, 1), labels = c(0, 1))
  axis(1, at = x, labels = FALSE)
  axis(2, at = c(0, 1), labels = c(0, 1), las = 1)
  axis(2, at = x, labels = FALSE)

  points(x, y, pch = 16, cex = 0.6)

  # txt = bquote(rho^"*" == .(format(round(rho_star, 3), nsmall = 3)))
  # text(0.65, 0.20, txt, cex = 1.2)
  
}


mod = "
  f1 =~ V1 + V2 + V3 + V4 + V5
  f2 =~ V6 + V7 + V8 + V9 + V10
  f1 ~~ f2
"

plotSemDiagramSynth2 = function(fit, shades, cols, J = 10){
  estimates = extractEstimatesSynth(fit)
  
  factors = rev(expression(F[1], F[2]))
  select_cols = paste("Y", 1:J)
  
  plot(NA, NA, xlim = c(0, 1.2), ylim = c(0, 2), axes = FALSE, 
       xlab = "", ylab = "")

  task_pos.y = seq(0.2, 1.8, length.out = J + 1)
  task_pos.y = task_pos.y[-6]
  task_pos.x = rep(0.2, J)
  task.length = c(0.05, 0.05)
  cex_val = 0.7
  
  rect(task_pos.x + task.length[1], task_pos.y + task.length[2],
       task_pos.x - task.length[1], task_pos.y - task.length[2], lwd = 2)
  text(task_pos.x, task_pos.y, rev(select_cols), cex = .8)
  
  n_factors = 2
  factor_pos.y = c(task_pos.y[3], task_pos.y[8])
  factor_pos.x = rep(0.7, n_factors)
  radius_x = 0.1
  
  for (i in 1:n_factors) {
    .drawEllipse(factor_pos.x[i], factor_pos.y[i], 
                 radius_y = 0.1, radius_x = radius_x)
  }
  text(factor_pos.x, factor_pos.y, factors, cex = 1.15, col = rev(unique(cols)))
  
  first_con = cbind(task_pos.x + task.length[1], task_pos.y)
  second_con = cbind(
    rep(factor_pos.x[1] - radius_x, J),
    rep(factor_pos.y, each = J/n_factors)
  )
  
  cols_line_first = sapply(1:J, function(x) shades[estimates$vals1_ind[x]])
  
  for (i in 1:J) {
    lines(c(first_con[i, 1], second_con[i, 1]),
          c(first_con[i, 2], second_con[i, 2]),
          lwd = 4, col = cols_line_first[i])
    
    # mid_x = (first_con[i, 1] + second_con[i, 1]) / 2 + 0.05
    # mid_y = (first_con[i, 2] + second_con[i, 2]) / 2 
    # rect(mid_x - 0.019, mid_y + 0.05,
    #      mid_x - 0.082, mid_y - 0.05, col = "white", border = "white")
    # text(mid_x - 0.05, mid_y, labels = round(estimates$vals1[i], 3),
    #      cex = .8*cex_val, col = "black")
  }
  
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
        # text(mid_x, mid_y + 0.08,
        #      labels = sprintf("%.3f", estimates$vals2[count]),
        #      cex = 1*cex_val, col = "black")
        rho_star = estimates$vals2[count]
        txt = bquote(rho^"*" == .(format(round(rho_star, 3), nsmall = 3)))
        text(mid_x, mid_y + 0.08, txt, cex = 1.2)
        count = count + 1
      }
    }
  }
}


plotBifactorDiagram = function(B, shades, cols, J = 10) {
  general_loads = B[, 1]
  specific1_loads = B[1:(J/2), 2]
  specific2_loads = B[(J/2 + 1):J, 3]
  
  factors = c(expression(G), expression(S[1]), expression(S[2]))
  select_cols = paste("Y", J:1)
  
  plot(NA, NA, xlim = c(0, 1.4), ylim = c(0, 2), axes = FALSE, 
       xlab = "", ylab = "")
  
  task_pos.y = seq(0.2, 1.8, length.out = J + 1)
  task_pos.y = task_pos.y[-6]
  task_pos.x = rep(0.6, J)
  task.length = c(0.05, 0.05)
  
  rect(task_pos.x + task.length[1], task_pos.y + task.length[2],
       task_pos.x - task.length[1], task_pos.y - task.length[2], lwd = 2)
  text(task_pos.x, task_pos.y, rev(select_cols), cex = .8)
  
  general_pos.x = 0.2
  general_pos.y = 1.0
  radius_x = 0.1
  
  .drawEllipse(general_pos.x, general_pos.y, 
               radius_y = 0.1, radius_x = radius_x)
  text(general_pos.x, general_pos.y, factors[1], cex = 1.15, col = cols[1])
  
  specific_pos.x = rep(1.0, 2)
  specific_pos.y = c(task_pos.y[3], task_pos.y[8])
  
  for (i in 1:2) {
    .drawEllipse(specific_pos.x[i], specific_pos.y[i], 
                 radius_y = 0.1, radius_x = radius_x)
  }
  text(specific_pos.x, specific_pos.y, factors[2:3], 
       cex = 1.15, col = cols[2:3])
  
  general_con = cbind(
    rep(general_pos.x + radius_x, J),
    rep(general_pos.y, J)
  )
  task_con_left = cbind(task_pos.x - task.length[1], task_pos.y)
  
  vals1_ind = .map2idx(abs(general_loads))
  general_cols = sapply(1:J, function(x) shades[vals1_ind])

  
  for (i in 1:J) {
    lines(c(general_con[i, 1], task_con_left[i, 1]),
          c(general_con[i, 2], task_con_left[i, 2]),
          lwd = 4, col = general_cols[i])
  }
  
  task_con_right = cbind(task_pos.x + task.length[1], task_pos.y)
  
  specific1_con = cbind(
    rep(specific_pos.x[1] - radius_x, J/2),
    rep(specific_pos.y[1], J/2)
  )


  vals2_ind = .map2idx(abs(specific1_loads))
  specific1_cols = sapply(1:(J/2), function(x) shades[vals2_ind[x]])
  
  
  for (i in 1:(J/2)) {
    lines(c(task_con_right[i, 1], specific1_con[i, 1]),
          c(task_con_right[i, 2], specific1_con[i, 2]),
          lwd = 4, col = specific1_cols[i])
  }
  
  specific2_con = cbind(
    rep(specific_pos.x[2] - radius_x, J/2),
    rep(specific_pos.y[2], J/2)
  )
  
  vals3_ind = .map2idx(abs(specific2_loads))
  specific2_cols = sapply(1:(J/2), function(x) shades[vals3_ind[x]])
  
  for (i in 1:(J/2)) {
    j = i + J/2
    lines(c(task_con_right[j, 1], specific2_con[i, 1]),
          c(task_con_right[j, 2], specific2_con[i, 2]),
          lwd = 4, col = specific2_cols[i])
  }
}
getARs = function(Sigma, indA, indB){
  A  = Sigma[indA, indA, drop = FALSE]
  B  = Sigma[indB, indB, drop = FALSE]
  AB = Sigma[indA, indB, drop = FALSE]

  eps = 1e-10

  pairsA = combn(seq_along(indA), 2)  # unordered
  pairsB = combn(seq_along(indB), 2)  # unordered

  ars = c()
  for (pa in 1:ncol(pairsA)){
    i1 = pairsA[1, pa]
    i2 = pairsA[2, pa]

    for (pb in 1:ncol(pairsB)){
      j1 = pairsB[1, pb]
      j2 = pairsB[2, pb]

      denom = A[i1, i2] * B[j1, j2]
      if (is.na(denom) || abs(denom) < eps) next

      num1 = AB[i1, j1] * AB[i2, j2]
      num2 = AB[i1, j2] * AB[i2, j1]

      if (!is.na(num1)) ars = c(ars, num1 / denom)
      if (!is.na(num2)) ars = c(ars, num2 / denom)
    }
  }
  ars
}
plotTetradRatio = function(Sigma, indA, indB, rho_star){
  tetrads = getARs(Sigma, indA, indB)
  target_rho = prod(tetrads)^(1/length(tetrads))
  plot(sort(tetrads), ylim = c(0, 1), ylab = "", xlab = "", pch = 16, cex = 1.2, axes = FALSE)
  box()
  mtext(expression(Q[1]), side = 2, line = 2.5)
  axis(2, las = 1)
  axis(1, at = c(1,length(tetrads)), labels = c(1,length(tetrads)))   
  mtext("Ordered Alignment Ratios", side = 1, line = .5, cex = .7)
  abline(h = rho_star^2, col = "red3", lty = 2, lwd = 2)
  abline(h = target_rho, col = "gray40", lty = 2, lwd = 2)
}




genKappa = function(m, kappa0, sd_log_kappa, eps = 1e-12){
  if (sd_log_kappa < eps) return(rep(kappa0, m))
  z = rnorm(m)
  z = (z - mean(z)) / sd(z)
  exp(log(kappa0) + sd_log_kappa * z)
}



getTargetPhi1 = function(kappa0){
  1 / (1 + kappa0^2)
}
getTargetPhi2 = function(kappaA, kappaB){
  1 / sqrt((1 + kappaA^2) * (1 + kappaB^2))
}

genBifactorOmega = function(m,
                           kappa0,
                           sd_log_kappa,
                           h2 = 0.5){

  J = 2 * m
  indA = 1:m
  indB = (m + 1):J

  kappa_A = genKappa(m, kappa0, sd_log_kappa)
  kappa_B = genKappa(m, kappa0, sd_log_kappa)
  kappa = c(kappa_A, kappa_B)
  kappa_mean = prod(kappa)^(1/J)
  if (length(h2) == 1) h2 = rep(h2, J)
  stopifnot(length(h2) == J, all(h2 > 0), all(h2 < 1))

  lambda_g = sqrt(h2 / (1 + kappa^2))
  lambda_s = kappa * lambda_g
  psi = 1 - h2

  Lambda = matrix(0, nrow = J, ncol = 3)
  colnames(Lambda) = c("g", "sA", "sB")
  Lambda[, "g"] = lambda_g
  Lambda[indA, "sA"] = lambda_s[indA]
  Lambda[indB, "sB"] = lambda_s[indB]

  Omega = tcrossprod(Lambda) + diag(psi, J)


  list(
    Omega = Omega,
    Lambda = Lambda,
    kappa = kappa,
    psi = psi,
    indA = indA,
    indB = indB,
    kappa_mean = kappa_mean,
    sd_log_kappa_emp = sd(log(kappa))

  )
}