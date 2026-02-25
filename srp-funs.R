runJags = function(mod, dat, pars, nchains = 2, niter = 3000, nburnin = 1000){
  result = jags(model.file = textConnection(mod),
                data = dat,
                n.chains = nchains,
                n.iter = niter,
                n.burnin = nburnin,
                n.thin = 1,
                parameters.to.save = pars)
  return(result)
}

quickCorPlot = function(mat,
                        title = NULL,
                        type = c("upper","full"),
                        show_numbers = TRUE,
                        digits = 2,
                        tl.cex = 0.85,
                        number.cex = 0.7,
                        mar = c(0,0,2,0),
                        reorder = FALSE) {
  if (!requireNamespace("corrplot", quietly = TRUE))
    stop("Package 'corrplot' is required.")

  type = match.arg(type)

  pal = rev(grDevices::colorRampPalette(
    c("#3B4CC0","#84A7F2","#FFFFFF","#F7A889","#B40426")
  )(200))

  corrplot(as.matrix(mat),
    method       = "color",
    type         = type,
    col          = pal,
    tl.col       = "black",
    tl.cex       = tl.cex,
    tl.srt       = 45,
    addCoef.col  = if (show_numbers) "black" else NULL,
    number.cex   = number.cex,
    number.digits= digits,
    diag         = FALSE,
    cl.pos       = "r",        
    cl.ratio     = 0.15,
    addgrid.col  = rgb(0,0,0,0.08),
    mar          = mar,
    order        = if (reorder) "hclust" else "original",
    hclust.method= "complete")
  if (!is.null(title)) mtext(title, 3, line = 0.5, cex = 1.1)
}

getSigmaFullLavan = function(fit){
  S = lavInspect(fit, "sigma")
  J = nrow(S) 
  Out = matrix(0, J, J)
  Out[lower.tri(Out)] = S[lower.tri(S)]
  diag(Out) = diag(S)
  Out = Out + t(Out) - diag(diag(Out))
  return(list(Sigma_hat = Out, R_hat = cov2cor(Out)))
}

.getMeInd = function(value, cols) {
  if (is.na(value)) return("black")
  if (value < 0.2) return(cols[1])
  else if (value < 0.4) return(cols[2])
  else if (value < 0.6) return(cols[3])
  else if (value < 0.8) return(cols[4])
  else if (value < .99) return(cols[5])
  else return(cols[6])
}

.drawEllipse = function(x_center, y_center, radius_x = 0.25, radius_y = 0.2, l = 2){
  theta = seq(0, 2 * pi, length.out = 100)
  x = radius_x * cos(theta) + x_center
  y = radius_y * sin(theta) + y_center
  polygon(x, y, border = "black", col = "white", lwd = l)
}

.curveFunction = function(start, end, distance, mid_y=0, col = "gray", lty = 1, lwd = 4) {
  if (mid_y == 0){
    mid_y = (start[2] + end[2]) / 2
  }
  control_points = rbind(start, c(start[1] + distance, mid_y), end)
  
  bezier_curve = function(t, control_points) {
    (1 - t)^2 * control_points[1,] + 2 * (1 - t) * t * control_points[2,] + t^2 * control_points[3,]
  }
  
  t_values = seq(0, 1, length.out = 100)
  curve_points = t(sapply(t_values, function(t) bezier_curve(t, control_points)))
  
  lines(curve_points, col = col, lwd = lwd, lty = lty)
}


.myCorPlot = function(cormat, limits = c(-1,1), main = "", coef = FALSE, tol = 0.1, mar = c(0,0,2,0), col = NULL, cl = "n", outline = FALSE){
  lims = sort(limits)
  rng = range(cormat, na.rm = TRUE)
  if (rng[1] < lims[1] - tol || rng[2] > lims[2] + tol) lims = range(lims, rng) else cormat = pmin(pmax(cormat, lims[1]), lims[2])
  J = nrow(cormat)
  rownames(cormat) = rep("", J); colnames(cormat) = rep("", J)
  if (!coef) corrplot::corrplot(cormat, method = "color", cl.pos = cl, col.lim = lims, main = main, mar = mar, col = col, outline = outline) else
    corrplot::corrplot(cormat, method = "color", addCoef.col = "black", cl.pos = cl, col.lim = lims, main = main, mar = mar, col = col, outline = outline)
}


readTsukStudy1 = function(flip_pos = T, scale_dat = T){
  dat = read.csv("_data/Tsukahara-2020-study1-data.csv")
  select_cols = c("RAPM", "LetterSets", "NumberSeries", 
                "OSpan", "SymSpan", "RotSpan",
                "Antisaccade", "FlankerEffect", "StroopEffect",
                "PitchThreshold.avg", "LoudThreshold.avg", 
                "LineThreshold.avg", "CircleThreshold.avg")
  select_cols = c("RAPM", "LetterSets", "NumberSeries", 
              "OSpan", "SymSpan", "RotSpan",
              "Antisaccade", "FlankerEffect", "StroopEffect",
              "PitchThreshold.last4rev", "LoudThreshold.last4rev", 
              "LineThreshold.last4rev", "CircleThreshold.last4rev")
  tdat = dat[,select_cols]
  tdat = tdat[complete.cases(tdat), ]
  if (flip_pos == T){tdat[, 8:13] = tdat[, 8:13]*(-1)}
  if (scale_dat == T){tdat = scale(tdat)}
  return(tdat)
}

modelMiss = function(CM1, CM2, truth = T, lns  = 0){
  par(mfrow = c(1,3))
  lab1 = "Sample"
  if (truth == T){lab1 = "Truth"}
  .myCorPlot(CM1)
  mtext(lab1, line = lns)
  .myCorPlot(CM2)
  mtext("Model", line = lns)
  .myCorPlot(CM1-CM2, coef=T)
  mtext("Miss", line = lns)
}


makeEmptyDat = function(I, J, L){
  sub   = rep(1:I, each = J*L)
  task  = rep(rep(1:J, each = L), I)
  trial = rep(1:L, I*J)
  data.frame(sub = sub, task = task, trial = trial)
}

makeTheta = function(I, J, t.b, t.c = 1, t.nu = 0, jitter = 1e-8){
  t.lam = as.matrix(t.b)                  
  if(nrow(t.lam) != J) stop("t.b must be J x K with nrow = J.")
  
  if(length(t.c) == 1) t.c = rep(t.c, J)  
  if(length(t.c) != J) stop("t.c must be length 1 or J.")
  
  comm = rowSums(t.lam^2)                     
  psi  = pmax(t.c - comm, 0)                
  Psi  = diag(psi, J, J)              
  
  Cov = t.lam %*% t(t.lam) + Psi
  
  ev_min = min(eigen(Cov, symmetric = TRUE, only.values = TRUE)$values)
  if(ev_min <= 0){                          
    Cov = Cov + diag(abs(ev_min) + jitter, J, J)
    psi = pmax(diag(Cov) - diag(t.lam %*% t(t.lam)), 0)
    Psi = diag(psi, J, J)
  }
  
  rho   = cov2cor(Cov)
  Nu    = rep(t.nu, J)
  theta = mvtnorm::rmvnorm(I, Nu, Cov)
  
  list(Nu = Nu, lam = t.lam, comm = comm, psi = psi,
       Psi = Psi, Cov = Cov, rho = rho, val = theta)
}


addTrialNoise = function(dat, Theta, sigma_trial){
  idx = cbind(dat$sub, dat$task)
  mu  = Theta[idx]                    
  dat$y = mu + rnorm(nrow(dat), 0, sigma_trial)
  dat
}

makeThetaFromRho = function(I, Rho, sd = 1, nu = 0, jitter = 1e-8){
  J = nrow(Rho)
  if (ncol(Rho) != J) stop("Rho must be JxJ.")
  if (any(abs(Rho - t(Rho)) > 1e-8)) stop("Rho not symmetric.")
  if (any(abs(diag(Rho) - 1) > 1e-8)) stop("diag(Rho) must be 1.")
  ev_min = min(eigen(Rho, symmetric = TRUE, only.values = TRUE)$values)
  if (ev_min <= 0){
    Rho = Rho + diag(abs(ev_min) + jitter, J)
    D = diag(1 / sqrt(diag(Rho)), J)
    Rho = D %*% Rho %*% D
  }
  if (length(sd) == 1) sd = rep(sd, J)
  if (length(sd) != J) stop("sd length mismatch.")
  if (length(nu) == 1) nu = rep(nu, J)
  if (length(nu) != J) stop("nu length mismatch.")
  Cov = diag(sd, J) %*% Rho %*% diag(sd, J)
  theta = mvtnorm::rmvnorm(I, mean = nu, sigma = Cov)
  list(Nu = nu, Cov = Cov, rho = Rho, val = theta)
}

popThetaCors = function(theta, probs = c(0.025, 0.975), include_post = FALSE){
  M = dim(theta)[1]
  J = dim(theta)[3]
  theta_cors = array(NA_real_, dim = c(M, J, J))
  for (m in 1:M){
    theta_cors[m,,] = cor(theta[m,,], use = "pairwise.complete.obs")
  }
  out = list(
    avg_cor = apply(theta_cors, c(2,3), mean,   na.rm = TRUE),
    med_cor = apply(theta_cors, c(2,3), median, na.rm = TRUE),
    sd_cor  = apply(theta_cors, c(2,3), sd,     na.rm = TRUE),
    ci_lo   = apply(theta_cors, c(2,3), quantile, probs = probs[1], na.rm = TRUE),
    ci_hi   = apply(theta_cors, c(2,3), quantile, probs = probs[2], na.rm = TRUE)
  )
  if (include_post) out$post_cor = theta_cors
  return(out)
}

popSigCors = function(Sig, p = F, probs = c(0.025, 0.975), include_post = F){
  M = dim(Sig)[1]
  J = dim(Sig)[3]
  cor_vals = array(dim=c(M,J,J))
  for (m in 1:M){
    Cov = Sig[m,,]
    if (p == T){
      Cov = solve(Cov)
    }
    cor_vals[m,,]=cov2cor(Cov)
  }
  out = list(
    avg_cor = apply(cor_vals, c(2,3), mean, na.rm = TRUE),
    med_cor = apply(cor_vals, c(2,3), median, na.rm = TRUE),
    sd_cor  = apply(cor_vals, c(2,3), sd, na.rm = TRUE),
    ci_lo   = apply(cor_vals, c(2,3), quantile, probs = probs[1], na.rm = TRUE),
    ci_hi   = apply(cor_vals, c(2,3), quantile, probs = probs[2], na.rm = TRUE)
  )
  if (include_post) out$post_cor = cor_vals
  return(out)
}

 
reportProp = function(lam_std, del2_std, var_unit = F, some_names = NULL){
  dims = dim(lam_std)
  M = dims[1]
  J = dims[2]
  D = dims[3]
  if (is.null(some_names)) some_names = 1:J
  lam_sum = matrix(0, J, D)
  uni_prop_task_mat = matrix(0, M, J)
  fac_prop_total_mat = matrix(0, M, D)
  uni_prop_total_vec = numeric(M)
  
  for (m in 1:M) {
    tlam = lam_std[m,,]
    tdel = del2_std[m,]
    
    tot_com_task = rowSums(tlam^2)
    tot_task_j   = tot_com_task + tdel
    tot_com_fac  = colSums(tlam^2)
    tot_uni      = sum(tdel)
    tot_var      = sum(tot_task_j)
    
    lam_sum = lam_sum + tlam
    uni_prop_task_mat[m,] = tdel / tot_task_j
    fac_prop_total_mat[m,] = tot_com_fac / tot_var
    uni_prop_total_vec[m] = tot_uni / tot_var
  }
  lam_mean = lam_sum / M
  
  if (var_unit==T){lam_mean = lam_mean^2}
  
  uni_prop_task = round(colMeans(uni_prop_task_mat), 3)
  fac_prop_total = round(colMeans(fac_prop_total_mat), 3)
  uni_prop_total = round(mean(uni_prop_total_vec), 3)
  
  out = cbind(lam_mean, uni_prop_task)
  out = rbind(out, c(fac_prop_total, uni_prop_total))
  out = round(out, 3)
  
  rownames(out) = c(some_names, "Prop. Var.")
  colnames(out) = c(paste0("F", 1:D), "Unique Var.")
  out
}

prepCrossPlot = function(lambda, del2 = NULL){
  dims = dim(lambda)
  M = dims[1]
  J = dims[2]
  D = dims[3]
  CMs_len = D + 2
  CMs = lapply(1:CMs_len, function(x) array(data = NA, dim = c(M, J, J)))
  
  for (m in 1:M){
    if (!is.null(del2)){
      tdels = diag(del2[m,])
    } else{
      tdels = diag(1-diag(tcrossprod(lambda[m,,])))
    }
    CMs[[1]][m,,] = tcrossprod(lambda[m,,]) + tdels
    for (d in 1:D){
      CMs[[(d+1)]][m,,] = tcrossprod(lambda[m,,d])
    }
    CMs[[CMs_len]][m,,] = tdels
  }
  
  CM_avg = list()
  CM_avg[[1]] = apply(CMs[[1]], 2:3, mean) 
  for (d in 1:D){
    CM_avg[[(d+1)]] = apply(CMs[[(d+1)]], 2:3, mean)
  }
  CM_avg[[CMs_len]] = apply(CMs[[CMs_len]], 2:3, mean) 
  
  return(CM_avg)
}


plotCompFac = function(CMs){
  com_num = length(CMs)


  diag(CMs[[1]]) = 1
  
  vlim = range(CMs, na.rm = TRUE)
    
  
  
  par(mfrow = c(1,com_num))
  .myCorPlot(CMs[[1]], main = "")
  ln = -2
  mtext("Correlation", side = 3, line = -ln, cex = 0.75, font = 2, adj = .5)
  mtext("=", side = 3, line = -ln, cex = 0.75, font = 2, adj = 1)
  for (i in 2:(com_num-1)){
    .myCorPlot(CMs[[i]], main = "")
    mtext(paste0("Factor", i-1), side = 3, line = -ln, cex = 0.75, font = 2, adj = .5)
    mtext("+", side = 3, line = -ln, cex = 0.75, font = 2, adj = 1)
  }
  .myCorPlot(CMs[[com_num]], main = "")
  mtext("Residuals", side = 3, line = -ln, cex = 0.75, font = 2, adj = .5)
}

itemRMS = function(R_true, R_hat){
  Res = R_true - R_hat
  sapply(1:nrow(Res), function(i) sqrt(mean(Res[i, -i]^2)))
}


.skipOnWarning = function(expr){
  tryCatch(withCallingHandlers(force(expr), warning = function(w) stop(w)), error = function(e) NULL)
}


rmseUpperTri = function(matrix_a, matrix_b) {
  inds = upper.tri(matrix_a, diag = FALSE)
  sqrt(mean((matrix_a[inds] - matrix_b[inds])^2, na.rm = TRUE))
}


.par = function(x = 1, y = 1, mar = c(4,4,2,2)) {
  par(mfrow = c(x,y), mar = mar)
  invisible(NULL)
}

.alpha = function(col, alpha) {
  rgb_vals = col2rgb(col) / 255
  rgb(rgb_vals[1, ], rgb_vals[2, ], rgb_vals[3, ], alpha = alpha)
}

.map2idx = function(vals) {
  idx = round((vals + 1) / 0.05) + 1      
  as.integer(pmax(1, pmin(41, idx)))      
}