library(lavaan)
library(psych)
library(parallel)


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


makeIcmCfaModel = function(J){
  if (J %% 2 != 0) stop("J must be even.")
  y = paste0("y", 1:J)
  f1 = paste(y[1:(J/2)], collapse = " + ")
  f2 = paste(y[(J/2 + 1):J], collapse = " + ")
  paste0(
    "f1 =~ ", f1, "\n",
    "f2 =~ ", f2, "\n"
  )
}




safeCfaLite = function(mod,
                       Sigma,
                       N,
                       tol = 1e-8,
                       keep_R_hat = FALSE,
                       engine = c("lavaan", "psych"),
                       estimator = c("ML", "ULS", "GLS"),
                       rhat_tol = 1.05,
                       blavaan_args = list(target = "stan", n.chains = 3, burnin = 500, sample = 500, adapt = 500),
                       psych_args = list()){
  engine = match.arg(engine)
  getPhi12 = function(Phi){
    if (is.null(Phi)) return(NA_real_)
    rn = rownames(Phi); cn = colnames(Phi)
    if (!is.null(rn) && !is.null(cn) && all(c("f1","f2") %in% rn) && all(c("f1","f2") %in% cn)) {
      return(unname(Phi["f1","f2"]))
    }
    if (nrow(Phi) >= 2 && ncol(Phi) >= 2) return(unname(Phi[1,2]))
    NA_real_
  }

  out = tryCatch({

    if (engine == "lavaan") {

      fit = lavaan::cfa(
        mod,
        sample.cov = Sigma,
        sample.nobs = N,
        std.lv = TRUE,
        estimator = estimator,
        se = "none",
        test = "none",
        warn = FALSE
      )

      converged = tryCatch(lavaan::lavInspect(fit, "converged"), error = function(e) FALSE)

      theta = tryCatch(lavaan::lavInspect(fit, "theta"), error = function(e) NULL)
      psi   = tryCatch(lavaan::lavInspect(fit, "psi"),   error = function(e) NULL)

      heywood = FALSE
      if (!is.null(theta)) heywood = heywood || any(diag(theta) < -tol)
      if (!is.null(psi))   heywood = heywood || any(diag(psi)   < -tol)

      corlv = tryCatch(lavaan::lavInspect(fit, "cor.lv"), error = function(e) NULL)
      phi12 = getPhi12(corlv)

      std = tryCatch(lavaan::lavInspect(fit, "std"), error = function(e) NULL)
      lambda_std_mat = if (is.null(std)) NULL else std$lambda

      R_hat = NULL
      if (isTRUE(keep_R_hat)) {
        Sig_hat = tryCatch(lavaan::lavInspect(fit, "sigma.hat"), error = function(e) NULL)
        if (!is.null(Sig_hat)) R_hat = stats::cov2cor(Sig_hat)
      }

      list(
        ok = TRUE,
        err = NA_character_,
        converged = isTRUE(converged),
        heywood = isTRUE(heywood),
        phi = phi12,
        lambda_std_mat = lambda_std_mat,
        R_hat = R_hat,
        engine = "lavaan"
      )

    } else if (engine == "psych") {

      R = tryCatch(cov2cor(Sigma), error = function(e) Sigma)

      args = c(
      list(model = mod, r = R, n.obs = N, all = FALSE, cor = "cor", orthog = FALSE),
        psych_args
      )

      CFA_fun = getFromNamespace("CFA", "psych")
      fit = do.call(CFA_fun, args)


      Phi = if (!is.null(fit$Phi)) as.matrix(fit$Phi) else NULL
      phi12 = getPhi12(Phi)

      Lambda = if (!is.null(fit$loadings)) as.matrix(fit$loadings) else NULL

      h2 = NULL
      if (!is.null(fit$communalities)) h2 = as.numeric(fit$communalities)
      if (is.null(h2) && !is.null(fit$communality)) h2 = as.numeric(fit$communality)

      uniq = NULL
      if (!is.null(fit$uniquenesses)) uniq = as.numeric(fit$uniquenesses)
      if (is.null(uniq) && !is.null(h2)) uniq = 1 - h2

      heywood = FALSE
      if (!is.null(h2))  heywood = heywood || any(h2 > 1 + tol)
      if (!is.null(uniq)) heywood = heywood || any(uniq < -tol)

      R_hat = NULL
      if (isTRUE(keep_R_hat) && !is.null(Lambda) && !is.null(Phi)) {
        if (is.null(uniq)) {
          h2_hat = diag(Lambda %*% Phi %*% t(Lambda))
          uniq = 1 - h2_hat
        }
        R_hat = Lambda %*% Phi %*% t(Lambda) + diag(uniq, nrow(Lambda))
        rownames(R_hat) = rownames(Lambda)
        colnames(R_hat) = rownames(Lambda)
      }

      list(
        ok = TRUE,
        err = NA_character_,
        converged = TRUE,                 
        heywood = isTRUE(heywood),
        phi = phi12,
        lambda_std_mat = Lambda,
        R_hat = R_hat,
        engine = "psych"
      )

    } else {

      args = c(
        list(model = mod, sample.cov = Sigma, sample.nobs = N, std.lv = TRUE),
        blavaan_args
      )
      fit = do.call(blavaan::bcfa, args)

      psrf = tryCatch(blavaan::blavInspect(fit, "psrf"), error = function(e) NULL)
      converged = if (is.null(psrf)) NA else (max(psrf, na.rm = TRUE) < rhat_tol)

      theta = tryCatch(lavaan::lavInspect(fit, "theta"), error = function(e) NULL)
      psi   = tryCatch(lavaan::lavInspect(fit, "psi"),   error = function(e) NULL)

      heywood = FALSE
      if (!is.null(theta)) heywood = heywood || any(diag(theta) < -tol)
      if (!is.null(psi))   heywood = heywood || any(diag(psi)   < -tol)

      corlv = tryCatch(lavaan::lavInspect(fit, "cor.lv"), error = function(e) NULL)
      phi12 = getPhi12(corlv)

      std = tryCatch(lavaan::lavInspect(fit, "std"), error = function(e) NULL)
      lambda_std_mat = if (is.null(std)) NULL else std$lambda

      R_hat = NULL
      if (isTRUE(keep_R_hat)) {
        Sig_hat = tryCatch(lavaan::lavInspect(fit, "sigma.hat"), error = function(e) NULL)
        if (!is.null(Sig_hat)) R_hat = stats::cov2cor(Sig_hat)
      }

      list(
        ok = TRUE,
        err = NA_character_,
        converged = converged,
        heywood = isTRUE(heywood),
        phi = phi12,
        lambda_std_mat = lambda_std_mat,
        R_hat = R_hat,
        engine = "blavaan",
        psrf_max = if (is.null(psrf)) NA_real_ else max(psrf, na.rm = TRUE)
      )
    }

  }, error = function(e){
    list(
      ok = FALSE,
      err = conditionMessage(e),
      converged = FALSE,
      heywood = NA,
      phi = NA_real_,
      lambda_std_mat = NULL,
      R_hat = NULL,
      engine = engine
    )
  })

  out
}


getKappaGrid = function(phi_lo, phi_hi, n){
  phi_seq = seq(phi_lo, phi_hi, length.out = n)
  sqrt((1/phi_seq) - 1)
}


runSimulation = function(kappa_target = getKappaGrid(0.01, .99, 100),
                         kappa_sd = seq(0.01, 3, length.out = 100),
                         h2_vals = 0.5,
                         Js = 10,
                         Ns = 2000,
                         n_sims = 10,
                         ncore = 2,
                         seed = 1234){

  set.seed(seed)

  grid_lav = expand.grid(
    kappa0 = kappa_target,
    sd_log_kappa = kappa_sd,
    h2 = h2_vals,
    J = Js,
    N = Ns,
    nsim = seq_len(n_sims),
    engine = c("lavaan"),
    estimator = c("ML", "ULS", "GLS"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  grid_psych = expand.grid(
    kappa0 = kappa_target,
    sd_log_kappa = kappa_sd,
    h2 = h2_vals,
    J = Js,
    N = Ns,
    nsim = seq_len(n_sims),
    engine = c("psych"),
    estimator = NA_character_,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  grid = rbind(grid_lav, grid_psych)

  unique_J = sort(unique(grid$J))
  mods = setNames(lapply(unique_J, makeIcmCfaModel), as.character(unique_J))
  var_names = setNames(lapply(unique_J, function(J) paste0("y", seq_len(J))), as.character(unique_J))

  seeds = sample.int(1e9, nrow(grid))

  sims = mclapply(seq_len(nrow(grid)), function(i){

    out = tryCatch({

      set.seed(seeds[i])

      J_i = grid$J[i]
      N_i = grid$N[i]
      m = J_i / 2
      indA = seq_len(m)
      indB = (m + 1):J_i

      kappa_0 = grid$kappa0[i]
      sd_log_kappa = grid$sd_log_kappa[i]
      h2 = grid$h2[i]
      engine_i = grid$engine[i]
      estimator_i = grid$estimator[i]

      gen = genBifactorOmega(
        m = m,
        kappa0 = kappa_0,
        sd_log_kappa = sd_log_kappa,
        h2 = h2
      )

      Omega = gen$Omega
      vn = var_names[[as.character(J_i)]]
      dimnames(Omega) = list(vn, vn)

      Lambda = gen$Lambda
      kappa = gen$kappa

      target_phi_0 = getTargetPhi1(kappa_0)
      target_phi_emp = getTargetPhi2(
        prod(kappa[indA])^(1/length(indA)),
        prod(kappa[indB])^(1/length(indB))
      )

      fit =
        if (engine_i == "lavaan") {
          safeCfaLite(mods[[as.character(J_i)]], Omega, N_i, engine = "lavaan", estimator = estimator_i)
      } else {
          safeCfaLite(mods[[as.character(J_i)]], Omega, N_i, engine = "psych")
      }
      phi_hat = fit$phi
      Lambda_ICM = fit$lambda_std_mat

      Lambdas = list(
        Lambda_true = Lambda,
        Lambda_ICM = Lambda_ICM
      )

      list(
        ok = isTRUE(fit$ok),
        err = fit$err,
        converged = isTRUE(fit$converged),
        heywood = isTRUE(fit$heywood),

        phi_hat = phi_hat,
        target_phi_0 = target_phi_0,
        target_phi_emp = target_phi_emp,

        sd_log_kappa = sd_log_kappa,
        sd_log_kappa_emp = gen$sd_log_kappa_emp,
        kappa_mean_emp = gen$kappa_mean,

        engine = fit$engine,
        Lambdas = Lambdas,

        J = J_i,
        N = N_i,
        kappa_0 = kappa_0,
        h2 = h2,
        omega = Omega
      )

    }, error = function(e){

      list(
        ok = FALSE,
        err = conditionMessage(e),
        converged = FALSE,
        heywood = NA,

        phi_hat = NA_real_,
        target_phi_0 = NA_real_,
        target_phi_emp = NA_real_,

        sd_log_kappa = grid$sd_log_kappa[i],
        sd_log_kappa_emp = NA_real_,
        kappa_mean_emp = NA_real_,

        engine = grid$engine[i],

        J = grid$J[i],
        N = grid$N[i],
        kappa_0 = grid$kappa0[i],
        h2 = grid$h2[i],
        omega = NA
      )
    })

    out

    }, mc.cores = ncore, mc.set.seed = FALSE)

  ok = vapply(sims, function(x) isTRUE(x$ok), logical(1))
  converged = vapply(sims, function(x) isTRUE(x$converged), logical(1))
  heywood = vapply(sims, function(x) isTRUE(x$heywood), logical(1))
  err = vapply(sims, function(x) x$err, character(1))

  phi_hat = vapply(sims, function(x) x$phi_hat, numeric(1))
  target_phi_0 = vapply(sims, function(x) x$target_phi_0, numeric(1))
  target_phi_emp = vapply(sims, function(x) x$target_phi_emp, numeric(1))

  sd_log_kappa = vapply(sims, function(x) x$sd_log_kappa, numeric(1))
  sd_log_kappa_emp = vapply(sims, function(x) x$sd_log_kappa_emp, numeric(1))
  kappa_mean_emp = vapply(sims, function(x) x$kappa_mean_emp, numeric(1))

  kappa_0 = vapply(sims, function(x) x$kappa_0, numeric(1))
  h2 = vapply(sims, function(x) x$h2, numeric(1))
  J = vapply(sims, function(x) x$J, numeric(1))
  N = vapply(sims, function(x) x$N, numeric(1))

  keep = ok & converged & !heywood & is.finite(phi_hat) & is.finite(target_phi_emp)

  Lambdas = lapply(sims, function(x) x$Lambdas)
  omega = lapply(sims, function(x) x$omega)

  list(
    grid = grid,
    seeds = seeds,
    sims = sims,
    Lambda_list = Lambdas,

    ok = ok,
    converged = converged,
    heywood = heywood,
    err = err,
    keep = keep,

    phi_hat = phi_hat,
    target_phi_0 = target_phi_0,
    target_phi_emp = target_phi_emp,
    deviation = phi_hat - target_phi_emp,

    sd_log_kappa = sd_log_kappa,
    sd_log_kappa_emp = sd_log_kappa_emp,
    kappa_mean_emp = kappa_mean_emp,

    kappa_0 = kappa_0,
    h2 = h2,
    J = J,
    N = N,
    omega = omega
  )
}



runSimulationCached = function(result_file = "_results/res_sim_vio.RDS",
                               force_rerun = FALSE,
                               ...){

  if (!force_rerun && file.exists(result_file)) {
    message("Loading cached results from: ", result_file)
    return(readRDS(result_file))
  }

  message("Running simulation...")
  results = runSimulation(...)

  result_dir = dirname(result_file)
  if (!dir.exists(result_dir)) dir.create(result_dir, recursive = TRUE)

  saveRDS(results, result_file)
  message("Results saved to: ", result_file)

  results
}

