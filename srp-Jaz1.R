
readJazStudy1 = function(){
  dat = read.csv("_data/Jastrzebski-study1-data.csv")
  
  keep = c("FD_SIZ", "FD_ORIENT", "FD_BRI",
           "FM_SIZ", "FM_ORIENT", "FM_BRI",
           "IT_CIR", "IT_LIN", "IT_SQU",
           "TOJ_DIS", "TOJ_APP", "DD_EMPTY", "DD_FILLED")
  stopifnot(all(keep %in% names(dat)))
  tdat = dat[keep]
  scale(tdat)
}

semFitJaz = function(tdat){
  model = '
    FD =~ FD_SIZ + FD_ORIENT + FD_BRI
    FM =~ FM_SIZ + FM_ORIENT + FM_BRI
    IT =~ IT_CIR + IT_LIN + IT_SQU
    TD =~ TOJ_DIS + TOJ_APP + DD_EMPTY + DD_FILLED
  '
  
  fit = cfa(model, data = tdat,
            mimic   = "lavaan",
            missing = "ML",
            std.lv  = TRUE,
            std.ov  = FALSE,
            se      = "standard",
            bootstrap = 1000)
  fit
}

getJazLabels = function(){
  c(
    expression(FD[1], FD[2], FD[3],
               FM[1], FM[2], FM[3],
               IT[1], IT[2], IT[3],
               TD[1], TD[2], TD[3], TD[4])
  )
}

getJazColors = function(){
  list(
    cols = c(rep("black", 12), "black"),
    reds = c("#670D20", "#781121", "#891523", "#9F1A27", "#AE1E29",
             "#BB2B33", "#C43C3C", "#CD4E44", "#D5604C", "#DC6E58",
             "#E58469", "#ED9475", "#F3A483", "#F4B295", "#F6C1A4",
             "#F7CDB7", "#F9DCC7", "#FAE3D6", "#FCEDE4", "#FDF6F2"),
    whites = "#FDFFFE",
    blues = c("#F2F7FA", "#E6F1F8", "#DCEAF4", "#CFE3EF", "#BFDCEA",
              "#AFD3E6", "#A0CBE1", "#8FC3DD", "#7AB6D6", "#66AACF",
              "#529DC7", "#4492C2", "#3F84BD", "#3B7AB6", "#366FB0",
              "#3263A8", "#2A5595", "#234882", "#1D3A71", "#172F61")
  )
}


extractEstimatesJaz = function(fit){
  values_first_layer = subset(standardizedSolution(fit), op == "=~")
  vals1 = rev(values_first_layer$est.std)
  values_first_layer_ind = .map2idx(vals1)
  
  values_second_layer = subset(standardizedSolution(fit), 
                               op == "~~" & lhs != rhs & 
                               lhs %in% c("FD", "FM", "IT", "TD"))
  vals2 = rev(values_second_layer$est.std)
  vals2 = as.vector(sapply(vals2, function(x) min(x, 1)))
  values_second_layer_ind = .map2idx(vals2)
  
  list(
    vals1 = vals1,
    vals2 = vals2,
    vals1_ind = values_first_layer_ind,
    vals2_ind = values_second_layer_ind
  )
}

plotSemDiagramJaz = function(fit, shades, cols, labs, flag = F){
  estimates = extractEstimatesJaz(fit)  # FIXED: use Jaz-specific function
  
  factors = rev(c("Feature \nDiscrimination", "Feature \nMemory", 
                  "Inspection \nTime", "Temporal \nDiscrimination"))
  
  select_cols = c("FD Size", "FD Orient", "FD Bright", 
                  "FM Size", "FM Orient", "FM Bright", 
                  "IT Line", "IT Circle", "IT Square",
                  "Dissapear", "Appear", "Filled", "Empty")

  plot(NA, NA, xlim = c(0, 1.2), ylim = c(0, 2), axes = FALSE, 
       xlab = "", ylab = "")

  J = length(select_cols)
  task_pos.y = seq(0, 2, length.out = J + 3)
  task_pos.y = task_pos.y[-c(5, 9, 13)]
  task_pos.x = rep(0.2, J)
  task.length = c(0.1, 0.05)
  cex_val = 0.7
  
  rect(task_pos.x + task.length[1], task_pos.y + task.length[2],
       task_pos.x - task.length[1], task_pos.y - task.length[2], lwd = 2)
  text(task_pos.x, task_pos.y, rev(select_cols), cex = .8)
  
  labs_rev = rev(labs)
  for (i in 1:J) {
    mtext(labs_rev[i], side = 2, line = -1.75, at = task_pos.y[i],
          cex = .6, adj = 0.5, las = 1, col = cols[J - i + 1])
  }
  
  n_factors = 4
  factor_pos.y = c(task_pos.y[2] + (task_pos.y[3] - task_pos.y[2])/2,
                   task_pos.y[6], task_pos.y[9], task_pos.y[12])
  factor_pos.x = rep(0.7, n_factors)
  radius_x = 0.1
  
  for (i in 1:n_factors) {
    .drawEllipse(factor_pos.x[i], factor_pos.y[i], 
                 radius_y = 0.1, radius_x = radius_x)
  }
  text(factor_pos.x, factor_pos.y, factors, cex = .6, col = rev(unique(cols)))
  
  first_con = cbind(task_pos.x + task.length[1], task_pos.y)
  second_con = cbind(rep(factor_pos.x[1] - radius_x, 13),
                     c(factor_pos.y[1], rep(factor_pos.y, each = 3)))
  
  cols_line_first = sapply(1:13, function(x) shades[estimates$vals1_ind[x]])
  
  for (i in 1:13) {
    lines(c(first_con[i, 1], second_con[i, 1]),
          c(first_con[i, 2], second_con[i, 2]),
          lwd = 4, col = cols_line_first[i])
    
    mid_x = (first_con[i, 1] + second_con[i, 1]) / 2 + 0.05
    mid_y = (first_con[i, 2] + second_con[i, 2]) / 2 
    rect(mid_x - 0.019, mid_y + 0.05,
         mid_x - 0.082, mid_y - 0.05, col = "white", border = "white")
    text(mid_x - 0.05, mid_y, label = round(estimates$vals1[i], 3),
         cex = .9*cex_val, col = "black")
  }
  
  cols_line_second = sapply(1:6, function(x) shades[estimates$vals2_ind[x]])
  dists = c(0.15, 0.4, 0.15, 0.7, 0.4, 0.15)
  count = 1
  
  for (i in 1:n_factors) {
    for (j in 1:n_factors) {
      if (i > j) {
        .curveFunction(
          c(factor_pos.x[i] + radius_x, factor_pos.y[i]),
          c(factor_pos.x[j] + radius_x, factor_pos.y[j]),
          dists[count],
          col = cols_line_second[count]
        )
        
        mid_x = factor_pos.x[i] + dists[count]/2 + 0.1
        mid_y = (factor_pos.y[i] + factor_pos.y[j]) / 2 - 0.08
        rect(mid_x - 0.05, mid_y + 0.05,
             mid_x + 0.05, mid_y + 0.1, col = "white", border = "white")
        text(mid_x, mid_y + 0.08,
             labels = sprintf("%.3f", estimates$vals2[count]),
             cex = cex_val, col = "black")
        count = count + 1
      }
    }
  }
  
  if (flag == T){
    rect(0.82, 0.475, 0.925, 0.525, border = "indianred4", lwd = 3)
    rect(0.82, 1.575, 0.925, 1.625, border = "indianred1", lwd = 3)
  }
}

plotCorrelationMatrixJaz = function(R_hat, shades, cols, labs, flag = F){
  h = 14
  corrplot(unname(R_hat), 
           method = "color", 
           cl.pos = "b", 
           col = shades,
           mar = c(0, 2, 1, 0),
           tl.pos = "n",
           cl.cex = 1)
  
  for (i in seq_along(cols)) {
    mtext(labs[i], side = 2, line = 2.5, at = h - i,
          cex = .6, adj = 0.5, las = 1, col = cols[i])
    mtext(labs[14-i], side = 3, line = 0, at = h - i,
          cex = .6, adj = 0.5, las = 1, col = cols[14-i])
  }
  

  if (flag == T){
    rect(6.5, 0.5, 9.5, 4.5, border = "indianred4", lwd = 3)
    mean_outter_blk_IT_TD = round(mean(R_hat[7:9, 10:13], na.rm = TRUE), 3)
    lbl = bquote(bar(italic(r)) == .(mean_outter_blk_IT_TD))
    text(8, 2.5, label = lbl, cex = 1, col = "black")
    
    rect(3.5, 10.5, 6.5, 13.5, border = "indianred1", lwd = 3)
    mean_outter_DF_FM = round(mean(R_hat[1:3, 4:6], na.rm = TRUE), 3)
    lbl2 = bquote(bar(italic(r)) == .(mean_outter_DF_FM))
    text(5, 12, label = lbl2, cex = 1, col = "black")
  }
}

