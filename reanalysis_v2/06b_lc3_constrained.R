# =====================================================================
# 06b  Three-class model with a cost-RESPONSIVE trader class (repair of a degenerate solution in the wind sample)
# =====================================================================
# In step 06 the wind LC3 solution let the 'trader' class drift to b_cost ~ 0.0002 with WTAs in the millions: in WTA space this is a class that
# ignores cost (utility-space weights b_cost * w stay finite while w -> infinity), i.e. a second cost-ignorer class. Leeds [A] 2A (slide 17):
# "with more than two classes the risk increases of some classes collapsing together ... fine tuning can help, e.g. merging or restricting parameters".
# Remedy: keep the trader class cost-responsive by fixing its cost coefficient at the value of the two-class model (LC2-DA); everything else is free.
# For solar the unrestricted three-class model is valid; it is re-estimated here with the same restriction so that both samples are treated alike.
e06 <- new.env(); load(file.path(V2, "state_06_both.RData"), envir = e06); lc <- e06$lc; lc3 <- e06$lc3
prob_lc3 <- e06$prob_lc3; lcp3 <- e06$lcp3
fx_cost <- function(s) lc[[s]]$lDA$estimate[["b_cost"]]
res <- list()
for (s in c("solar", "wind")) {
  st <- lc3[[s]]$estimate
  wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view")
  st[wn] <- lc[[s]]$lDA$estimate[wn]; st[["b_cost"]] <- fx_cost(s)
  res[[s]] <- fit_model(paste0("lc3c_", s), dat[[s]], prob_lc3, st, fixed = "b_cost", lcpars = lcp3)
}
cat("Three classes with a cost-responsive trader class (b_cost fixed at its two-class value)\n")
for (s in c("solar", "wind")) {
  m <- res[[s]]; l2 <- lc[[s]]$lDA
  cat(sprintf("\n===== %s =====\nLC2-DA: LL %.1f, k %d, BIC %.1f   |   LC3 (restricted): LL %.1f, k %d, BIC %.1f   (unrestricted LC3: LL %.1f)\n", toupper(s),
              l2$maximum, nfree(l2), bic_of(l2), m$maximum, nfree(m), bic_of(m), lc3[[s]]$maximum))
  cat(sprintf("LR test LC3 vs LC2-DA: %.1f on %d df, p = %s\n", 2 * (m$maximum - l2$maximum), nfree(m) - nfree(l2), format(lrp(m, l2), digits = 2)))
  e <- m$estimate; se <- rse(m); r <- dat[[s]][!duplicated(dat[[s]]$id), ]
  v2 <- e[["d2"]] + e[["h2_age"]] * r$age_c + e[["h2_fem"]] * r$female + e[["h2_att"]] * r$att_c + e[["h2_conc"]] * r$conc_c + e[["h2_natp"]] * r$natp_c + e[["h2_cons"]] * r$conseq + e[["h2_icost"]] * r$imp_cost_c
  v3 <- e[["d3"]] + e[["h3_age"]] * r$age_c + e[["h3_fem"]] * r$female + e[["h3_att"]] * r$att_c + e[["h3_conc"]] * r$conc_c + e[["h3_natp"]] * r$natp_c + e[["h3_cons"]] * r$conseq + e[["h3_icost"]] * r$imp_cost_c
  den <- 1 + exp(v2) + exp(v3); sh <- cbind(trader = 1 / den, cost_ignorer = exp(v2) / den, sq_only = exp(v3) / den)
  cat("Class shares (unweighted | population-weighted):\n")
  wtv <- as.numeric(r$wt)
  print(round(rbind(unweighted = colMeans(sh), weighted = colSums(sh * wtv) / sum(wtv)), 3))
  wn2 <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view"); sg <- c(100, rep(-100, 10))
  cat("Trader-class WTA, NOK (robust s.e.):\n")
  print(data.frame(attribute = sub("w_", "", wn2), WTA = round(sg * e[wn2]), se = round(100 * se[wn2])), row.names = FALSE)
  en <- grep("^e_", names(e), value = TRUE)
  cat("Cost-ignorer class, utility weights (robust t):\n"); print(data.frame(par = sub("^e_", "", en), weight = round(e[en], 2), t = round(e[en] / se[en], 1)), row.names = FALSE)
  cat(sprintf("Status-quo-only class: P(status quo) = %.3f.  Cost-ignorers: P(status quo) at the reference scenario = %.3f\n", plogis(e[["a_sq2"]]), plogis(e[["a_sq_c"]])))
  cat(sprintf("Membership of the cost-ignorer class, odds ratio per +1 point of STATED cost importance (Q25): %.2f (z %.1f);  attitude Q7: %.2f (z %.1f)\n",
              exp(e[["h2_icost"]]), e[["h2_icost"]] / se[["h2_icost"]], exp(e[["h2_att"]]), e[["h2_att"]] / se[["h2_att"]]))
}
