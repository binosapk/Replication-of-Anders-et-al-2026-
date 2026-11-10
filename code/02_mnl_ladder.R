# =====================================================================
# 02  Deterministic heterogeneity ladder (MNL in WTA space)  -  exogenous covariates first, attitudes second
# =====================================================================
# V_dev = bc * ( sum_k w_k x_k + cost/100 ),   V_sq = bc * ( w_sq + shifters ),   bc = b_cost * exp(g_cinc * income)
# Rung 0  baseline
# Rung 1  status quo depends on EXOGENOUS characteristics: age, gender, university degree, household income, exposure, env. organisation
# Rung 2  + cost sensitivity depends on income (a priori the most natural cost interaction)
# Rung 3  + ATTITUDES: general attitude (Q7), concern (Q10), visits nature (Q31), preservation (Q32)      <- endogenous, diagnostic
# Rung 4  + believes the results will be used (Q41)
prob_ladder <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  sq <- w_sq + g_age * age_c + g_fem * female + g_univ * univ + g_inc * inc_c + g_incm * inc_missing + g_exp * exposed + g_org * org +
        g_att * att_c + g_conc * conc_c + g_natv * natv_c + g_natp * natp_c + g_cons * conseq
  bc <- b_cost * exp(g_cinc * inc_c + g_cincm * inc_missing)
  V <- list()
  V[["dev"]] <- bc * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
                      w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
                      w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- bc * sq
  P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
  P <- apollo_panelProd(P, apollo_inputs, functionality); P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view", "b_cost")
set_exo  <- c("g_age", "g_fem", "g_univ", "g_inc", "g_incm", "g_exp", "g_org")
set_cost <- c("g_cinc", "g_cincm")
set_att  <- c("g_att", "g_conc", "g_natv", "g_natp")
set_cons <- "g_cons"
set_all  <- c(set_exo, set_cost, set_att, set_cons)
st0 <- c(setNames(c(rep(0, 11), 0.005), wn), setNames(rep(0, length(set_all)), set_all))

lad <- list()
for (s in c("solar", "wind")) {
  d <- dat[[s]]
  m0 <- fit_model(paste0("lad0_", s), d, prob_ladder, st0, fixed = set_all)
  st <- c(m0$estimate[wn], setNames(rep(0, length(set_all)), set_all))
  m1 <- fit_model(paste0("lad1_", s), d, prob_ladder, st, fixed = c(set_cost, set_att, set_cons))
  m2 <- fit_model(paste0("lad2_", s), d, prob_ladder, st, fixed = c(set_att, set_cons))
  m3 <- fit_model(paste0("lad3_", s), d, prob_ladder, st, fixed = set_cons)
  m4 <- fit_model(paste0("lad4_", s), d, prob_ladder, st)
  lad[[s]] <- list(m0 = m0, m1 = m1, m2 = m2, m3 = m3, m4 = m4)
}
cat("Model ladder\n")
for (s in c("solar", "wind")) {
  M <- lad[[s]]; nm <- c("0 baseline", "1 + exogenous shifters on status quo", "2 + cost x income", "3 + attitudes and values", "4 + believes results used")
  cat("\n=====", toupper(s), "=====\n")
  print(data.frame(rung = nm, LL = sapply(M, function(m) round(m$maximum, 1)), k = sapply(M, nfree), BIC = sapply(M, function(m) round(bic_of(m), 1)),
                   p_vs_previous = c(NA, lrp(M$m1, M$m0), lrp(M$m2, M$m1), lrp(M$m3, M$m2), lrp(M$m4, M$m3)), row.names = NULL), row.names = FALSE)
}

cat("\nStatus quo shifters: effect on WTA in NOK per household-year (positive = more compensation needed), robust t in brackets. Rung 4.\n")
eff <- function(m) { e <- m$estimate; se <- rse(m); k <- c(set_exo, set_att, set_cons); sg <- 100
  sprintf("%7.0f (%5.1f)", sg * e[k], e[k] / se[k]) }
lbl <- c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "exposed", "env. organisation",
         "attitude Q7 (+1)", "concern Q10 (+1)", "visit nature Q31 (+1)", "preserve nature Q32 (+1)", "believes results used")
print(data.frame(shifter = lbl, solar = eff(lad$solar$m4), wind = eff(lad$wind$m4)), row.names = FALSE)
cat("\nCost sensitivity: elasticity of the cost coefficient to +1 M NOK income (robust t):\n")
print(sapply(lad, function(M) { e <- M$m4$estimate; se <- rse(M$m4); sprintf("%.3f (%.1f)", e["g_cinc"], e["g_cinc"] / se["g_cinc"]) }))

