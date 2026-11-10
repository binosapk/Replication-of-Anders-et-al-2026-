# =====================================================================
# 07  Hybrid choice (ICLV) with FIVE indicators and demographics in the structural equation
# =====================================================================
# Leeds [A] 3C/3D.  The survey measures one underlying "concern about nature / resistance to development" with five answers:
#   Q7  general attitude to development (reverse), Q10 concern about nature effects, Q31 importance of visiting untouched nature,
#   Q32 should untouched nature be kept at today's level (3 levels), Q33 perceived negative landscape impact of the technology in question.
# Structural equation:  LV = gamma' z + eta,  eta ~ N(0,1),   z = age, gender, degree, income, exposure, organisation (EXOGENOUS)
# Measurement:          each indicator ~ ordered logit( zeta_i * LV ),  zeta_conc > 0 anchors the sign (higher LV = more concern)
# Choice:               V_sq = b_cost * ( w_sq + lambda * LV ),   V_dev has  w_nop + lambda2 * LV  on the no-protection dummy
# Following the slides (3C, p. 24-29): the log-likelihood includes the indicators and must NOT be compared with choice-only models;
# the safe way to use the model is through the covariates z (reweighting a sample, a future population), not by moving LV by "+-1 SD".
d <- dat[[SAMPLE]]
d$nat_pres_f <- ifelse(is.na(d$nat_pres), 2, d$nat_pres); d$nat_pres_ok <- as.numeric(!is.na(d$nat_pres))
d$impact_f <- ifelse(is.na(d$impact_own), 4, d$impact_own); d$impact_ok <- as.numeric(!is.na(d$impact_own))
d$att_f <- ifelse(is.na(d$att), 4, d$att); d$conc_f <- ifelse(is.na(d$conc), 4, d$conc); d$natv_f <- ifelse(is.na(d$nat_visit), 4, d$nat_visit)

prob_H <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  # ---- thresholds built from positive increments (they stay ordered)
  ta2 <- ta1 + exp(ta2i); ta3 <- ta2 + exp(ta3i); ta4 <- ta3 + exp(ta4i); ta5 <- ta4 + exp(ta5i); ta6 <- ta5 + exp(ta6i)
  tc2 <- tc1 + exp(tc2i); tc3 <- tc2 + exp(tc3i); tc4 <- tc3 + exp(tc4i); tc5 <- tc4 + exp(tc5i); tc6 <- tc5 + exp(tc6i)
  tv2 <- tv1 + exp(tv2i); tv3 <- tv2 + exp(tv3i); tv4 <- tv3 + exp(tv4i); tv5 <- tv4 + exp(tv5i); tv6 <- tv5 + exp(tv6i)
  ti2 <- ti1 + exp(ti2i); ti3 <- ti2 + exp(ti3i); ti4 <- ti3 + exp(ti4i); ti5 <- ti4 + exp(ti5i); ti6 <- ti5 + exp(ti6i)
  tp2 <- tp1 + exp(tp2i)
  first <- (task == 1)
  P[["i_att"]]  <- apollo_ol(list(outcomeOrdered = att_f,        V = z_att * LV,  tau = list(ta1, ta2, ta3, ta4, ta5, ta6), coding = 1:7, rows = first, componentName = "i_att"), functionality)
  P[["i_conc"]] <- apollo_ol(list(outcomeOrdered = conc_f,       V = exp(lz_conc) * LV, tau = list(tc1, tc2, tc3, tc4, tc5, tc6), coding = 1:7, rows = first, componentName = "i_conc"), functionality)
  P[["i_natv"]] <- apollo_ol(list(outcomeOrdered = natv_f,       V = z_natv * LV, tau = list(tv1, tv2, tv3, tv4, tv5, tv6), coding = 1:7, rows = first, componentName = "i_natv"), functionality)
  P[["i_imp"]]  <- apollo_ol(list(outcomeOrdered = impact_f,     V = z_imp * LV,  tau = list(ti1, ti2, ti3, ti4, ti5, ti6), coding = 1:7, rows = (first & impact_ok == 1), componentName = "i_imp"), functionality)
  P[["i_natp"]] <- apollo_ol(list(outcomeOrdered = nat_pres_f,   V = z_natp * LV, tau = list(tp1, tp2), coding = 1:3, rows = (first & nat_pres_ok == 1), componentName = "i_natp"), functionality)
  # ---- choice component
  V <- list()
  V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
                          (w_nop + lambda2 * LV) * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
                          w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- b_cost * (w_sq + lambda * LV)
  P[["choice"]] <- apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V, componentName = "choice"), functionality)
  P <- apollo_combineModels(P, apollo_inputs, functionality)
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
rc_H <- function(apollo_beta, apollo_inputs) {
  r <- list(); r[["LV"]] <- g_age * age_c + g_fem * female + g_univ * univ + g_inc * inc_c + g_incm * inc_missing + g_exp * exposed + g_org * org + eta; return(r)
}
dr_H <- list(interDrawsType = "halton", interNDraws = 120, interNormDraws = "eta")

wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view", "b_cost")
m0 <- fit_model(paste0("hmnl_", SAMPLE), d, function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  V <- list(); V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
    w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) + w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- b_cost * w_sq
  P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
  P <- apollo_panelProd(P, apollo_inputs, functionality); P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P) },
  setNames(c(rep(0, 11), 0.005), wn))
thr <- function(p, n) setNames(c(-2, rep(0, n - 1)), c(paste0(p, "1"), paste0(p, 2:n, "i")))
st <- c(m0$estimate[wn], lambda = 40, lambda2 = 0,
        g_age = 0, g_fem = 0, g_univ = 0, g_inc = 0, g_incm = 0, g_exp = 0, g_org = 0,
        z_att = -1, lz_conc = 0, z_natv = 0.5, z_imp = 0.8, z_natp = 0.8,
        thr("ta", 6), thr("tc", 6), thr("tv", 6), thr("ti", 6), thr("tp", 2))
t0 <- Sys.time()
mH <- fit_model(paste0("hybrid_", SAMPLE), d, prob_H, st, draws = dr_H, randcoeff = rc_H, estimate_settings = list(estimationRoutine = "bfgs", maxIterations = 250))
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("Hybrid choice model, %s sample (%.0f min; %s)\nlog-likelihood %.1f (includes the 5 indicators), parameters %d\n", SAMPLE, mins, mH$message, mH$maximum, nfree(mH)))

e <- mH$estimate; V <- mH$robvarcov; se <- rse(mH)
cat("\nLoadings of the indicators on the latent variable 'concern about nature' (robust t)\n")
ld <- c(Q7_attitude = e[["z_att"]], Q10_concern = exp(e[["lz_conc"]]), Q31_visit = e[["z_natv"]], Q33_landscape_impact = e[["z_imp"]], Q32_preserve = e[["z_natp"]])
lt <- c(e[["z_att"]] / se[["z_att"]], e[["lz_conc"]] / se[["lz_conc"]], e[["z_natv"]] / se[["z_natv"]], e[["z_imp"]] / se[["z_imp"]], e[["z_natp"]] / se[["z_natp"]])
print(data.frame(indicator = names(ld), loading = round(ld, 2), t = round(lt, 1), row.names = NULL), row.names = FALSE)

cat("\nStructural equation: effect of respondent characteristics on the latent variable (LV has unit variance)\n")
gz <- c("g_age", "g_fem", "g_univ", "g_inc", "g_incm", "g_exp", "g_org")
lab <- c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "exposed", "env. organisation")
print(data.frame(characteristic = lab, gamma = round(e[gz], 3), t = round(e[gz] / se[gz], 1), row.names = NULL), row.names = FALSE)

cat("\nEffect of the latent variable on the choice (robust t): lambda (status quo, 100 NOK per LV unit)", round(e[["lambda"]], 1), sprintf("(%.1f)", e[["lambda"]] / se[["lambda"]]),
    "; lambda2 (no-protection dummy)", round(e[["lambda2"]], 1), sprintf("(%.1f)\n", e[["lambda2"]] / se[["lambda2"]]))
# effect of a characteristic on the status quo WTA through the LV = lambda * gamma  (delta method for a product of two parameters)
dm <- function(a, b) { g <- c(e[[b]], e[[a]]); v <- V[c(a, b), c(a, b)]; c(est = e[[a]] * e[[b]], se = sqrt(drop(t(g) %*% v %*% g))) }
cat("\nImplied effect on the status quo WTA through the latent variable, NOK per household-year  (delta method, robust)\n")
imp <- t(sapply(gz, function(g) dm("lambda", g) * c(100, 100)))
print(data.frame(characteristic = lab, effect_NOK = round(imp[, "est"]), se = round(imp[, "se"]), z = round(imp[, "est"] / imp[, "se"], 1), row.names = NULL), row.names = FALSE)
