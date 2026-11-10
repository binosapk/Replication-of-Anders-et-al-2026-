# =====================================================================
# 06  Latent class family: who never trades, who ignores cost, what do traders want?  (+ a formal solar-vs-wind test)
# =====================================================================
# Leeds [W] 5A / [A] 2A: class allocation with covariates (apollo_classAlloc), MNL kernel in each class.
# Anders: 2 classes with CONSTANT shares (traders / status-quo-only with a random constant).  Here:
#   LC2-0   2 classes, constant shares                                  (Anders's structure, MNL kernel)
#   LC2-D   2 classes, shares depend on EXOGENOUS characteristics        (age, gender, degree, income, exposure, organisation)
#   LC2-DA  2 classes, shares depend on exogenous characteristics + attitudes (Q7, Q10, Q31, Q32) + belief that results are used (Q41)
#   LC3     3 classes: traders / cost-IGNORERS (stated attribute non-attendance: utility without cost) / status-quo-only
# Class 1 (traders) is in WTA space; class "sq-only" ignores all attributes: P(status quo) = logistic(a).
zd <- c("age_c", "female", "univ", "inc_c", "inc_missing", "exposed", "org")
lcp2 <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  V <- list(trader = 0, sqonly = d_cls + gd_age * age_c + gd_fem * female + gd_univ * univ + gd_inc * inc_c + gd_incm * inc_missing + gd_exp * exposed + gd_org * org +
                                   ga_att * att_c + ga_conc * conc_c + ga_natv * natv_c + ga_natp * natp_c + gc_cons * conseq)
  lcpars[["pi_values"]] <- apollo_classAlloc(list(classes = c(trader = 1, sqonly = 2), utilities = V))
  return(lcpars)
}
trader_V <- function() NULL   # (kept simple: utilities are written out inline below)
prob_lc2 <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list(); mnl <- list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice)
  V <- list()
  V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
                          w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
                          w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- b_cost * w_sq
  mnl$utilities <- V; mnl$componentName <- "trader"
  P[["trader"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
  mnl$utilities <- list(dev = 0, sq = a_sq2); mnl$componentName <- "sqonly"
  P[["sqonly"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
  P[["model"]] <- apollo_lc(list(inClassProb = P, classProb = pi_values), apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
set_gd <- c("gd_age", "gd_fem", "gd_univ", "gd_inc", "gd_incm", "gd_exp", "gd_org")
set_ga <- c("ga_att", "ga_conc", "ga_natv", "ga_natp", "gc_cons")
wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view", "b_cost")

lc <- list()
for (s in c("solar", "wind")) {
  d <- dat[[s]]
  base <- c(setNames(c(rep(0, 11), 0.005), wn), a_sq2 = 3, d_cls = 0, setNames(rep(0, 12), c(set_gd, set_ga)))
  base[["w_sq"]] <- 0
  l0  <- fit_model(paste0("lc20_", s), d, prob_lc2, base, fixed = c(set_gd, set_ga), lcpars = lcp2)
  st1 <- l0$estimate; st1[c(set_gd, set_ga)] <- 0
  lD  <- fit_model(paste0("lc2D_", s), d, prob_lc2, st1, fixed = set_ga, lcpars = lcp2)
  lDA <- fit_model(paste0("lc2DA_", s), d, prob_lc2, lD$estimate, lcpars = lcp2)
  lc[[s]] <- list(l0 = l0, lD = lD, lDA = lDA)
}

# ---------------- 3 classes -------------
zc <- c("age_c", "female", "att_c", "conc_c", "natp_c", "conseq", "imp_cost_c")
lcp3 <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  V <- list(trader = 0,
            costign = d2 + h2_age * age_c + h2_fem * female + h2_att * att_c + h2_conc * conc_c + h2_natp * natp_c + h2_cons * conseq + h2_icost * imp_cost_c,
            sqonly  = d3 + h3_age * age_c + h3_fem * female + h3_att * att_c + h3_conc * conc_c + h3_natp * natp_c + h3_cons * conseq + h3_icost * imp_cost_c)
  lcpars[["pi_values"]] <- apollo_classAlloc(list(classes = c(trader = 1, costign = 2, sqonly = 3), utilities = V))
  return(lcpars)
}
prob_lc3 <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list(); mnl <- list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice)
  V <- list()
  V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
                          w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
                          w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- b_cost * w_sq
  mnl$utilities <- V; mnl$componentName <- "trader"
  P[["trader"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
  # cost-ignorers: the same attributes, but NO cost term (stated attribute non-attendance to cost); utility units
  V <- list()
  V[["dev"]] <- e_el38 * (el_alt1 == 38) + e_el57 * (el_alt1 == 57) + e_pl3 * (plants_alt1 == 3) + e_pl6 * (plants_alt1 == 6) +
                e_nop * (wild_alt1 == 0) + e_wil * (wild_alt1 == 1) + e_d1 * (dist_alt1 == 1) + e_d2 * (dist_alt1 == 2) +
                e_d3 * (dist_alt1 == 3) + e_view * (view_alt1 == 2)
  V[["sq"]] <- a_sq_c
  mnl$utilities <- V; mnl$componentName <- "costign"
  P[["costign"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
  mnl$utilities <- list(dev = 0, sq = a_sq2); mnl$componentName <- "sqonly"
  P[["sqonly"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
  P[["model"]] <- apollo_lc(list(inClassProb = P, classProb = pi_values), apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
hn <- c(paste0("h2_", c("age", "fem", "att", "conc", "natp", "cons", "icost")), paste0("h3_", c("age", "fem", "att", "conc", "natp", "cons", "icost")))
en <- c("e_el38", "e_el57", "e_pl3", "e_pl6", "e_nop", "e_wil", "e_d1", "e_d2", "e_d3", "e_view", "a_sq_c")
lc3 <- list()
for (s in c("solar", "wind")) {
  d <- dat[[s]]; e0 <- lc[[s]]$lDA$estimate
  mnl_pref <- c(e_el38 = -0.1, e_el57 = -0.1, e_pl3 = -0.1, e_pl6 = -0.2, e_nop = -0.5, e_wil = -0.2, e_d1 = -0.3, e_d2 = -0.2, e_d3 = -0.1, e_view = -0.3, a_sq_c = 0.5)
  st <- c(e0[c(wn, "a_sq2")], mnl_pref, d2 = -1, d3 = 0, setNames(rep(0, 14), hn))
  lc3[[s]] <- try(fit_model(paste0("lc3_", s), d, prob_lc3, st, lcpars = lcp3), silent = TRUE)
}

# ---------------- report -------------
cat("Model comparison\n")
for (s in c("solar", "wind")) {
  M <- lc[[s]]; m3 <- lc3[[s]]; ok3 <- !inherits(m3, "try-error")
  nm <- c("LC2-0  constant shares", "LC2-D  shares: exogenous covariates", "LC2-DA shares: + attitudes, Q41", if (ok3) "LC3    + cost-ignorer class")
  ms <- c(list(M$l0, M$lD, M$lDA), if (ok3) list(m3))
  cat("\n=====", toupper(s), "=====\n")
  print(data.frame(model = nm, LL = sapply(ms, function(m) round(m$maximum, 1)), k = sapply(ms, nfree), BIC = sapply(ms, function(m) round(bic_of(m), 1)),
                   p_vs_previous = c(NA, lrp(M$lD, M$l0), lrp(M$lDA, M$lD), if (ok3) lrp(m3, M$lDA)), row.names = NULL), row.names = FALSE)
}
cat("\nPopulation-weighted vs unweighted share in the status-quo-only class, model LC2-DA\n")
shr <- function(m, s) { e <- m$estimate; r <- dat[[s]][!duplicated(dat[[s]]$id), ]
  p <- plogis(e[["d_cls"]] + e[["gd_age"]] * r$age_c + e[["gd_fem"]] * r$female + e[["gd_univ"]] * r$univ + e[["gd_inc"]] * r$inc_c + e[["gd_incm"]] * r$inc_missing +
              e[["gd_exp"]] * r$exposed + e[["gd_org"]] * r$org + e[["ga_att"]] * r$att_c + e[["ga_conc"]] * r$conc_c + e[["ga_natv"]] * r$natv_c + e[["ga_natp"]] * r$natp_c + e[["gc_cons"]] * r$conseq)
  c(unweighted = mean(p), weighted = sum(p * r$wt) / sum(r$wt)) }
print(round(sapply(c(solar = "solar", wind = "wind"), function(s) shr(lc[[s]]$lDA, s)), 3))
cat("(raw share of respondents who ALWAYS chose the status quo: solar 0.359, wind 0.520)\n")

cat("\nMembership of the status-quo-only class (LC2-DA): odds ratios, robust z in brackets\n")
gn <- c(set_gd, set_ga)
lab <- c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "exposed", "env. organisation",
         "attitude Q7 (+1)", "concern Q10 (+1)", "visit nature Q31 (+1)", "preserve nature Q32 (+1)", "believes results used")
fo <- function(m) { e <- m$estimate; se <- rse(m); sprintf("%5.2f (%5.1f)", exp(e[gn]), e[gn] / se[gn]) }
print(data.frame(variable = lab, solar = fo(lc$solar$lDA), wind = fo(lc$wind$lDA)), row.names = FALSE)

cat("\nWTA of the TRADER class, NOK per household-year (robust s.e.); LC2-DA\n")
wt <- function(m) { e <- m$estimate; se <- rse(m); k <- setdiff(wn, "b_cost"); sg <- ifelse(k == "w_sq", 100, -100)
  data.frame(attr = sub("^w_", "", k), WTA = round(sg * e[k]), se = round(100 * se[k])) }
a <- wt(lc$solar$lDA); b <- wt(lc$wind$lDA)
print(data.frame(attribute = a$attr, solar = a$WTA, se_s = a$se, wind = b$WTA, se_w = b$se), row.names = FALSE)
cat(sprintf("\nStatus-quo-only class: P(status quo) = %.3f (solar), %.3f (wind)\n", plogis(lc$solar$lDA$estimate[["a_sq2"]]), plogis(lc$wind$lDA$estimate[["a_sq2"]])))

for (s in c("solar", "wind")) if (!inherits(lc3[[s]], "try-error")) {
  m3 <- lc3[[s]]; e <- m3$estimate; r <- dat[[s]][!duplicated(dat[[s]]$id), ]
  v2 <- e[["d2"]] + e[["h2_age"]] * r$age_c + e[["h2_fem"]] * r$female + e[["h2_att"]] * r$att_c + e[["h2_conc"]] * r$conc_c + e[["h2_natp"]] * r$natp_c + e[["h2_cons"]] * r$conseq + e[["h2_icost"]] * r$imp_cost_c
  v3 <- e[["d3"]] + e[["h3_age"]] * r$age_c + e[["h3_fem"]] * r$female + e[["h3_att"]] * r$att_c + e[["h3_conc"]] * r$conc_c + e[["h3_natp"]] * r$natp_c + e[["h3_cons"]] * r$conseq + e[["h3_icost"]] * r$imp_cost_c
  den <- 1 + exp(v2) + exp(v3)
  cat(sprintf("\n[%s LC3] class shares: traders %.2f, cost-ignorers %.2f, status-quo-only %.2f;  P(status quo | cost-ignorers, reference scenario) = %.2f\n",
              s, mean(1 / den), mean(exp(v2) / den), mean(exp(v3) / den), plogis(e[["a_sq_c"]])))
  cat("  cost-ignorers, utility weights: ", paste(sub("^e_", "", names(e)[grepl("^e_", names(e))]), round(e[grepl("^e_", names(e))], 2), collapse = "; "), "\n")
  cat("  membership of the cost-ignorer class, odds ratio (z): stated cost importance (Q25) ", sprintf("%.2f (%.1f)", exp(e[["h2_icost"]]), e[["h2_icost"]] / rse(m3)[["h2_icost"]]),
      "; attitude ", sprintf("%.2f (%.1f)", exp(e[["h2_att"]]), e[["h2_att"]] / rse(m3)[["h2_att"]]), "\n")
}
