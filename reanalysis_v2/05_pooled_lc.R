# =====================================================================
# 05  Is wind different from solar?  A formal test inside the latent class model
# =====================================================================
# Article, conclusion: "this result is based on separate models ... further analysis can include testing whether the differences in
# preferences between the technologies are statistically significant."   Leeds [W] 4B (joint estimation, scale) + [A] 2A (latent class).
# Pooled two-class model; the TRADER class has preferences w (WTA space) and the status-quo-only class chooses the status quo with fixed probability.
#   A  traders have the SAME attribute WTAs and the same status quo constant in both technologies (wind may have its own scale, lambda);
#      the class shares depend on the same covariates, with a wind intercept only
#   B  = A + technology-specific status quo constants (traders and status-quo-only class)
#   C  = B + technology-specific effects of ALL covariates on the class shares
#   U  = the two separate LC2-DA models of step 06 (all trader WTAs technology-specific)
# Tests:  A -> B   constants;   B -> C   who ends up in the status-quo-only class;   C -> U   do the ten attribute WTAs differ?
e06 <- new.env(); load(file.path(V2, "state_06_both.RData"), envir = e06); lc <- e06$lc
zc12 <- c("age_c", "female", "univ", "inc_c", "inc_missing", "exposed", "org", "att_c", "conc_c", "natv_c", "natp_c", "conseq")
lcp_pool <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  base <- d_cls + g_wind * wind + gd_age * age_c + gd_fem * female + gd_univ * univ + gd_inc * inc_c + gd_incm * inc_missing + gd_exp * exposed + gd_org * org +
          ga_att * att_c + ga_conc * conc_c + ga_natv * natv_c + ga_natp * natp_c + gc_cons * conseq
  inter <- wd_age * wind * age_c + wd_fem * wind * female + wd_univ * wind * univ + wd_inc * wind * inc_c + wd_incm * wind * inc_missing + wd_exp * wind * exposed + wd_org * wind * org +
           wa_att * wind * att_c + wa_conc * wind * conc_c + wa_natv * wind * natv_c + wa_natp * wind * natp_c + wc_cons * wind * conseq
  lcpars[["pi_values"]] <- apollo_classAlloc(list(classes = c(trader = 1, sqonly = 2), utilities = list(trader = 0, sqonly = base + inter)))
  return(lcpars)
}
prob_pool <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  sc <- 1 + (lambda - 1) * wind                                  # solar scale = 1, wind scale = lambda
  P <- list(); mnl <- list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice)
  V <- list()
  V[["dev"]] <- sc * b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
                               w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
                               w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- sc * b_cost * (w_sq + dsq_w * wind)
  mnl$utilities <- V; mnl$componentName <- "trader"
  P[["trader"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
  mnl$utilities <- list(dev = 0, sq = a_sq2 + da_w * wind); mnl$componentName <- "sqonly"
  P[["sqonly"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
  P[["model"]] <- apollo_lc(list(inClassProb = P, classProb = pi_values), apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view", "b_cost")
gn <- c("gd_age", "gd_fem", "gd_univ", "gd_inc", "gd_incm", "gd_exp", "gd_org", "ga_att", "ga_conc", "ga_natv", "ga_natp", "gc_cons")
wi <- c("wd_age", "wd_fem", "wd_univ", "wd_inc", "wd_incm", "wd_exp", "wd_org", "wa_att", "wa_conc", "wa_natv", "wa_natp", "wc_cons")
es <- lc$solar$lDA$estimate
st <- c(es[c(wn, "a_sq2", "d_cls", gn)], lambda = 1, dsq_w = 0, da_w = 0, g_wind = 0.3, setNames(rep(0, 12), wi))
fA <- c("dsq_w", "da_w", wi); fB <- wi
mA <- fit_model("pool_A", dat$pooled, prob_pool, st, fixed = fA, lcpars = lcp_pool)
mB <- fit_model("pool_B", dat$pooled, prob_pool, mA$estimate[names(st)] |> replace(c("dsq_w", "da_w"), c(0, 0)), fixed = fB, lcpars = lcp_pool)
mC <- fit_model("pool_C", dat$pooled, prob_pool, mB$estimate[names(st)], lcpars = lcp_pool)
LL_U <- lc$solar$lDA$maximum + lc$wind$lDA$maximum; k_U <- nfree(lc$solar$lDA) + nfree(lc$wind$lDA)
n_all <- lc$solar$lDA$nObs + lc$wind$lDA$nObs
ms <- list(mA, mB, mC)
tb <- data.frame(model = c("A  common WTA, common constants", "B  + technology-specific constants", "C  + technology-specific class-share effects", "U  all separate (two LC2-DA models)"),
                 LL = round(c(sapply(ms, function(m) m$maximum), LL_U), 1), k = c(sapply(ms, nfree), k_U))
tb$BIC <- round(-2 * tb$LL + tb$k * log(n_all), 1); print(tb, row.names = FALSE)
lr <- function(l1, k1, l0, k0) { x <- 2 * (l1 - l0); data.frame(LR = round(x, 1), df = k1 - k0, p = signif(pchisq(x, k1 - k0, lower.tail = FALSE), 2)) }
cat("\nLikelihood-ratio tests\n")
print(rbind("A -> B  constants differ between technologies?" = lr(mB$maximum, nfree(mB), mA$maximum, nfree(mA)),
            "B -> C  covariates act differently on class shares?" = lr(mC$maximum, nfree(mC), mB$maximum, nfree(mB)),
            "C -> U  do the 10 attribute WTAs of the traders differ?" = lr(LL_U, k_U, mC$maximum, nfree(mC)),
            "A -> U  everything" = lr(LL_U, k_U, mA$maximum, nfree(mA))))
e <- mB$estimate; se <- rse(mB)
cat(sprintf("\nModel B: wind scale lambda = %.2f (robust s.e. %.2f); wind intercept in the class-share equation %.2f (z %.1f)\n", e[["lambda"]], se[["lambda"]], e[["g_wind"]], e[["g_wind"]] / se[["g_wind"]]))
cat(sprintf("Status quo constant of traders: solar %.0f NOK, wind %.0f NOK (difference %.0f, z %.1f);  P(status quo | status-quo-only class): solar %.3f, wind %.3f\n",
            100 * e[["w_sq"]], 100 * (e[["w_sq"]] + e[["dsq_w"]]), 100 * e[["dsq_w"]], e[["dsq_w"]] / se[["dsq_w"]], plogis(e[["a_sq2"]]), plogis(e[["a_sq2"]] + e[["da_w"]])))
