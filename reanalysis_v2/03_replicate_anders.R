# =====================================================================
# 03  Independent replication of Anders's mixed logit (his Table A.3), with my own code
# =====================================================================
# Specification (from the article text and his replication script):
#   WTA space, panel; 11 Normal coefficients; correlation only WITHIN an attribute (Cholesky blocks);
#   lognormal cost coefficient that multiplies the whole utility.   30 parameters.
# Starting values are the published estimates, so this checks whether the published optimum is reproduced
# (different code, 500 MLHS draws instead of 2,000).
A3 <- list(
  wind = data.frame(par = c("sq", "el38", "el57", "pl3", "pl6", "nop", "wil", "d1", "d2", "d3", "view", "cost"),
    mean = c(137.256, 2.863, 9.602, -7.261, -7.949, -52.835, -12.876, -95.514, -66.840, -25.801, -28.000, -1.860),
    mean_se = c(0.703, 0.271, 0.334, 0.254, 0.293, 0.248, 0.234, 0.364, 0.351, 0.251, 0.210, 0.316),
    sd = c(261.896, 23.445, 25.249, 11.505, 39.353, 71.033, 3.703, 63.241, 7.169, 13.783, 32.003, 2.696),
    sd_se = c(1.437, 0.155, 0.299, 0.141, 0.242, 0.376, 0.134, 0.442, 0.189, 0.107, 0.226, 0.334), LL = -2065.2014, BIC = 4392, k = 30),
  solar = data.frame(par = c("sq", "el38", "el57", "pl3", "pl6", "nop", "wil", "d1", "d2", "d3", "view", "cost"),
    mean = c(38.722, -9.485, -13.704, -7.853, -17.439, -74.290, -17.546, -43.593, -25.899, -5.919, -46.991, -2.612),
    mean_se = c(0.449, 1.039, 0.680, 0.625, 0.803, 1.601, 0.732, 1.732, 3.193, 1.115, 0.685, 0.262),
    sd = c(206.727, 3.946, 3.806, 29.425, 1.854, 80.416, 50.478, 19.739, 8.780, 7.753, 50.374, 2.178),
    sd_se = c(5.690, 0.853, 0.345, 0.718, 0.364, 1.188, 1.652, 0.284, 0.502, 0.304, 1.551, 0.319), LL = -2536.36, BIC = 5335, k = 30))

prob_A <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  V <- list()
  V[["dev"]] <- b_cost * (b_el38 * (el_alt1 == 38) + b_el57 * (el_alt1 == 57) + b_pl3 * (plants_alt1 == 3) + b_pl6 * (plants_alt1 == 6) +
                          b_nop * (wild_alt1 == 0) + b_wil * (wild_alt1 == 1) + b_d1 * (dist_alt1 == 1) + b_d2 * (dist_alt1 == 2) +
                          b_d3 * (dist_alt1 == 3) + b_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- b_cost * b_sq
  P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
  P <- apollo_panelProd(P, apollo_inputs, functionality); P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
rc_A <- function(apollo_beta, apollo_inputs) {
  r <- list()
  r[["b_sq"]]   <- m_sq   + s_sq   * xi_sq
  r[["b_el38"]] <- m_el38 + s_el38 * xi_el38
  r[["b_el57"]] <- m_el57 + s_el57 * xi_el57 + c_el * xi_el38
  r[["b_pl3"]]  <- m_pl3  + s_pl3  * xi_pl3
  r[["b_pl6"]]  <- m_pl6  + s_pl6  * xi_pl6  + c_pl * xi_pl3
  r[["b_nop"]]  <- m_nop  + s_nop  * xi_nop
  r[["b_wil"]]  <- m_wil  + s_wil  * xi_wil  + c_wil * xi_nop
  r[["b_d1"]]   <- m_d1   + s_d1   * xi_d1
  r[["b_d2"]]   <- m_d2   + s_d2   * xi_d2   + c_d12 * xi_d1
  r[["b_d3"]]   <- m_d3   + s_d3   * xi_d3   + c_d13 * xi_d1 + c_d23 * xi_d2
  r[["b_view"]] <- m_view + s_view * xi_view
  r[["b_cost"]] <- exp(m_cost + s_cost * xi_cost)
  return(r)
}
dr_A <- list(interDrawsType = "mlhs", interNDraws = 500,
             interNormDraws = paste0("xi_", c("sq", "el38", "el57", "pl3", "pl6", "nop", "wil", "d1", "d2", "d3", "view", "cost")))

a <- A3[[SAMPLE]]
pn <- c("sq", "el38", "el57", "pl3", "pl6", "nop", "wil", "d1", "d2", "d3", "view")
st <- c(setNames(a$mean[1:11], paste0("m_", pn)), m_cost = a$mean[12],
        setNames(a$sd[1:11], paste0("s_", pn)), s_cost = a$sd[12],
        c_el = 0, c_pl = 0, c_wil = 0, c_d12 = 0, c_d13 = 0, c_d23 = 0)
t0 <- Sys.time()
mA <- fit_model(paste0("repA_", SAMPLE), dat[[SAMPLE]], prob_A, st, draws = dr_A, randcoeff = rc_A,
                estimate_settings = list(estimationRoutine = "bfgs", maxIterations = 300))
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

cat(sprintf("Replication of Anders's model, %s sample (%.0f min, converged: %s)\n", SAMPLE, mins, mA$message))
cat(sprintf("\n%-28s %10s %10s %6s\n", "", "published", "mine", "diff"))
cat(sprintf("%-28s %10.1f %10.1f %6.1f\n", "log-likelihood", a$LL[1], mA$maximum, mA$maximum - a$LL[1]))
cat(sprintf("%-28s %10.0f %10.1f %6.1f\n", "BIC", a$BIC[1], bic_of(mA), bic_of(mA) - a$BIC[1]))
cat(sprintf("%-28s %10d %10d\n", "parameters", a$k[1], nfree(mA)))

e <- mA$estimate; sr <- rse(mA); sc <- cse(mA)
sdv <- c(sq = e[["s_sq"]], el38 = e[["s_el38"]], el57 = sqrt(e[["s_el57"]]^2 + e[["c_el"]]^2), pl3 = e[["s_pl3"]], pl6 = sqrt(e[["s_pl6"]]^2 + e[["c_pl"]]^2),
         nop = e[["s_nop"]], wil = sqrt(e[["s_wil"]]^2 + e[["c_wil"]]^2), d1 = e[["s_d1"]], d2 = sqrt(e[["s_d2"]]^2 + e[["c_d12"]]^2),
         d3 = sqrt(e[["s_d3"]]^2 + e[["c_d13"]]^2 + e[["c_d23"]]^2), view = e[["s_view"]])
tab <- data.frame(parameter = pn, mean_published = a$mean[1:11], se_published = a$mean_se[1:11],
                  mean_mine = round(e[paste0("m_", pn)], 2), se_classical = round(sc[paste0("m_", pn)], 2), se_robust = round(sr[paste0("m_", pn)], 2),
                  sd_published = a$sd[1:11], sd_mine_chol = round(abs(e[paste0("s_", pn)]), 1), row.names = NULL)
cat("\nMeans in 100 NOK (a coefficient of -52.8 = 5,284 NOK of compensation); sd_mine_chol = diagonal Cholesky term\n")
print(tab, row.names = FALSE)
cat(sprintf("\nCost: mean of log b_cost  published %.2f  mine %.2f (robust s.e. %.2f);  sd of log b_cost  published %.2f  mine %.2f\n",
            a$mean[12], e[["m_cost"]], sr[["m_cost"]], a$sd[12], abs(e[["s_cost"]])))
cat("\nRatio of published standard error to my robust standard error for the coefficient MEANS (published / mine):\n")
print(round(a$mean_se[1:11] / sr[paste0("m_", pn)], 3))
