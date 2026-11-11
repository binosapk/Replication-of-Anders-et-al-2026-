################################################################################
#                                                                              #
#   SOLAR versus WIND POWER                                                    #
#   Re-analysis of Dugstad, Hoel and Lindhjem (2026), "Kraft med mindre        #
#   konflikt? Preferanser for sol- vs. vindkraft", Samfunnsokonomen 138(2)     #
#                                                                              #
#   Contributor: Binod                                                         #
#   Data: sol_complete.csv, vind_complete.csv, CE_design_final.csv (original   #
#   study), ssb_pop.csv (population margins, Statistics Norway)                #
#                                                                              #
################################################################################
#
# HOW TO RUN (from the folder that holds this file; nothing to edit):
#
#     Rscript Solar_vs_Wind_Reanalysis.R <step> <sample>
#
# For example   Rscript Solar_vs_Wind_Reanalysis.R 00 all      (builds the data; run this first)
#               Rscript Solar_vs_Wind_Reanalysis.R 02 all
#               Rscript Solar_vs_Wind_Reanalysis.R 04 solar
#               Rscript Solar_vs_Wind_Reanalysis.R 13 all
# One run executes ONE step.  The output goes to the console (redirect it to a file).  Each step saves its results in
# state_<step>_<sample>.RData, which later steps load.  Some steps are slow (hours; the pooled mixed logit takes 1 to 3 days per job).
#
# STEPS
#   00  Data: design, covariates, raking weights                      sample: any
#   01  Descriptive statistics                                        all
#   02  MNL ladder (shifters, income, attitudes)                      all
#   03  Replication of the published mixed logit                      solar | wind
#   04  Full-covariance mixed logit                                   solar | wind
#   05  Pooled two-class latent class                                 both
#   06  Latent classes LC2 and LC3                                    both
#   06b LC3 with the trader cost coefficient fixed                    both
#   07  Hybrid choice model                                           solar | wind
#   08  Bayesian full-covariance model (not converged)                solar | wind
#   09  Cost per kWh and acceptance (policy)                          all
#   10  Exclusion rules (robustness)                                  all
#   11  Pooled mixed logit, all restrictions                          job: U2 | scale2 | cov2 | attr2 | sq2 | all2 (first round: U | scale | cov | attr | sq | g | all)
#   12  Pooled MNL scale test                                         all
#   13  Pooled mixed logit: tests (LR, Holm) and WTA differences      all
#   14  Figure: acceptance curves                                     all
#   15  Pooled mixed logit with a scaled optimiser                    job: U3 | scale3 | cov3 | attr3 | sq3 | all3
#   16  Collinearity check of the covariates                          all
#   17  Pooled mixed logit: draw-seed check (another seed)            job: U4 | all4
#
# QUICK TEST: set the environment variable SW_QUICK=1 to run any step with 10 draws and 2 iterations. It writes only quick_* files and checks
# that the code runs; its numbers mean nothing.
#
# Units: 100 NOK in the estimates; reported WTA is in NOK per household and year.
#
################################################################################


# ==============================================================================
# 0  SETTINGS, FOLDERS, PACKAGES AND HELPER FUNCTIONS
# ==============================================================================

# ---- folder of this script (V2) and its parent (ROOT, with the data files)
.f   <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
V2   <- normalizePath(if (length(.f)) dirname(.f[1]) else getwd(), winslash = "/")
ROOT <- dirname(V2)

# ---- step and sample (from the command line)
args    <- commandArgs(trailingOnly = TRUE)
.step   <- if (length(args) >= 1) args[1] else "00"
.sample <- if (length(args) >= 2) args[2] else "all"

# ---- results of earlier steps (state_base.RData holds the data built by step 00)
if (.step != "00") {
  .e <- new.env()
  load(file.path(V2, "state_base.RData"), envir = .e)
  for (.n in ls(.e, all.names = TRUE)) assign(.n, get(.n, .e), globalenv())
  rm(.e, .n)
}

# ---- settings (set after loading, so that older saved files cannot overwrite them)
STEP   <- .step
SAMPLE <- .sample
QUICK  <- nzchar(Sys.getenv("SW_QUICK"))
QK     <- function(full, quick) if (QUICK) quick else full
V2     <- normalizePath(if (length(.f)) dirname(.f[1]) else getwd(), winslash = "/")
ROOT   <- dirname(V2)
setwd(ROOT)

outdir <- file.path(V2, if (QUICK) "apollo_output_quick" else "apollo_output")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(apollo)
  library(dplyr)
  library(tidyr)
})

# ---- helper functions
quiet <- function(expr) {
  out <- NULL
  invisible(capture.output(out <- suppressMessages(expr)))
  out
}
num1   <- function(x) suppressWarnings(as.numeric(substr(x, 1, 1)))        # "4 Neutral" -> 4
nfree  <- function(m) length(m$estimate) - length(m$apollo_fixed)
bic_of <- function(m) -2 * m$maximum + nfree(m) * log(m$nObs)
aic_of <- function(m) -2 * m$maximum + 2 * nfree(m)
rse    <- function(m) {
  v <- sqrt(diag(m$robvarcov))
  names(v) <- rownames(m$robvarcov)
  v
}   # robust s.e. with names
cse    <- function(m) {
  v <- sqrt(diag(m$varcov))
  names(v) <- rownames(m$varcov)
  v
}   # classical s.e.
lrp    <- function(big, small) signif(pchisq(2 * (big$maximum - small$maximum), nfree(big) - nfree(small), lower.tail = FALSE), 2)

# one wrapper that estimates any Apollo model silently
fit_model <- function(name, data, probs, beta, fixed = c(), draws = NULL, randcoeff = NULL, lcpars = NULL,
                      ncores = 1, estimate_settings = list(), seed = NULL) {
  database <<- data
  quiet(apollo_initialise())
  apollo_control <<- list(modelName = name, indivID = "id", mixing = !is.null(draws), nCores = ncores,
    noDiagnostics = TRUE, outputDirectory = outdir)
  if (!is.null(seed)) apollo_control$seed <<- seed
  apollo_beta <<- beta
  apollo_fixed <<- fixed
  if (!is.null(draws)) apollo_draws     <<- draws
  if (!is.null(randcoeff)) apollo_randCoeff <<- randcoeff
  if (!is.null(lcpars)) apollo_lcPars    <<- lcpars
  on.exit(suppressWarnings(rm(apollo_draws, apollo_randCoeff, apollo_lcPars, envir = globalenv())))
  apollo_inputs <<- apollo_validateInputs(silent = TRUE)
  es <- modifyList(list(silent = TRUE, writeIter = FALSE), estimate_settings)
  if (QUICK) es$maxIterations <- 2
  quiet(apollo_estimate(beta, fixed, probs, apollo_inputs, estimate_settings = es))
}

# run the code of a step only if it is the requested one; top-level values are printed as when the code is sourced
section <- function(steps, code) {
  if (!(STEP %in% steps)) return(invisible(NULL))
  ex <- substitute(code)
  if (is.symbol(ex)) ex <- eval(ex, globalenv())
  withAutoprint(ex, evaluated = TRUE, local = FALSE, echo = FALSE)
  invisible(NULL)
}


# ==============================================================================
# PUBLISHED RESULTS: Table A.3 of the article (units: 100 NOK), transcribed from the article image
# ==============================================================================
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



.t0 <- Sys.time()
.before <- ls(all.names = TRUE)

################################################################################
# 00  Data preparation: choice data, respondent covariates, population weights
################################################################################
section("00", {
  set.seed(2026)

  # ---------------------------------------------------------------- design (one row per card)
  des <- read.csv("CE_design_final.csv")
  des <- des[des$alt == 1, c("choice_set", "el", "plants", "wild", "dist", "view", "cost")]
  names(des)[-1] <- paste0(names(des)[-1], "_alt1")

  # ---------------------------------------------------------------- income: category text -> midpoint (million NOK)
  income_mid <- function(x) {
    y <- gsub("[[:space:]\u00a0]", "", x)
    sapply(y, function(s) {
      if (is.na(s) || s == "" || grepl("nsker|Vetikke", s)) return(NA_real_)
      nums <- as.numeric(regmatches(s, gregexpr("[0-9]+", s))[[1]])
      if (grepl("^Inntil", s)) return(nums[1] / 2 / 1e6)
      if (grepl("^Mer", s)) return(4.5)
      if (length(nums) >= 2) return(mean(nums[1:2]) / 1e6)
      NA_real_
    })
  }

  # ---------------------------------------------------------------- one survey file -> one row per answered card
  prepare <- function(file, tech, id_start) {
    r <- read.csv(file, encoding = "UTF-8", stringsAsFactors = FALSE)
    r$id <- id_start + seq_len(nrow(r))
    r$tech <- tech
    r$wind <- as.numeric(tech == "wind")
    # --- demographics (exogenous)
    r$female <- as.numeric(r$gender == "Kvinne")
    r$age_c  <- (r$age - 50) / 10
    r$agegrp <- cut(r$age, c(17, 29, 39, 49, 59, 120), labels = c("18-29", "30-39", "40-49", "50-59", "60+"))
    r$educ4  <- ifelse(grepl("^Grunnskole", r$Q38), 1,
      ifelse(grepl("^Videreg|^Fagbrev", r$Q38), 2,
        ifelse(grepl("^3-4", r$Q38), 3,
          ifelse(grepl("^5-|^PhD", r$Q38), 4, NA))))
    r$educ_missing <- as.numeric(is.na(r$educ4))
    r$educ4[is.na(r$educ4)] <- 2     # 1% "unsure": modal category
    r$univ   <- as.numeric(r$educ4 >= 3)
    inc <- income_mid(r$Q40)
    r$inc_missing <- as.numeric(is.na(inc))
    r$inc_c  <- ifelse(is.na(inc), 0, inc - 0.8)                                    # household income, million NOK, centred at 0.8
    r$org    <- as.numeric(r$Q37 == "Ja")
    r$near   <- as.numeric(grepl("bor i n", r$Q39))
    r$owner <- as.numeric(grepl("eier arealer", r$Q39))
    r$exposed <- as.numeric(r$near == 1 | r$owner == 1)                             # lives near / owns land
    r$sees   <- as.numeric(r$Q13 == "Ja")
    r$own_pv <- as.numeric(r$Q36 == "Ja")
    # --- attitudes and values (potentially endogenous; used as indicators or as diagnostics)
    r$att    <- num1(r$Q7)
    r$conc <- num1(r$Q10)
    r$nat_visit <- num1(r$Q31)
    r$nat_pres <- ifelse(grepl("ikke redusere", r$Q32), 3, ifelse(grepl("noe reduksjon", r$Q32), 2, ifelse(grepl("^Nei", r$Q32), 1, NA)))
    r$impact_own <- if (tech == "wind") num1(r$Q33r1) else num1(r$Q33r2)           # perceived landscape impact of own technology
    r$imp_el <- num1(r$Q25r1)
    r$imp_pl <- num1(r$Q25r2)
    r$imp_wild <- num1(r$Q25r3)
    r$imp_dist <- num1(r$Q25r4)
    r$imp_view <- num1(r$Q25r5)
    r$imp_cost <- num1(r$Q25r6)
    r$difficult <- num1(r$Q26)
    r$conseq <- as.numeric(r$Q41 %in% c("Veldig sikker", "Ganske sikker"))           # believes results will be used
    r$consq_missing <- 0
    keep <- c("id", "tech", "wind", "female", "age", "age_c", "agegrp", "educ4", "educ_missing", "univ", "inc_c", "inc_missing", "org", "near", "owner",
      "exposed", "sees", "own_pv", "att", "conc", "nat_visit", "nat_pres", "impact_own", "imp_el", "imp_pl", "imp_wild", "imp_dist", "imp_view",
      "imp_cost", "difficult", "conseq", "qtime", "Q27", "Q27b")
    r[, c(keep, grep("^Q17_", names(r), value = TRUE))] %>%
      pivot_longer(starts_with("Q17_"), names_to = "choice_set", names_prefix = "Q17_", values_to = "answer") %>%
      mutate(choice_set = as.integer(choice_set), choice = case_when(answer == "Utbygging" ~ 1L, answer == "Ingen utbygging" ~ 2L)) %>%
      filter(!is.na(choice)) %>% select(-answer) %>% left_join(des, by = "choice_set") %>% mutate(cost_alt2 = 0) %>% arrange(id, choice_set)
  }
  dat <- list(solar = prepare("sol_complete.csv", "solar", 0), wind = prepare("vind_complete.csv", "wind", 10000))

  # ---------------------------------------------------------------- population weights (raking to Statistics Norway, Innlandet 2025)
  pop <- read.csv(file.path(V2, "ssb_pop.csv"), stringsAsFactors = FALSE)
  names(pop) <- c("region", "sex", "age_txt", "n")
  pop$age <- as.numeric(sub("^([0-9]+).*", "\\1", pop$age_txt))
  pop <- pop[pop$age >= 18, ]
  pop$agegrp <- cut(pop$age, c(17, 29, 39, 49, 59, 120), labels = c("18-29", "30-39", "40-49", "50-59", "60+"))
  pop$female <- as.numeric(pop$sex == "Females")
  pop_cell <- aggregate(n ~ agegrp + female, pop, sum)
  pop_cell$share <- pop_cell$n / sum(pop_cell$n)
  edu_pop <- c(`1` = 83586, `2` = 127937 + 11217, `3` = 72522, `4` = 24563)          # basic; upper secondary (+ tertiary vocational); higher short; higher long
  edu_share <- edu_pop / sum(edu_pop)

  rake <- function(r, iter = 100) {            # r: one row per respondent
    w <- rep(1, nrow(r))
    cell <- paste(r$agegrp, r$female)
    pcell <- setNames(pop_cell$share, paste(pop_cell$agegrp, pop_cell$female))
    for (k in 1:iter) {
      s <- tapply(w, cell, sum) / sum(w)
      w <- w * (pcell[cell] / s[cell])
      s <- tapply(w, r$educ4, sum) / sum(w)
      w <- w * (edu_share[as.character(r$educ4)] / s[as.character(r$educ4)])
    }
    w / mean(w)
  }
  for (s in c("solar", "wind")) {
    d <- dat[[s]]
    r <- d[!duplicated(d$id), ]
    w <- setNames(rake(r), r$id)
    d$wt <- unname(w[as.character(d$id)])
    dat[[s]] <- d
  }

  # ---------------------------------------------------------------- derived variables shared by all models
  dat <- lapply(dat, function(d) {
    d %>% group_by(id) %>% mutate(always_sq = as.numeric(all(choice == 2)), task = row_number()) %>% ungroup() %>%
      mutate(att_c = ifelse(is.na(att), 4, att) - 4, conc_c = ifelse(is.na(conc), 4, conc) - 4,
        natv_c = ifelse(is.na(nat_visit), 4, nat_visit) - 4, natp_c = ifelse(is.na(nat_pres), 2, nat_pres) - 2,
        imp_cost_c = ifelse(is.na(imp_cost), 4, imp_cost) - 4) %>% as.data.frame()
  })
  dat$pooled <- rbind(dat$solar, dat$wind)
  rm(s, d, r, w)

  # ---------------------------------------------------------------- report: sample vs population
  cat("Sample vs population (Statistics Norway, Innlandet 2025, age 18+)\n")
  cmp <- function(s) {
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    a <- prop.table(table(r$agegrp))
    e <- prop.table(table(r$educ4))
    f <- mean(r$female)
    data.frame(variable = c(paste("age", levels(r$agegrp)), "women", "basic school", "upper secondary", "higher short (<=4y)", "higher long (5y+)"),
      sample = round(100 * c(a, f, e), 1))
  }
  pa <- c(tapply(pop$n, pop$agegrp, sum) / sum(pop$n), women = sum(pop$n[pop$female == 1]) / sum(pop$n), edu_share)
  tab <- cbind(cmp("solar"), wind = cmp("wind")$sample, population = round(100 * pa, 1))
  names(tab)[2] <- "solar"
  print(tab, row.names = FALSE)
  cat("\nWeights after raking on age x sex and education: summary and effective sample size\n")
  for (s in c("solar", "wind")) {
    w <- dat[[s]]$wt[!duplicated(dat[[s]]$id)]
    cat(sprintf("%-6s min %.2f  max %.2f  effective n = %.0f of %d\n", s, min(w), max(w), sum(w)^2 / sum(w^2), length(w)))
  }
  cat("\nShare of respondents with missing information: income", round(mean(dat$pooled$inc_missing[!duplicated(dat$pooled$id)]), 3),
    " education", round(mean(dat$pooled$educ_missing[!duplicated(dat$pooled$id)]), 3), "\n")
})


################################################################################
# 01  Study the data first (Leeds [A] 1B: "start with the data, even before a simple MNL, analyse them non-parametrically")
################################################################################
section("01", {
  # (a) acceptance by compensation level, with confidence intervals and population weights
  # (b) who never trades? binary logit of "chose the status quo on all 8 cards"
  # (c) does belief that the results will be used (consequentiality, Q41) matter?
  wmean <- function(x, w) sum(x * w) / sum(w)
  cat("(a) Share choosing the development option by bill reduction (%), unweighted [95% CI] and weighted\n")
  acc <- do.call(rbind, lapply(c("solar", "wind"), function(s) {
    d <- dat[[s]]
    do.call(rbind, lapply(sort(unique(d$cost_alt1)), function(cst) {
      z <- d[d$cost_alt1 == cst, ]
      p <- mean(z$choice == 1)
      n <- nrow(z)
      # cluster-robust s.e. (respondents answer several cards): use respondent-level means
      zi <- tapply(z$choice == 1, z$id, mean)
      se <- sd(zi) / sqrt(length(zi))
      data.frame(sample = s, cost = cst, n = n, share = round(100 * p, 1), lo = round(100 * (p - 1.96 * se), 1), hi = round(100 * (p + 1.96 * se), 1),
        weighted = round(100 * wmean(z$choice == 1, z$wt), 1))
    }))
  }))
  print(reshape(acc[, c("sample", "cost", "share", "lo", "hi", "weighted")], idvar = "cost", timevar = "sample", direction = "wide"), row.names = FALSE)

  # slope test: logit of development choice on cost, cluster-robust by respondent (cards differ in other attributes, so this is only descriptive)
  cat("\nDescriptive slope: change in P(development) per +1,000 NOK (linear probability model, respondent-clustered s.e.)\n")
  for (s in c("solar", "wind")) {
    d <- dat[[s]]
    d$y <- as.numeric(d$choice == 1)
    d$c1000 <- d$cost_alt1 / 1000
    f <- lm(y ~ c1000, d)
    X <- model.matrix(f)
    e <- resid(f)
    meat <- crossprod(rowsum(X * e, d$id))
    bread <- solve(crossprod(X))
    V <- bread %*% meat %*% bread
    cat(sprintf("%-6s slope %.4f (cluster s.e. %.4f), t = %.2f\n", s, coef(f)[2], sqrt(V[2, 2]), coef(f)[2] / sqrt(V[2, 2])))
  }

  cat("\n(b) Who NEVER trades?  Binary logit, odds ratios (z in brackets), respondent level\n")
  resp <- lapply(c("solar", "wind"), function(s) {
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    r
  })
  names(resp) <- c("solar", "wind")
  form <- always_sq ~ I(age_c) + female + univ + inc_c + inc_missing + org + exposed + att_c + conc_c + natv_c + natp_c + conseq
  lab <- c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "env. organisation", "exposed (lives near / owns land)",
    "attitude Q7 (+1)", "concern Q10 (+1)", "visit nature Q31 (+1)", "preserve nature Q32 (+1)", "believes results are used")
  fits <- lapply(resp, function(r) glm(form, binomial, r))
  tab <- sapply(fits, function(f) {
    s <- summary(f)$coefficients[-1, ]
    sprintf("%5.2f (%5.1f)", exp(s[, 1]), s[, 3])
  })
  print(data.frame(variable = lab, solar = tab[, 1], wind = tab[, 2]), row.names = FALSE)
  cat(sprintf("McFadden pseudo-R2: solar %.3f, wind %.3f\n", 1 - fits$solar$deviance / fits$solar$null.deviance, 1 - fits$wind$deviance / fits$wind$null.deviance))

  cat("\n(c) Consequentiality (Q41): 'very/fairly sure the authorities will use the results'\n")
  print(do.call(rbind, lapply(c("solar", "wind"), function(s) {
    r <- resp[[s]]
    data.frame(sample = s, share_sure = round(mean(r$conseq), 3),
      always_SQ_if_sure = round(mean(r$always_sq[r$conseq == 1]), 3), always_SQ_if_unsure = round(mean(r$always_sq[r$conseq == 0]), 3))
  })), row.names = FALSE)

  cat("\nWeighted vs unweighted share of respondents who ALWAYS chose the status quo\n")
  print(round(sapply(resp, function(r) c(unweighted = mean(r$always_sq), weighted = wmean(r$always_sq, r$wt))), 3))
})


################################################################################
# 02  Deterministic heterogeneity ladder (MNL in WTA space)  -  exogenous covariates first, attitudes second
################################################################################
section("02", {
  # V_dev = bc * ( sum_k w_k x_k + cost/100 ),   V_sq = bc * ( w_sq + shifters ),   bc = b_cost * exp(g_cinc * income)
  # Rung 0  baseline
  # Rung 1  status quo depends on EXOGENOUS characteristics: age, gender, university degree, household income, exposure, env. organisation
  # Rung 2  + cost sensitivity depends on income (a priori the most natural cost interaction)
  # Rung 3  + ATTITUDES: general attitude (Q7), concern (Q10), visits nature (Q31), preservation (Q32)      <- endogenous, diagnostic
  # Rung 4  + believes the results will be used (Q41)
  prob_ladder <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    sq <- w_sq + g_age * age_c + g_fem * female + g_univ * univ + g_inc * inc_c + g_incm * inc_missing + g_exp * exposed + g_org * org +
      g_att * att_c + g_conc * conc_c + g_natv * natv_c + g_natp * natp_c + g_cons * conseq
    bc <- b_cost * exp(g_cinc * inc_c + g_cincm * inc_missing)
    V <- list()
    V[["dev"]] <- bc * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
      w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
      w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- bc * sq
    P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
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
    M <- lad[[s]]
    nm <- c("0 baseline", "1 + exogenous shifters on status quo", "2 + cost x income", "3 + attitudes and values", "4 + believes results used")
    cat("\n=====", toupper(s), "=====\n")
    print(data.frame(rung = nm, LL = sapply(M, function(m) round(m$maximum, 1)), k = sapply(M, nfree), BIC = sapply(M, function(m) round(bic_of(m), 1)),
      p_vs_previous = c(NA, lrp(M$m1, M$m0), lrp(M$m2, M$m1), lrp(M$m3, M$m2), lrp(M$m4, M$m3)), row.names = NULL), row.names = FALSE)
  }

  cat("\nStatus quo shifters: effect on WTA in NOK per household-year (positive = more compensation needed), robust t in brackets. Rung 4.\n")
  eff <- function(m) {
    e <- m$estimate
    se <- rse(m)
    k <- c(set_exo, set_att, set_cons)
    sg <- 100
    sprintf("%7.0f (%5.1f)", sg * e[k], e[k] / se[k])
  }
  lbl <- c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "exposed", "env. organisation",
    "attitude Q7 (+1)", "concern Q10 (+1)", "visit nature Q31 (+1)", "preserve nature Q32 (+1)", "believes results used")
  print(data.frame(shifter = lbl, solar = eff(lad$solar$m4), wind = eff(lad$wind$m4)), row.names = FALSE)
  cat("\nCost sensitivity: elasticity of the cost coefficient to +1 M NOK income (robust t):\n")
  print(sapply(lad, function(M) {
    e <- M$m4$estimate
    se <- rse(M$m4)
    sprintf("%.3f (%.1f)", e["g_cinc"], e["g_cinc"] / se["g_cinc"])
  }))
})


################################################################################
# 03  Independent replication of Anders's mixed logit (his Table A.3), with my own code
################################################################################
section("03", {
  # Specification (from the article text and his replication script):
  #   WTA space, panel; 11 Normal coefficients; correlation only WITHIN an attribute (Cholesky blocks);
  #   lognormal cost coefficient that multiplies the whole utility.   30 parameters.
  # Starting values are the published estimates, so this checks whether the published optimum is reproduced
  # (different code, 500 MLHS draws instead of 2,000).
  prob_A <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    V <- list()
    V[["dev"]] <- b_cost * (b_el38 * (el_alt1 == 38) + b_el57 * (el_alt1 == 57) + b_pl3 * (plants_alt1 == 3) + b_pl6 * (plants_alt1 == 6) +
      b_nop * (wild_alt1 == 0) + b_wil * (wild_alt1 == 1) + b_d1 * (dist_alt1 == 1) + b_d2 * (dist_alt1 == 2) +
      b_d3 * (dist_alt1 == 3) + b_view * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- b_cost * b_sq
    P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
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
  dr_A <- list(interDrawsType = "mlhs", interNDraws = QK(500, 10),
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

  e <- mA$estimate
  sr <- rse(mA)
  sc <- cse(mA)
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
})


################################################################################
# 04  Mixed logit with a FULL covariance matrix (Cholesky), random cost and exogenous shifters  -  estimated in STAGES
################################################################################
section("04", {
  # Leeds [A] 1E: (i) fixing the cost coefficient imposes homogeneous scale; (ii) uncorrelated coefficients overstate heterogeneity
  # if scale varies; (iii) a Cholesky factor gives the full covariance of any set of Normals.  Anders restricts correlation to levels of the same
  # attribute and writes (footnote 16) that a fully correlated model did not converge.
  # Five random dimensions in WTA space with a FULL 5 x 5 Cholesky factor:  status quo, no protection, < 1 km, visible (Normal) and the log of the
  # cost / scale coefficient (Normal on the log scale).  Other coefficients are fixed.  The status quo constant depends on exogenous characteristics.
  # A first attempt that started the cost dimension at the published lognormal values ran away to an invalid solution
  # (log-likelihood -2232 vs -2065 published, mean of log cost 50).  Leeds [A] 1B (slide 93-96, starting values and local optima): build up in STAGES:
  #   stage 1  four Normals, full 4 x 4 Cholesky, cost fixed (b_cost estimated),             covariates on the status quo
  #   stage 2  add the cost dimension, starting at the stage-1 solution with a small, free log-cost spread
  prob_M <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    sq <- b_sq + g_age * age_c + g_fem * female + g_univ * univ + g_inc * inc_c + g_incm * inc_missing + g_exp * exposed + g_org * org
    V <- list()
    V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
      b_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + b_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
      w_d3 * (dist_alt1 == 3) + b_view * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- b_cost * sq
    P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  # stage 1: cost enters as a fixed parameter exp(m_cost) (no draw), stage 2: with a draw on the 5th dimension
  rc_M1 <- function(apollo_beta, apollo_inputs) {
    r <- list()
    r[["b_sq"]]   <- m_sq   + l11 * xi1
    r[["b_nop"]]  <- m_nop  + l21 * xi1 + l22 * xi2
    r[["b_d1"]]   <- m_d1   + l31 * xi1 + l32 * xi2 + l33 * xi3
    r[["b_view"]] <- m_view + l41 * xi1 + l42 * xi2 + l43 * xi3 + l44 * xi4
    r[["b_cost"]] <- exp(m_cost) + 0 * xi5
    return(r)
  }
  rc_M2 <- function(apollo_beta, apollo_inputs) {
    r <- list()
    r[["b_sq"]]   <- m_sq   + l11 * xi1
    r[["b_nop"]]  <- m_nop  + l21 * xi1 + l22 * xi2
    r[["b_d1"]]   <- m_d1   + l31 * xi1 + l32 * xi2 + l33 * xi3
    r[["b_view"]] <- m_view + l41 * xi1 + l42 * xi2 + l43 * xi3 + l44 * xi4
    r[["b_cost"]] <- exp(m_cost + l51 * xi1 + l52 * xi2 + l53 * xi3 + l54 * xi4 + l55 * xi5)
    return(r)
  }
  dr_M <- list(interDrawsType = "mlhs", interNDraws = QK(250, 10), interNormDraws = paste0("xi", 1:5))
  a <- A3[[SAMPLE]]
  av <- setNames(a$mean, a$par)
  cn <- c("l11", "l21", "l22", "l31", "l32", "l33", "l41", "l42", "l43", "l44", "l51", "l52", "l53", "l54", "l55")
  sf <- c("w_el38", "w_el57", "w_pl3", "w_pl6", "w_wil", "w_d2", "w_d3")
  gg <- c("g_age", "g_fem", "g_univ", "g_inc", "g_incm", "g_exp", "g_org")
  # sensible starts from the MNL baseline (step 02), NOT from the published lognormal values
  m0 <- (function() {
    e <- new.env()
    load(file.path(V2, "state_02_all.RData"), envir = e)
    e$lad[[SAMPLE]]$m0
  })()
  e0 <- m0$estimate
  st1 <- c(m_sq = e0[["w_sq"]], m_nop = e0[["w_nop"]], m_d1 = e0[["w_d1"]], m_view = e0[["w_view"]], m_cost = log(e0[["b_cost"]]),
    setNames(e0[paste0("w_", c("el38", "el57", "pl3", "pl6", "wil", "d2", "d3"))], sf),
    setNames(rep(0, 15), cn), setNames(rep(0, 7), gg))
  st1[c("l11", "l22", "l33", "l44")] <- c(30, 15, 10, 10)
  f1 <- c(cn[11:15])                                                        # stage 1: no cost dimension
  t0 <- Sys.time()
  M1 <- fit_model(paste0("mvn1_", SAMPLE), dat[[SAMPLE]], prob_M, st1, fixed = f1, draws = dr_M, randcoeff = rc_M1,
    estimate_settings = list(maxIterations = 200))
  t1 <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("Stage 1 (%s): four correlated Normals, fixed cost: LL %.1f, %d parameters, BIC %.1f (%.0f min; %s)\n", SAMPLE, M1$maximum, nfree(M1), bic_of(M1), t1, M1$message))
  st2 <- M1$estimate
  st2[c("l51", "l52", "l53", "l54")] <- 0
  st2[["l55"]] <- 0.3
  t0 <- Sys.time()
  M2 <- fit_model(paste0("mvn2_", SAMPLE), dat[[SAMPLE]], prob_M, st2, draws = dr_M, randcoeff = rc_M2, estimate_settings = list(maxIterations = 200))
  t2 <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("Stage 2 (%s): + random cost / scale dimension: LL %.1f, %d parameters, BIC %.1f (%.0f min; %s)\n", SAMPLE, M2$maximum, nfree(M2), bic_of(M2), t2, M2$message))
  cat(sprintf("Reference: Anders's published model  LL %.1f, 30 parameters, BIC %.0f\n", a$LL[1], a$BIC[1]))
  cat(sprintf("LR test stage 2 vs stage 1 (random cost): %.1f on %d df, p = %s\n", 2 * (M2$maximum - M1$maximum), nfree(M2) - nfree(M1), format(lrp(M2, M1), digits = 2)))

  mM <- if (M2$maximum > M1$maximum + 1) M2 else M1       # keep the better model, but report which
  cat("Model reported below:", if (identical(mM, M2)) "stage 2 (random cost)" else "stage 1 (cost fixed)", "\n")
  e <- mM$estimate
  sr <- rse(mM)
  L <- matrix(0, 5, 5)
  idx <- c("l11", "l21", "l31", "l41", "l51", "l22", "l32", "l42", "l52", "l33", "l43", "l53", "l44", "l54", "l55")
  L[lower.tri(L, diag = TRUE)] <- e[idx]
  Sig <- L %*% t(L)
  nm <- c("status quo", "no protection", "< 1 km", "visible", "log cost (scale)")
  dimnames(Sig) <- list(nm, nm)
  cat("\nMeans of the random coefficients (NOK; compensation = sign-flipped) and standard deviations implied by the Cholesky factor\n")
  mn <- c(sq = 100 * e[["m_sq"]], nop = -100 * e[["m_nop"]], d1 = -100 * e[["m_d1"]], view = -100 * e[["m_view"]])
  print(data.frame(coefficient = nm[1:4], mean_NOK = round(mn), robust_se = round(100 * sr[c("m_sq", "m_nop", "m_d1", "m_view")]),
    sd_NOK = round(100 * sqrt(diag(Sig)[1:4])), row.names = NULL), row.names = FALSE)
  cat(sprintf("log cost: mean %.2f (robust s.e. %.2f), sd %.2f\n", e[["m_cost"]], sr[["m_cost"]], sqrt(Sig[5, 5])))
  cat("\nImplied correlation matrix of the random coefficients\n")
  cc <- cov2cor(Sig + diag(1e-12, 5))
  print(round(cc, 2))
  cat("\nFixed coefficients, WTA in NOK (robust t)\n")
  print(data.frame(attribute = c("production 38k", "production 57k", "plants moderate", "plants many", "wilderness only", "1-3 km", "3-6 km"),
    WTA = round(-100 * e[sf]), t = round(e[sf] / sr[sf], 1), row.names = NULL), row.names = FALSE)
  cat("\nStatus quo shifters in NOK (robust t)\n")
  print(data.frame(shifter = c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "exposed", "env. organisation"),
    effect = round(100 * e[gg]), t = round(e[gg] / sr[gg], 1), row.names = NULL), row.names = FALSE)
})


################################################################################
# 05  Is wind different from solar?  A formal test inside the latent class model
################################################################################
section("05", {
  # Article, conclusion: "this result is based on separate models ... further analysis can include testing whether the differences in
  # preferences between the technologies are statistically significant."   Leeds [W] 4B (joint estimation, scale) + [A] 2A (latent class).
  # Pooled two-class model; the TRADER class has preferences w (WTA space) and the status-quo-only class chooses the status quo with fixed probability.
  #   A  traders have the SAME attribute WTAs and the same status quo constant in both technologies (wind may have its own scale, lambda);
  #      the class shares depend on the same covariates, with a wind intercept only
  #   B  = A + technology-specific status quo constants (traders and status-quo-only class)
  #   C  = B + technology-specific effects of ALL covariates on the class shares
  #   U  = the two separate LC2-DA models of step 06 (all trader WTAs technology-specific)
  # Tests:  A -> B   constants;   B -> C   who ends up in the status-quo-only class;   C -> U   do the ten attribute WTAs differ?
  e06 <- new.env()
  load(file.path(V2, "state_06_both.RData"), envir = e06)
  lc <- e06$lc
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
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    sc <- 1 + (lambda - 1) * wind                                  # solar scale = 1, wind scale = lambda
    P <- list()
    mnl <- list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice)
    V <- list()
    V[["dev"]] <- sc * b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
      w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
      w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- sc * b_cost * (w_sq + dsq_w * wind)
    mnl$utilities <- V
    mnl$componentName <- "trader"
    P[["trader"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
    mnl$utilities <- list(dev = 0, sq = a_sq2 + da_w * wind)
    mnl$componentName <- "sqonly"
    P[["sqonly"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
    P[["model"]] <- apollo_lc(list(inClassProb = P, classProb = pi_values), apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view", "b_cost")
  gn <- c("gd_age", "gd_fem", "gd_univ", "gd_inc", "gd_incm", "gd_exp", "gd_org", "ga_att", "ga_conc", "ga_natv", "ga_natp", "gc_cons")
  wi <- c("wd_age", "wd_fem", "wd_univ", "wd_inc", "wd_incm", "wd_exp", "wd_org", "wa_att", "wa_conc", "wa_natv", "wa_natp", "wc_cons")
  es <- lc$solar$lDA$estimate
  st <- c(es[c(wn, "a_sq2", "d_cls", gn)], lambda = 1, dsq_w = 0, da_w = 0, g_wind = 0.3, setNames(rep(0, 12), wi))
  fA <- c("dsq_w", "da_w", wi)
  fB <- wi
  mA <- fit_model("pool_A", dat$pooled, prob_pool, st, fixed = fA, lcpars = lcp_pool)
  mB <- fit_model("pool_B", dat$pooled, prob_pool, mA$estimate[names(st)] |> replace(c("dsq_w", "da_w"), c(0, 0)), fixed = fB, lcpars = lcp_pool)
  mC <- fit_model("pool_C", dat$pooled, prob_pool, mB$estimate[names(st)], lcpars = lcp_pool)
  LL_U <- lc$solar$lDA$maximum + lc$wind$lDA$maximum
  k_U <- nfree(lc$solar$lDA) + nfree(lc$wind$lDA)
  n_all <- lc$solar$lDA$nObs + lc$wind$lDA$nObs
  ms <- list(mA, mB, mC)
  tb <- data.frame(model = c("A  common WTA, common constants", "B  + technology-specific constants", "C  + technology-specific class-share effects", "U  all separate (two LC2-DA models)"),
    LL = round(c(sapply(ms, function(m) m$maximum), LL_U), 1), k = c(sapply(ms, nfree), k_U))
  tb$BIC <- round(-2 * tb$LL + tb$k * log(n_all), 1)
  print(tb, row.names = FALSE)
  lr <- function(l1, k1, l0, k0) {
    x <- 2 * (l1 - l0)
    data.frame(LR = round(x, 1), df = k1 - k0, p = signif(pchisq(x, k1 - k0, lower.tail = FALSE), 2))
  }
  cat("\nLikelihood-ratio tests\n")
  print(rbind("A -> B  constants differ between technologies?" = lr(mB$maximum, nfree(mB), mA$maximum, nfree(mA)),
    "B -> C  covariates act differently on class shares?" = lr(mC$maximum, nfree(mC), mB$maximum, nfree(mB)),
    "C -> U  do the 10 attribute WTAs of the traders differ?" = lr(LL_U, k_U, mC$maximum, nfree(mC)),
    "A -> U  everything" = lr(LL_U, k_U, mA$maximum, nfree(mA))))
  e <- mB$estimate
  se <- rse(mB)
  cat(sprintf("\nModel B: wind scale lambda = %.2f (robust s.e. %.2f); wind intercept in the class-share equation %.2f (z %.1f)\n", e[["lambda"]], se[["lambda"]], e[["g_wind"]], e[["g_wind"]] / se[["g_wind"]]))
  cat(sprintf("Status quo constant of traders: solar %.0f NOK, wind %.0f NOK (difference %.0f, z %.1f);  P(status quo | status-quo-only class): solar %.3f, wind %.3f\n",
    100 * e[["w_sq"]], 100 * (e[["w_sq"]] + e[["dsq_w"]]), 100 * e[["dsq_w"]], e[["dsq_w"]] / se[["dsq_w"]], plogis(e[["a_sq2"]]), plogis(e[["a_sq2"]] + e[["da_w"]])))
})


################################################################################
# 06  Latent class family: who never trades, who ignores cost, what do traders want?  (+ a formal solar-vs-wind test)
################################################################################
section("06", {
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
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list()
    mnl <- list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice)
    V <- list()
    V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
      w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
      w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- b_cost * w_sq
    mnl$utilities <- V
    mnl$componentName <- "trader"
    P[["trader"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
    mnl$utilities <- list(dev = 0, sq = a_sq2)
    mnl$componentName <- "sqonly"
    P[["sqonly"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
    P[["model"]] <- apollo_lc(list(inClassProb = P, classProb = pi_values), apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
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
    st1 <- l0$estimate
    st1[c(set_gd, set_ga)] <- 0
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
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list()
    mnl <- list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice)
    V <- list()
    V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
      w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
      w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- b_cost * w_sq
    mnl$utilities <- V
    mnl$componentName <- "trader"
    P[["trader"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
    # cost-ignorers: the same attributes, but NO cost term (stated attribute non-attendance to cost); utility units
    V <- list()
    V[["dev"]] <- e_el38 * (el_alt1 == 38) + e_el57 * (el_alt1 == 57) + e_pl3 * (plants_alt1 == 3) + e_pl6 * (plants_alt1 == 6) +
      e_nop * (wild_alt1 == 0) + e_wil * (wild_alt1 == 1) + e_d1 * (dist_alt1 == 1) + e_d2 * (dist_alt1 == 2) +
      e_d3 * (dist_alt1 == 3) + e_view * (view_alt1 == 2)
    V[["sq"]] <- a_sq_c
    mnl$utilities <- V
    mnl$componentName <- "costign"
    P[["costign"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
    mnl$utilities <- list(dev = 0, sq = a_sq2)
    mnl$componentName <- "sqonly"
    P[["sqonly"]] <- apollo_panelProd(apollo_mnl(mnl, functionality), apollo_inputs, functionality)
    P[["model"]] <- apollo_lc(list(inClassProb = P, classProb = pi_values), apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  hn <- c(paste0("h2_", c("age", "fem", "att", "conc", "natp", "cons", "icost")), paste0("h3_", c("age", "fem", "att", "conc", "natp", "cons", "icost")))
  en <- c("e_el38", "e_el57", "e_pl3", "e_pl6", "e_nop", "e_wil", "e_d1", "e_d2", "e_d3", "e_view", "a_sq_c")
  lc3 <- list()
  for (s in c("solar", "wind")) {
    d <- dat[[s]]
    e0 <- lc[[s]]$lDA$estimate
    mnl_pref <- c(e_el38 = -0.1, e_el57 = -0.1, e_pl3 = -0.1, e_pl6 = -0.2, e_nop = -0.5, e_wil = -0.2, e_d1 = -0.3, e_d2 = -0.2, e_d3 = -0.1, e_view = -0.3, a_sq_c = 0.5)
    st <- c(e0[c(wn, "a_sq2")], mnl_pref, d2 = -1, d3 = 0, setNames(rep(0, 14), hn))
    lc3[[s]] <- try(fit_model(paste0("lc3_", s), d, prob_lc3, st, lcpars = lcp3), silent = TRUE)
  }

  # ---------------- report -------------
  cat("Model comparison\n")
  for (s in c("solar", "wind")) {
    M <- lc[[s]]
    m3 <- lc3[[s]]
    ok3 <- !inherits(m3, "try-error")
    nm <- c("LC2-0  constant shares", "LC2-D  shares: exogenous covariates", "LC2-DA shares: + attitudes, Q41", if (ok3) "LC3    + cost-ignorer class")
    ms <- c(list(M$l0, M$lD, M$lDA), if (ok3) list(m3))
    cat("\n=====", toupper(s), "=====\n")
    print(data.frame(model = nm, LL = sapply(ms, function(m) round(m$maximum, 1)), k = sapply(ms, nfree), BIC = sapply(ms, function(m) round(bic_of(m), 1)),
      p_vs_previous = c(NA, lrp(M$lD, M$l0), lrp(M$lDA, M$lD), if (ok3) lrp(m3, M$lDA)), row.names = NULL), row.names = FALSE)
  }
  cat("\nPopulation-weighted vs unweighted share in the status-quo-only class, model LC2-DA\n")
  shr <- function(m, s) {
    e <- m$estimate
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    p <- plogis(e[["d_cls"]] + e[["gd_age"]] * r$age_c + e[["gd_fem"]] * r$female + e[["gd_univ"]] * r$univ + e[["gd_inc"]] * r$inc_c + e[["gd_incm"]] * r$inc_missing +
      e[["gd_exp"]] * r$exposed + e[["gd_org"]] * r$org + e[["ga_att"]] * r$att_c + e[["ga_conc"]] * r$conc_c + e[["ga_natv"]] * r$natv_c + e[["ga_natp"]] * r$natp_c + e[["gc_cons"]] * r$conseq)
    c(unweighted = mean(p), weighted = sum(p * r$wt) / sum(r$wt))
  }
  print(round(sapply(c(solar = "solar", wind = "wind"), function(s) shr(lc[[s]]$lDA, s)), 3))
  cat("(raw share of respondents who ALWAYS chose the status quo: solar 0.359, wind 0.520)\n")

  cat("\nMembership of the status-quo-only class (LC2-DA): odds ratios, robust z in brackets\n")
  gn <- c(set_gd, set_ga)
  lab <- c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "exposed", "env. organisation",
    "attitude Q7 (+1)", "concern Q10 (+1)", "visit nature Q31 (+1)", "preserve nature Q32 (+1)", "believes results used")
  fo <- function(m) {
    e <- m$estimate
    se <- rse(m)
    sprintf("%5.2f (%5.1f)", exp(e[gn]), e[gn] / se[gn])
  }
  print(data.frame(variable = lab, solar = fo(lc$solar$lDA), wind = fo(lc$wind$lDA)), row.names = FALSE)

  cat("\nWTA of the TRADER class, NOK per household-year (robust s.e.); LC2-DA\n")
  wt <- function(m) {
    e <- m$estimate
    se <- rse(m)
    k <- setdiff(wn, "b_cost")
    sg <- ifelse(k == "w_sq", 100, -100)
    data.frame(attr = sub("^w_", "", k), WTA = round(sg * e[k]), se = round(100 * se[k]))
  }
  a <- wt(lc$solar$lDA)
  b <- wt(lc$wind$lDA)
  print(data.frame(attribute = a$attr, solar = a$WTA, se_s = a$se, wind = b$WTA, se_w = b$se), row.names = FALSE)
  cat(sprintf("\nStatus-quo-only class: P(status quo) = %.3f (solar), %.3f (wind)\n", plogis(lc$solar$lDA$estimate[["a_sq2"]]), plogis(lc$wind$lDA$estimate[["a_sq2"]])))

  for (s in c("solar", "wind")) if (!inherits(lc3[[s]], "try-error")) {
    m3 <- lc3[[s]]
    e <- m3$estimate
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    v2 <- e[["d2"]] + e[["h2_age"]] * r$age_c + e[["h2_fem"]] * r$female + e[["h2_att"]] * r$att_c + e[["h2_conc"]] * r$conc_c + e[["h2_natp"]] * r$natp_c + e[["h2_cons"]] * r$conseq + e[["h2_icost"]] * r$imp_cost_c
    v3 <- e[["d3"]] + e[["h3_age"]] * r$age_c + e[["h3_fem"]] * r$female + e[["h3_att"]] * r$att_c + e[["h3_conc"]] * r$conc_c + e[["h3_natp"]] * r$natp_c + e[["h3_cons"]] * r$conseq + e[["h3_icost"]] * r$imp_cost_c
    den <- 1 + exp(v2) + exp(v3)
    cat(sprintf("\n[%s LC3] class shares: traders %.2f, cost-ignorers %.2f, status-quo-only %.2f;  P(status quo | cost-ignorers, reference scenario) = %.2f\n",
      s, mean(1 / den), mean(exp(v2) / den), mean(exp(v3) / den), plogis(e[["a_sq_c"]])))
    cat("  cost-ignorers, utility weights: ", paste(sub("^e_", "", names(e)[grepl("^e_", names(e))]), round(e[grepl("^e_", names(e))], 2), collapse = "; "), "\n")
    cat("  membership of the cost-ignorer class, odds ratio (z): stated cost importance (Q25) ", sprintf("%.2f (%.1f)", exp(e[["h2_icost"]]), e[["h2_icost"]] / rse(m3)[["h2_icost"]]),
      "; attitude ", sprintf("%.2f (%.1f)", exp(e[["h2_att"]]), e[["h2_att"]] / rse(m3)[["h2_att"]]), "\n")
  }
})


################################################################################
# 06b  Three-class model with a cost-RESPONSIVE trader class (repair of a degenerate solution in the wind sample)
################################################################################
section("06b", {
  # In step 06 the wind LC3 solution let the 'trader' class drift to b_cost ~ 0.0002 with WTAs in the millions: in WTA space this is a class that
  # ignores cost (utility-space weights b_cost * w stay finite while w -> infinity), i.e. a second cost-ignorer class. Leeds [A] 2A (slide 17):
  # "with more than two classes the risk increases of some classes collapsing together ... fine tuning can help, e.g. merging or restricting parameters".
  # Remedy: keep the trader class cost-responsive by fixing its cost coefficient at the value of the two-class model (LC2-DA); everything else is free.
  # For solar the unrestricted three-class model is valid; it is re-estimated here with the same restriction so that both samples are treated alike.
  e06 <- new.env()
  load(file.path(V2, "state_06_both.RData"), envir = e06)
  lc <- e06$lc
  lc3 <- e06$lc3
  prob_lc3 <- e06$prob_lc3
  lcp3 <- e06$lcp3
  fx_cost <- function(s) lc[[s]]$lDA$estimate[["b_cost"]]
  res <- list()
  for (s in c("solar", "wind")) {
    st <- lc3[[s]]$estimate
    wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view")
    st[wn] <- lc[[s]]$lDA$estimate[wn]
    st[["b_cost"]] <- fx_cost(s)
    res[[s]] <- fit_model(paste0("lc3c_", s), dat[[s]], prob_lc3, st, fixed = "b_cost", lcpars = lcp3)
  }
  cat("Three classes with a cost-responsive trader class (b_cost fixed at its two-class value)\n")
  for (s in c("solar", "wind")) {
    m <- res[[s]]
    l2 <- lc[[s]]$lDA
    cat(sprintf("\n===== %s =====\nLC2-DA: LL %.1f, k %d, BIC %.1f   |   LC3 (restricted): LL %.1f, k %d, BIC %.1f   (unrestricted LC3: LL %.1f)\n", toupper(s),
      l2$maximum, nfree(l2), bic_of(l2), m$maximum, nfree(m), bic_of(m), lc3[[s]]$maximum))
    cat(sprintf("LR test LC3 vs LC2-DA: %.1f on %d df, p = %s\n", 2 * (m$maximum - l2$maximum), nfree(m) - nfree(l2), format(lrp(m, l2), digits = 2)))
    e <- m$estimate
    se <- rse(m)
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    v2 <- e[["d2"]] + e[["h2_age"]] * r$age_c + e[["h2_fem"]] * r$female + e[["h2_att"]] * r$att_c + e[["h2_conc"]] * r$conc_c + e[["h2_natp"]] * r$natp_c + e[["h2_cons"]] * r$conseq + e[["h2_icost"]] * r$imp_cost_c
    v3 <- e[["d3"]] + e[["h3_age"]] * r$age_c + e[["h3_fem"]] * r$female + e[["h3_att"]] * r$att_c + e[["h3_conc"]] * r$conc_c + e[["h3_natp"]] * r$natp_c + e[["h3_cons"]] * r$conseq + e[["h3_icost"]] * r$imp_cost_c
    den <- 1 + exp(v2) + exp(v3)
    sh <- cbind(trader = 1 / den, cost_ignorer = exp(v2) / den, sq_only = exp(v3) / den)
    cat("Class shares (unweighted | population-weighted):\n")
    wtv <- as.numeric(r$wt)
    print(round(rbind(unweighted = colMeans(sh), weighted = colSums(sh * wtv) / sum(wtv)), 3))
    wn2 <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view")
    sg <- c(100, rep(-100, 10))
    cat("Trader-class WTA, NOK (robust s.e.):\n")
    print(data.frame(attribute = sub("w_", "", wn2), WTA = round(sg * e[wn2]), se = round(100 * se[wn2])), row.names = FALSE)
    en <- grep("^e_", names(e), value = TRUE)
    cat("Cost-ignorer class, utility weights (robust t):\n")
    print(data.frame(par = sub("^e_", "", en), weight = round(e[en], 2), t = round(e[en] / se[en], 1)), row.names = FALSE)
    cat(sprintf("Status-quo-only class: P(status quo) = %.3f.  Cost-ignorers: P(status quo) at the reference scenario = %.3f\n", plogis(e[["a_sq2"]]), plogis(e[["a_sq_c"]])))
    cat(sprintf("Membership of the cost-ignorer class, odds ratio per +1 point of STATED cost importance (Q25): %.2f (z %.1f);  attitude Q7: %.2f (z %.1f)\n",
      exp(e[["h2_icost"]]), e[["h2_icost"]] / se[["h2_icost"]], exp(e[["h2_att"]]), e[["h2_att"]] / se[["h2_att"]]))
  }
})


################################################################################
# 07  Hybrid choice (ICLV) with FIVE indicators and demographics in the structural equation
################################################################################
section("07", {
  # The survey measures one underlying "concern about nature / resistance to development" with five answers:
  #   Q7  general attitude to development (reverse), Q10 concern about nature effects, Q31 importance of visiting untouched nature,
  #   Q32 should untouched nature be kept at today's level (3 levels), Q33 perceived negative landscape impact of the technology in question.
  # Structural equation:  LV = gamma' z + eta,  eta ~ N(0,1),   z = age, gender, degree, income, exposure, organisation (EXOGENOUS)
  # Measurement:          each indicator ~ ordered logit( zeta_i * LV ),  zeta_conc > 0 anchors the sign (higher LV = more concern)
  # Choice:               V_sq = b_cost * ( w_sq + lambda * LV ),   V_dev has  w_nop + lambda2 * LV  on the no-protection dummy
  # Following the Leeds slides (3C, p. 24-29): the log-likelihood includes the indicators and must NOT be compared with choice-only models;
  # the safe way to use the model is through the covariates z (reweighting a sample, a future population), not by moving LV by "+-1 SD".
  d <- dat[[SAMPLE]]
  d$nat_pres_f <- ifelse(is.na(d$nat_pres), 2, d$nat_pres)
  d$nat_pres_ok <- as.numeric(!is.na(d$nat_pres))
  d$impact_f <- ifelse(is.na(d$impact_own), 4, d$impact_own)
  d$impact_ok <- as.numeric(!is.na(d$impact_own))
  d$att_f <- ifelse(is.na(d$att), 4, d$att)
  d$conc_f <- ifelse(is.na(d$conc), 4, d$conc)
  d$natv_f <- ifelse(is.na(d$nat_visit), 4, d$nat_visit)

  prob_H <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    P <- list()
    # ---- thresholds built from positive increments (they stay ordered)
    ta2 <- ta1 + exp(ta2i)
    ta3 <- ta2 + exp(ta3i)
    ta4 <- ta3 + exp(ta4i)
    ta5 <- ta4 + exp(ta5i)
    ta6 <- ta5 + exp(ta6i)
    tc2 <- tc1 + exp(tc2i)
    tc3 <- tc2 + exp(tc3i)
    tc4 <- tc3 + exp(tc4i)
    tc5 <- tc4 + exp(tc5i)
    tc6 <- tc5 + exp(tc6i)
    tv2 <- tv1 + exp(tv2i)
    tv3 <- tv2 + exp(tv3i)
    tv4 <- tv3 + exp(tv4i)
    tv5 <- tv4 + exp(tv5i)
    tv6 <- tv5 + exp(tv6i)
    ti2 <- ti1 + exp(ti2i)
    ti3 <- ti2 + exp(ti3i)
    ti4 <- ti3 + exp(ti4i)
    ti5 <- ti4 + exp(ti5i)
    ti6 <- ti5 + exp(ti6i)
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
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  rc_H <- function(apollo_beta, apollo_inputs) {
    r <- list()
    r[["LV"]] <- g_age * age_c + g_fem * female + g_univ * univ + g_inc * inc_c + g_incm * inc_missing + g_exp * exposed + g_org * org + eta
    return(r)
  }
  dr_H <- list(interDrawsType = "halton", interNDraws = QK(120, 10), interNormDraws = "eta")

  wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view", "b_cost")
  m0 <- fit_model(paste0("hmnl_", SAMPLE), d, function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    V <- list()
    V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
      w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) + w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- b_cost * w_sq
    P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  },
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

  e <- mH$estimate
  V <- mH$robvarcov
  se <- rse(mH)
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
  dm <- function(a, b) {
    g <- c(e[[b]], e[[a]])
    v <- V[c(a, b), c(a, b)]
    c(est = e[[a]] * e[[b]], se = sqrt(drop(t(g) %*% v %*% g)))
  }
  cat("\nImplied effect on the status quo WTA through the latent variable, NOK per household-year  (delta method, robust)\n")
  imp <- t(sapply(gz, function(g) dm("lambda", g) * c(100, 100)))
  print(data.frame(characteristic = lab, effect_NOK = round(imp[, "est"]), se = round(imp[, "se"]), z = round(imp[, "est"] / imp[, "se"], 1), row.names = NULL), row.names = FALSE)
})


################################################################################
# 08  Bayesian (hierarchical) estimation of the mixed logit with a FULL 12 x 12 covariance matrix
################################################################################
section("08", {
  # Leeds [A] 3A.  Footnote 16 of the article: the fully correlated model "did not converge" with simulated maximum likelihood.
  # Bayesian estimation (Gibbs + Metropolis-Hastings, package RSGHB through Apollo) does not need to maximise a likelihood, so a full
  # covariance matrix is feasible. Leeds 1B (slide 3): "HB is NOT a model" - the model is the same mixed logit; only the estimator differs.
  # Units: coefficients in 1,000 NOK (so that the Metropolis-Hastings step sizes work; slides 69-70 warn about scaling).
  #   V_dev = b_cost * ( sum_k w_k x_k + cost / 1000 ),   V_sq = b_cost * w_sq,   w_k ~ Normal (full covariance),  b_cost ~ positive lognormal
  d <- dat[[SAMPLE]]
  d$agegrp <- as.character(d$agegrp)
  nm <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view")
  database <- d
  quiet(apollo_initialise())
  apollo_control <- list(modelName = paste0("hb_", SAMPLE), indivID = "id", HB = TRUE, nCores = 1, outputDirectory = outdir)
  apollo_beta <- c(setNames(c(3.9, -0.9, -1.4, -0.8, -1.7, -7.4, -1.8, -4.4, -2.6, -0.6, -4.7), nm), b_cost = 1)
  apollo_fixed <- c()
  apollo_HB <- list(hbDist = c(setNames(rep("N", 11), nm), b_cost = "LN+"), gNCREP = QK(40000, 400), gNEREP = QK(20000, 200), gINFOSKIP = 5000, gFULLCV = TRUE, gSeed = 2026)
  apollo_inputs <- quiet(apollo_validateInputs())
  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    V <- list()
    V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
      w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
      w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 1000)
    V[["sq"]] <- b_cost * w_sq
    P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)     # no apollo_panelProd under HB
  }
  t0 <- Sys.time()
  hb <- quiet(apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs))
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("Bayesian mixed logit, %s sample: %d burn-in + %d retained iterations (%.0f min)\n", SAMPLE, 40000, 20000, mins))

  rn <- hb$HB_names_random_params
  M  <- hb$HB_iterations_means
  cn <- colnames(M)
  cat("columns of the stored chain of means:", paste(cn, collapse = ", "), "\n")
  K  <- hb$HB_iterations_covar
  ck <- sub("^.*_", "", "")
  invisible(ck)
  pick <- function(x) M[, grep(paste0("^", x, "$"), cn)[1]]
  post <- sapply(c(nm, "b_cost"), function(v) {
    j <- which(cn == v | cn == paste0("mu_", v) | cn == paste0(v, "_mean"))
    if (length(j)) M[, j[1]] else rep(NA, nrow(M))
  })
  if (anyNA(post)) post <- M[, 1:12, drop = FALSE]
  colnames(post) <- c(nm, "b_cost")[seq_len(ncol(post))]
  sg <- c(w_sq = 1000, setNames(rep(-1000, 10), nm[-1]))
  tab <- data.frame(attribute = nm, WTA_post_mean = round(sg * colMeans(post[, nm])), post_sd = round(1000 * apply(post[, nm], 2, sd)),
    q025 = round(apply(sg * post[, nm], 2, function(z) quantile(z, if (TRUE) 0.025))), q975 = round(apply(sg * post[, nm], 2, function(z) quantile(z, 0.975))), row.names = NULL)
  # for the negative-signed attributes the percentile order flips; recompute properly
  lo <- sapply(nm, function(v) min(quantile(sg[[v]] * post[, v], c(.025, .975))))
  hi <- sapply(nm, function(v) max(quantile(sg[[v]] * post[, v], c(.025, .975))))
  tab$q025 <- round(lo)
  tab$q975 <- round(hi)
  a <- A3[[SAMPLE]]
  sg_a <- c(100, rep(-100, 10))
  tab$Anders_published <- round(sg_a * a$mean[1:11])
  cat("\nPosterior mean of the MEAN WTA (NOK per household-year), posterior sd and 95% credible interval; Anders's published value for comparison\n")
  print(tab, row.names = FALSE)

  # population standard deviations and correlations of the random coefficients
  cv <- apply(K, c(1, 2), mean)
  sdv <- sqrt(diag(cv))
  names(sdv) <- rn
  cat("\nPosterior mean standard deviation across respondents (NOK):\n")
  print(round(1000 * sdv[1:11]))
  cr <- cov2cor(cv)
  dimnames(cr) <- list(rn, rn)
  cat("\nPosterior mean correlations of the status quo constant with the other random coefficients (and with the cost / scale coefficient):\n")
  print(round(cr[rn[1], -1], 2))
  cat("\nPosterior mean correlation matrix, attribute block (rounded):\n")
  print(round(cr[rn[1:11], rn[1:11]], 1))
  g <- hb$HB_Geweke_test_means
  cat("\nConvergence: Geweke z of the chains of the means (|z| < 2 is fine):\n")
  print(round(g, 2))
})


################################################################################
# 09  09  Cost per kWh with the delta method (as in the article), across models, and a population-weighted acceptance curve
################################################################################
POLICY <- quote({
  # Anders: cost/kWh = 100 * ( w_sq - sum of the w's of the scenario ) NOK * 4,100 households / 600,000,000 kWh, in ore; delta method on the robust covariance.
  # (A) reproduce his Figure 6 from the published Table A.3                              (checks the formula and the units)
  # (B) the same formula, delta method, for every model of this re-analysis
  # (C) acceptance curves, post-stratified to the population (Leeds [A] 1B slides 54-56: weighting AFTER estimation, by sample enumeration)
  ld <- function(f) {
    e <- new.env()
    load(file.path(V2, f), envir = e)
    e
  }
  e02 <- ld("state_02_all.RData")
  lad <- e02$lad
  e06 <- ld("state_06_both.RData")
  lc <- e06$lc
  rep_ <- lapply(c(solar = "solar", wind = "wind"), function(s) ld(sprintf("state_03_%s.RData", s))$mA)
  mvn  <- lapply(c(solar = "solar", wind = "wind"), function(s) ld(sprintf("state_04_%s.RData", s))$mM)

  H <- 4100
  KWH <- 600e6
  scen <- list(P1 = c("el38", "d1"), P2 = c("el38", "pl6", "d1"), P3 = c("el38", "nop"), P4 = c("el38", "pl6", "nop"), P5 = c("el38", "pl3", "wil", "d3"), P6 = "el38")
  fig6 <- list(wind  = rbind(vis = c(17.6, 18.2, 14.7, 15.3, 14.2, 11.1), novis = c(15.7, 16.3, 12.8, 13.3, 12.3, 9.2)),
    solar = rbind(vis = c(9.5, 10.7, 11.6, 12.8, 8.6, 6.5),   novis = c(6.3, 7.5, 8.4, 9.6, 5.4, 3.3)))
  # generic cost: z = list(sq = status quo constant, w = named vector of attribute w's, view)
  ore <- function(z, lev, visible) 100 * (z$sq - sum(z$w[lev]) - if (visible) z$view else 0) * H / KWH * 100

  # ------------------------------------------------------------------ (A) Figure 6 from Table A.3
  cat("(A) Figure 6 recomputed from the published Table A.3 (ore/kWh) vs the figure in the article\n")
  chk <- do.call(rbind, lapply(c("wind", "solar"), function(s) {
    a <- A3[[s]]
    z <- list(sq = a$mean[1], w = setNames(a$mean[2:11], a$par[2:11]), view = a$mean[11])
    z$w <- setNames(a$mean[2:10], c("el38", "el57", "pl3", "pl6", "nop", "wil", "d1", "d2", "d3"))
    z$view <- a$mean[11]
    do.call(rbind, lapply(c(TRUE, FALSE), function(v) data.frame(sample = s, visible = v, scenario = names(scen),
      recomputed = round(sapply(scen, function(l) ore(z, l, v)), 1), article = fig6[[s]][ifelse(v, "vis", "novis"), ])))
  }))
  chk$diff <- round(chk$recomputed - chk$article, 1)
  print(reshape(chk[, c("sample", "visible", "scenario", "recomputed", "article")], idvar = c("scenario", "visible"), timevar = "sample", direction = "wide"), row.names = FALSE)
  cat("max absolute difference:", max(abs(chk$diff)), "ore/kWh\n")
  cat("\nWorked example, Table A.6 of the article (alternative 1 with visibility, i.e. < 1 km, all protected, few plants):\n")
  for (s in c("solar", "wind")) {
    a <- A3[[s]]
    el <- a$mean[2]
    el57 <- a$mean[3]
    d1 <- a$mean[8]
    vw <- a$mean[11]     # index 8 = < 1 km, 11 = visible
    s38 <- 100 * (a$mean[1] - el - d1 - vw)
    s57 <- 100 * (a$mean[1] - el57 - d1 - vw)
    cat(sprintf("  %-5s with the el38 coefficient (%.0f NOK): sum %.0f NOK -> %.2f ore/kWh;  with the el57 coefficient (%.0f NOK): sum %.0f NOK -> %.2f ore/kWh\n", s,
      -100 * el, s38, s38 * H / KWH * 100, -100 * el57, s57, s57 * H / KWH * 100))
  }

  # ------------------------------------------------------------------ delta method (numerical Jacobian; identical to car::deltaMethod for smooth functions)
  dm <- function(m, f) {
    th <- m$estimate
    V <- m$robvarcov
    nmv <- rownames(V)
    J <- sapply(nmv, function(n) {
      h <- 1e-5 * max(1, abs(th[[n]]))
      tp <- th
      tp[n] <- th[[n]] + h
      tn <- th
      tn[n] <- th[[n]] - h
      (f(tp) - f(tn)) / (2 * h)
    })
    c(est = f(th), se = sqrt(drop(J %*% V %*% J)))
  }
  wbar <- function(s, v) {
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    sum(r[[v]] * r$wt) / sum(r$wt)
  }
  zb <- function(s) sapply(c("age_c", "female", "univ", "inc_c", "inc_missing", "exposed", "org"), function(v) wbar(s, v))
  mk_z <- list(
    MNL    = function(th, s) list(sq = th[["w_sq"]], w = c(el38 = th[["w_el38"]], el57 = th[["w_el57"]], pl3 = th[["w_pl3"]], pl6 = th[["w_pl6"]], nop = th[["w_nop"]], wil = th[["w_wil"]], d1 = th[["w_d1"]], d2 = th[["w_d2"]], d3 = th[["w_d3"]]), view = th[["w_view"]]),
    LCtr   = function(th, s) list(sq = th[["w_sq"]], w = c(el38 = th[["w_el38"]], el57 = th[["w_el57"]], pl3 = th[["w_pl3"]], pl6 = th[["w_pl6"]], nop = th[["w_nop"]], wil = th[["w_wil"]], d1 = th[["w_d1"]], d2 = th[["w_d2"]], d3 = th[["w_d3"]]), view = th[["w_view"]]),
    REP    = function(th, s) list(sq = th[["m_sq"]], w = c(el38 = th[["m_el38"]], el57 = th[["m_el57"]], pl3 = th[["m_pl3"]], pl6 = th[["m_pl6"]], nop = th[["m_nop"]], wil = th[["m_wil"]], d1 = th[["m_d1"]], d2 = th[["m_d2"]], d3 = th[["m_d3"]]), view = th[["m_view"]]),
    MVN    = function(th, s) {
      z <- zb(s)
      g <- c(th[["g_age"]], th[["g_fem"]], th[["g_univ"]], th[["g_inc"]], th[["g_incm"]], th[["g_exp"]], th[["g_org"]])
      list(sq = th[["m_sq"]] + sum(g * z), w = c(el38 = th[["w_el38"]], el57 = th[["w_el57"]], pl3 = th[["w_pl3"]], pl6 = th[["w_pl6"]], nop = th[["m_nop"]], wil = th[["w_wil"]], d1 = th[["m_d1"]], d2 = th[["w_d2"]], d3 = th[["w_d3"]]), view = th[["m_view"]])
    })
  models <- function(s) list(MNL = list(lad[[s]]$m0, "MNL"), `LC traders (LC2-DA)` = list(lc[[s]]$lDA, "LCtr"), `MMNL replication of Anders` = list(rep_[[s]], "REP"), `MMNL full covariance` = list(mvn[[s]], "MVN"))

  cat("\n(B) Cost per kWh (ore), delta-method s.e. in brackets, for every model; scenario x visibility\n")
  tabB <- list()
  for (s in c("solar", "wind")) {
    rows <- list()
    for (mn in names(models(s))) {
      m <- models(s)[[mn]][[1]]
      ty <- models(s)[[mn]][[2]]
      for (v in c(FALSE, TRUE)) rows[[length(rows) + 1]] <- c(model = mn, visible = ifelse(v, "yes", "no"),
        sapply(names(scen), function(p) {
          r <- dm(m, function(th) ore(mk_z[[ty]](th, s), scen[[p]], v))
          sprintf("%5.1f (%.1f)", r["est"], r["se"])
        }))
    }
    tabB[[s]] <- as.data.frame(do.call(rbind, rows))
    cat("\n=====", toupper(s), "=====\n")
    print(tabB[[s]], row.names = FALSE)
  }

  cat("\nRatio of the status-quo term, wind / solar ('3.5 times higher' in the article), delta method (independent samples):\n")
  for (mn in c("MNL", "MMNL replication of Anders")) {
    ty <- models("solar")[[mn]][[2]]
    a <- dm(models("solar")[[mn]][[1]], function(th) mk_z[[ty]](th, "solar")$sq)
    b <- dm(models("wind")[[mn]][[1]], function(th) mk_z[[ty]](th, "wind")$sq)
    r <- b["est"] / a["est"]
    se <- r * sqrt((b["se"] / b["est"])^2 + (a["se"] / a["est"])^2)
    cat(sprintf("  %-28s solar %.0f NOK (se %.0f), wind %.0f NOK (se %.0f); ratio %.2f (delta se %.2f)\n", mn, 100 * a["est"], 100 * a["se"], 100 * b["est"], 100 * b["se"], r, se))
  }

  # ------------------------------------------------------------------ (C) acceptance curves from the exogenous-only latent class model, population-weighted
  cat("\n(C) Share of households accepting a scenario at a given annual bill reduction, LC2-D (class shares by exogenous characteristics only)\n")
  accept <- function(th, s, lev, visible, cvals, weighted) {
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    v2 <- th[["d_cls"]] + th[["gd_age"]] * r$age_c + th[["gd_fem"]] * r$female + th[["gd_univ"]] * r$univ + th[["gd_inc"]] * r$inc_c + th[["gd_incm"]] * r$inc_missing + th[["gd_exp"]] * r$exposed + th[["gd_org"]] * r$org
    pi2 <- plogis(v2)
    wts <- if (weighted) r$wt else rep(1, nrow(r))
    z <- mk_z$LCtr(th, s)
    Sw <- sum(z$w[lev]) + if (visible) z$view else 0
    sapply(cvals, function(cc) sum(wts * ((1 - pi2) * plogis(th[["b_cost"]] * (Sw + cc / 100 - z$sq)) + pi2 * (1 - plogis(th[["a_sq2"]])))) / sum(wts))
  }
  cv <- c(0, 2000, 5000, 9000, 15000, 25000, 50000)
  for (s in c("solar", "wind")) {
    m <- lc[[s]]$lD
    for (p in c("P6", "P1", "P4")) {
      out <- do.call(cbind, lapply(c(FALSE, TRUE), function(w) sapply(seq_along(cv), function(i) {
        r <- dm(m, function(th) accept(th, s, scen[[p]], FALSE, cv[i], w))
        c(r[["est"]], r[["se"]])
      })))
      cat(sprintf("\n[%s, %s, not visible]  compensation (NOK) | accepted share, population-weighted (delta s.e.) | unweighted\n", toupper(s), p))
      print(data.frame(NOK = cv, weighted = sprintf("%4.1f%% (%.1f)", 100 * out[1, seq_along(cv)], 100 * out[2, seq_along(cv)]),
        unweighted = sprintf("%4.1f%% (%.1f)", 100 * out[1, length(cv) + seq_along(cv)], 100 * out[2, length(cv) + seq_along(cv)])), row.names = FALSE)
    }
  }
  cat("\nShare of households in the status-quo-only class (LC2-D), weighted vs unweighted:\n")
  print(round(sapply(c(solar = "solar", wind = "wind"), function(s) {
    th <- lc[[s]]$lD$estimate
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    p <- plogis(th[["d_cls"]] + th[["gd_age"]] * r$age_c + th[["gd_fem"]] * r$female + th[["gd_univ"]] * r$univ + th[["gd_inc"]] * r$inc_c + th[["gd_incm"]] * r$inc_missing + th[["gd_exp"]] * r$exposed + th[["gd_org"]] * r$org)
    c(unweighted = mean(p), weighted = sum(p * r$wt) / sum(r$wt))
  }), 3))
})

section("09", POLICY)


################################################################################
# 10  Robustness: who is removed, and does belief in the survey's consequences matter?
################################################################################
section("10", {
  # Anders removes respondents who gave one of nine reasons for always choosing the status quo ("protest").  The reason that is by far the most
  # common ("Supplying Europe with renewable energy from Norway must not come at the expense of Norwegian nature") expresses a PROTECTED VALUE
  # for nature, not a rejection of the survey.  Leeds [A] 3E: non-compensatory (lexicographic) decision rules are a recognised behaviour.
  # Variants of the MNL baseline in WTA space (no weights in estimation - Leeds [A] 1B, slides 51-54):
  #   all | Anders's protest rule | protest rule WITHOUT the 'protected value' reasons | no always-status-quo respondents |
  #   only those who believe the results will be used (Q41) | without the 10% fastest | without those who found the cards difficult (6-7)
  pat_anders <- paste("bekostning av norsk", "^Andre", "ikke bidra tilstrekkelig", "gr.nt batteri", "uansett eksportert", "vil fungere etter", "^Vet ikke$",
    "for vanskelige", "visste ikke hvilket", sep = "|")
  pat_noval  <- paste("^Andre", "ikke bidra tilstrekkelig", "gr.nt batteri", "uansett eksportert", "vil fungere etter", "^Vet ikke$", "for vanskelige", "visste ikke hvilket", sep = "|")
  wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view", "b_cost")
  prob_b <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    V <- list()
    V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
      w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
      w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- b_cost * w_sq
    P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  st0 <- setNames(c(rep(0, 11), 0.005), wn)
  key <- c("w_sq", "w_nop", "w_d1", "w_view")
  rob <- list()
  for (s in c("solar", "wind")) {
    d <- dat[[s]]
    cut10 <- quantile(d$qtime[!duplicated(d$id)], 0.10)
    vv <- list(`all respondents` = d,
      `Anders's protest rule removed` = d[!grepl(pat_anders, d$Q27), ],
      `protest rule without protected-value reasons` = d[!grepl(pat_noval, d$Q27), ],
      `no always-status-quo respondents` = d[d$always_sq == 0, ],
      `only believers in consequences (Q41)` = d[d$conseq == 1, ],
      `without 10% fastest` = d[d$qtime > cut10, ],
      `without 'difficult' (6-7)` = d[is.na(d$difficult) | d$difficult < 6, ])
    out <- lapply(names(vv), function(v) {
      m <- fit_model(paste0("rb_", s, "_", abs(sum(utf8ToInt(v)))), vv[[v]], prob_b, st0)
      e <- m$estimate
      se <- rse(m)
      sgn <- c(100, -100, -100, -100)
      names(sgn) <- key
      data.frame(variant = v, respondents = length(unique(vv[[v]]$id)),
        status_quo = sprintf("%6.0f (%4.0f)", sgn[["w_sq"]] * e[["w_sq"]], 100 * se[["w_sq"]]),
        no_protection = sprintf("%6.0f (%4.0f)", sgn[["w_nop"]] * e[["w_nop"]], 100 * se[["w_nop"]]),
        lt_1km = sprintf("%6.0f (%4.0f)", sgn[["w_d1"]] * e[["w_d1"]], 100 * se[["w_d1"]]),
        visible = sprintf("%6.0f (%4.0f)", sgn[["w_view"]] * e[["w_view"]], 100 * se[["w_view"]]))
    })
    rob[[s]] <- do.call(rbind, out)
    cat("\n=====", toupper(s), " WTA in NOK per household-year (robust s.e.) =====\n")
    print(rob[[s]], row.names = FALSE)
  }
  cat("\nShare of respondents removed by the two protest rules (Anders's / without protected-value reasons):\n")
  print(round(sapply(c("solar", "wind"), function(s) {
    q <- dat[[s]][!duplicated(dat[[s]]$id), "Q27"]
    c(grepl(pat_anders, q) |> mean(), grepl(pat_noval, q) |> mean())
  }), 3))
  cat("Share of the always-status-quo respondents whose reason is the protected-value statement 'Supplying Europe ... at the expense of Norwegian nature':\n")
  print(round(sapply(c("solar", "wind"), function(s) {
    r <- dat[[s]][!duplicated(dat[[s]]$id), ]
    r <- r[r$always_sq == 1, ]
    mean(grepl("bekostning av norsk", r$Q27))
  }), 3))
})


################################################################################
# 11, 15, 17  Solar vs wind in ONE mixed logit: formal likelihood-ratio tests of scale, taste heterogeneity and means  (step 15: scaled optimiser; step 17: another draw seed)
################################################################################
section(c("11", "15", "17"), {
  # ---- options of the pooled mixed logit
  # Steps 11, 15 and 17 are the same model and the same code; they differ in two options:
  #   step 15: the optimiser works on scaled parameters (factor 10 for every parameter whose start value exceeds 5 in absolute value, 1 for the others).
  #            Model, draws and likelihood are unchanged; the returned estimates are in the original units.  Jobs: U3 | scale3 | cov3 | attr3 | sq3 | all3
  #   step 17: as step 15, with another seed for the simulation draws (2027; Apollo's default is 13).  Jobs: U4 | all4
  #   (jobs whose name ends in a digit restart from the best first-round solution 'g'; the unscaled restarts of step 11 are U2 ... all2)
  SCALED <- STEP %in% c("15", "17")
  SEED <- if (STEP == "17") 2027 else NULL
  SCf <- function(b) setNames(ifelse(abs(b) > 5, 10, 1), names(b))
  scaling_option <- function(b, fixed) if (SCALED) list(scaling = SCf(b[setdiff(names(b), fixed)])) else list()

  # Article: "further analysis can include testing whether the differences in preferences between the technologies are statistically significant".
  # Leeds [W] 4B (pooled estimation, scale), [A] 1E (taste vs scale heterogeneity), [A] 1B (starting values; nested tests).
  # The two samples are different respondents, so they are STACKED into one likelihood (id is unique).  The pooled model is parameterised as
  #   solar value = theta,      wind value = theta + x          (x = wind minus solar difference)
  # for all 34 parameters of the full-covariance mixed logit of step 04 (5 correlated dimensions: status quo, no protection, < 1 km, visible,
  # log cost = SCALE).  With every x free this is exactly the two separate models, estimated on the SAME draws; a hypothesis is imposed by fixing
  # a group of x at zero, so every test is a likelihood-ratio test between nested models (no reliance on the inflated random-coefficient t-ratios).
  #   scale : x for the mean and spread of the log-scale and its 4 correlations with the tastes                 (6 df)
  #   cov   : x for the Cholesky terms of the 4 taste dimensions = equal taste covariance                         (10 df)
  #   attr  : x for the means of no protection, < 1 km, visible and the 7 fixed WTAs                              (10 df)
  #   sq    : x for the mean status quo constant                                                                  (1 df)
  #   g     : x for the 7 shifters of the status quo                                                              (7 df)
  #   all   : every x = 0, i.e. wind and solar respondents share ONE preference distribution                      (34 df)
  # Job name = the 2nd argument:  U (unrestricted) | scale | cov | attr | sq | g | all
  job <- SAMPLE
  ld1 <- function(f) {
    e <- new.env()
    load(file.path(V2, f), envir = e)
    e
  }
  ms <- ld1("state_04_solar.RData")$mM
  mw <- ld1("state_04_wind.RData")$mM
  pn <- names(ms$estimate)
  stopifnot(setequal(pn, names(mw$estimate)), length(pn) == 34)
  es <- ms$estimate
  ew <- mw$estimate[pn]
  chol_cols <- list(c("l11", "l21", "l31", "l41", "l51"), c("l22", "l32", "l42", "l52"), c("l33", "l43", "l53"), c("l44", "l54"), "l55")
  norm_chol <- function(e) {
    for (cl in chol_cols) if (e[[cl[1]]] < 0) e[cl] <- -e[cl]
    e
  }   # the sign of a Cholesky column is arbitrary
  es <- norm_chol(es)
  ew <- norm_chol(ew)
  G <- list(scale = c("m_cost", "l55", "l51", "l52", "l53", "l54"),
    cov   = c("l11", "l21", "l22", "l31", "l32", "l33", "l41", "l42", "l43", "l44"),
    attr  = c("m_nop", "m_d1", "m_view", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_wil", "w_d2", "w_d3"),
    sq    = "m_sq", g = c("g_age", "g_fem", "g_univ", "g_inc", "g_incm", "g_exp", "g_org"))
  stopifnot(setequal(unlist(G), pn), !anyDuplicated(unlist(G)))
  # Jobs whose name ends in "2" are RESTARTS from the best solution found in the first round (job 'g', LL -4557.14).  The first-round models U, scale, attr and cov
  # stopped at LL -4568 to -4571, i.e. BELOW the 'g' model that is nested in U: they had stopped early at a poor point and their LR tests are NOT valid.
  restart <- grepl("[0-9]$", job)
  base <- sub("[0-9]$", "", job)
  restr <- switch(base, U = character(0), all = pn, unlist(G[[base]]))
  stopifnot(base %in% c("U", "all", names(G)))
  if (!restart) {
    start_th <- ifelse(pn %in% restr, (es[pn] + ew[pn]) / 2, es[pn])        # restricted: start between the two separate solutions
    start_x  <- ifelse(pn %in% restr, 0, ew[pn] - es[pn])
    st <- c(start_th, setNames(start_x, paste0("x_", pn)))
    names(st)[1:34] <- pn
  } else {
    gbest <- ld1("state_11_g.RData")$mP$estimate
    st <- gbest[c(pn, paste0("x_", pn))]
    for (p in restr) {
      st[[p]] <- gbest[[p]] + gbest[[paste0("x_", p)]] / 2
      st[[paste0("x_", p)]] <- 0
    }   # restricted: common value = midpoint of solar and wind
  }
  fixed <- if (length(restr)) paste0("x_", restr) else character(0)
  stopifnot(!anyDuplicated(dat$pooled$id[dat$pooled$task == 1]), all(tapply(dat$pooled$wind, dat$pooled$id, function(v) length(unique(v))) == 1))

  rc_P <- function(apollo_beta, apollo_inputs) {
    mk <- function(th) list(
      sq   = th[["m_sq"]] + th[["l11"]] * xi1,
      nop  = th[["m_nop"]] + th[["l21"]] * xi1 + th[["l22"]] * xi2,
      d1   = th[["m_d1"]] + th[["l31"]] * xi1 + th[["l32"]] * xi2 + th[["l33"]] * xi3,
      view = th[["m_view"]] + th[["l41"]] * xi1 + th[["l42"]] * xi2 + th[["l43"]] * xi3 + th[["l44"]] * xi4,
      cost = exp(th[["m_cost"]] + th[["l51"]] * xi1 + th[["l52"]] * xi2 + th[["l53"]] * xi3 + th[["l54"]] * xi4 + th[["l55"]] * xi5))
    E <- environment()
    gv <- function(n) get(n, envir = E)
    # Apollo attaches the parameters (apollo_beta[...] is not allowed); its internal probe calls this function with dummy parameters, hence tryCatch
    tryCatch(
      {
        ths <- sapply(pn, gv)
        thw <- ths + sapply(paste0("x_", pn), gv)
        a <- mk(ths)
        b <- mk(thw)
        r <- list()
        for (n in names(a)) {
          r[[paste0("rs_", n)]] <- a[[n]]
          r[[paste0("rw_", n)]] <- b[[n]]
        }
        r
      },
      error = function(e) list())
  }
  prob_P <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    sel <- function(a, b) (1 - wind) * a + wind * b
    fx  <- function(n) sel(get(n), get(n) + get(paste0("x_", n)))
    bc <- sel(rs_cost, rw_cost)
    nop <- sel(rs_nop, rw_nop)
    d1 <- sel(rs_d1, rw_d1)
    vw <- sel(rs_view, rw_view)
    sq <- sel(rs_sq, rw_sq) + fx("g_age") * age_c + fx("g_fem") * female + fx("g_univ") * univ + fx("g_inc") * inc_c + fx("g_incm") * inc_missing + fx("g_exp") * exposed + fx("g_org") * org
    V <- list()
    V[["dev"]] <- bc * (fx("w_el38") * (el_alt1 == 38) + fx("w_el57") * (el_alt1 == 57) + fx("w_pl3") * (plants_alt1 == 3) + fx("w_pl6") * (plants_alt1 == 6) +
      nop * (wild_alt1 == 0) + fx("w_wil") * (wild_alt1 == 1) + d1 * (dist_alt1 == 1) + fx("w_d2") * (dist_alt1 == 2) +
      fx("w_d3") * (dist_alt1 == 3) + vw * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- bc * sq
    P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  dr_P <- list(interDrawsType = "mlhs", interNDraws = QK(250, 10), interNormDraws = paste0("xi", 1:5))
  t0 <- Sys.time()
  mP <- fit_model(paste0("pooledmix_", job), dat$pooled, prob_P, st, fixed = fixed, draws = dr_P, randcoeff = rc_P,
    estimate_settings = c(list(maxIterations = 400), scaling_option(st, fixed)), seed = SEED)
  if (restart) {   # convergence check: restart once from the solution found; keep the better of the two
    ll1 <- mP$maximum
    mP_b <- fit_model(paste0("pooledmix_", job, "b"), dat$pooled, prob_P, mP$estimate, fixed = fixed, draws = dr_P, randcoeff = rc_P,
      estimate_settings = c(list(maxIterations = 400), scaling_option(mP$estimate, fixed)), seed = SEED)
    cat(sprintf("First fit LL %.2f; restart from it LL %.2f\n", ll1, mP_b$maximum))
    if (mP_b$maximum > mP$maximum) mP <- mP_b
  }
  cat(sprintf("Pooled mixed logit, job '%s': LL %.2f, %d free parameters, BIC %.1f, %d observations (%.0f min; %s)\n", job, mP$maximum, nfree(mP), bic_of(mP), mP$nObs,
    as.numeric(difftime(Sys.time(), t0, units = "mins")), mP$message))
  cat(sprintf("Reference: the two separate models of step 04 sum to LL %.2f with %d parameters (different draws, so only roughly comparable)\n", ms$maximum + mw$maximum, nfree(ms) + nfree(mw)))
  if (base == "U") {
    e <- mP$estimate
    sr <- rse(mP)
    xn <- paste0("x_", pn)
    cat("\nWind minus solar difference of every parameter (x), robust t  [coefficients in 100 NOK / Cholesky units]\n")
    print(data.frame(parameter = pn, solar = round(e[pn], 3), wind_minus_solar = round(e[xn], 3), robust_t = round(e[xn] / sr[xn], 1), row.names = NULL), row.names = FALSE)
  }
})


################################################################################
# 12  Pooled MNL in WTA space: classic scale test (Swait-Louviere style) before the mixed logit of step 11
################################################################################
section("12", {
  # Leeds [W] 4B.  Wind utility is multiplied by a scale lambda (solar = 1).  In WTA space lambda is the ratio of the cost coefficients.
  #   U  all WTAs and the status quo constant technology-specific, lambda free      (24 parameters)
  #   A  one set of WTAs for both technologies, lambda free                          -> U vs A: do the WTAs differ, once scale is allowed to differ?
  #   B  like A, lambda = 1                                                          -> A vs B: is the scale different?
  #   C  only the status quo constant technology-specific, lambda free               -> A vs C: status quo differs?  C vs U: do the ten other WTAs differ?
  #   D  like C, lambda = 1                                                          -> C vs D: scale difference once the status quo may differ
  wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view")
  ld1 <- function(f) {
    e <- new.env()
    load(file.path(V2, f), envir = e)
    e
  }
  l2 <- ld1("state_02_all.RData")$lad
  bs <- l2$solar$m0$estimate
  bw <- l2$wind$m0$estimate
  prob_S <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))
    sc <- 1 + (lambda - 1) * wind
    fx <- function(n) get(n) + get(paste0("x_", n)) * wind
    V <- list()
    V[["dev"]] <- sc * b_cost * (fx("w_el38") * (el_alt1 == 38) + fx("w_el57") * (el_alt1 == 57) + fx("w_pl3") * (plants_alt1 == 3) + fx("w_pl6") * (plants_alt1 == 6) +
      fx("w_nop") * (wild_alt1 == 0) + fx("w_wil") * (wild_alt1 == 1) + fx("w_d1") * (dist_alt1 == 1) + fx("w_d2") * (dist_alt1 == 2) +
      fx("w_d3") * (dist_alt1 == 3) + fx("w_view") * (view_alt1 == 2) + cost_alt1 / 100)
    V[["sq"]] <- sc * b_cost * fx("w_sq")
    P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_prepareProb(P, apollo_inputs, functionality)
    return(P)
  }
  sol <- c(bs[c(wn, "b_cost")])
  st <- c(sol, setNames(bw[wn] - bs[wn], paste0("x_", wn)), lambda = unname(bw[["b_cost"]] / bs[["b_cost"]]))
  xs <- paste0("x_", wn)
  runS <- function(nm, fixed, start) fit_model(nm, dat$pooled, prob_S, start, fixed = fixed)
  st0 <- function(zero) {
    s <- st
    s[zero] <- 0
    s
  }
  mU <- runS("mnl_U", c(), st)
  mA <- runS("mnl_A", xs, st0(xs))
  mB <- runS("mnl_B", c(xs, "lambda"), replace(st0(xs), "lambda", 1))
  mC <- runS("mnl_C", xs[-1], st0(xs[-1]))
  mD <- runS("mnl_D", c(xs[-1], "lambda"), replace(st0(xs[-1]), "lambda", 1))
  mods <- list(U = mU, A = mA, B = mB, C = mC, D = mD)
  nm <- c("U  everything technology-specific, scale free", "A  common WTAs, scale free", "B  common WTAs, scale = 1", "C  only status quo differs, scale free", "D  only status quo differs, scale = 1")
  print(data.frame(model = nm, LL = round(sapply(mods, function(m) m$maximum), 1), k = sapply(mods, nfree), BIC = round(sapply(mods, bic_of), 1)), row.names = FALSE)
  lr <- function(big, small) {
    x <- 2 * (big$maximum - small$maximum)
    d <- nfree(big) - nfree(small)
    data.frame(LR = round(x, 1), df = d, p = signif(pchisq(x, d, lower.tail = FALSE), 2))
  }
  cat("\nLikelihood-ratio tests\n")
  print(rbind("U vs A  Do the WTAs differ (scale allowed to differ)?" = lr(mU, mA), "A vs B  Is the scale different (WTAs common)?" = lr(mA, mB),
    "A vs C  Does the status quo constant differ (scale free)?" = lr(mC, mA), "C vs U  Do the ten other WTAs differ (SQ and scale free)?" = lr(mU, mC),
    "C vs D  Is the scale different (SQ may differ)?" = lr(mC, mD)))
  e <- mC$estimate
  se <- rse(mC)
  cat(sprintf("\nModel C: wind scale lambda = %.2f (robust s.e. %.2f; test of lambda = 1: z = %.2f);  status quo WTA difference wind - solar = %.0f NOK (z %.1f)\n",
    e[["lambda"]], se[["lambda"]], (e[["lambda"]] - 1) / se[["lambda"]], 100 * e[["x_w_sq"]], e[["x_w_sq"]] / se[["x_w_sq"]]))
  eU <- mU$estimate
  sU <- rse(mU)
  cat("\nModel U: WTA differences wind - solar (NOK; WTA = -100 x coefficient, status quo +100 x), robust t\n")
  print(data.frame(attribute = c("status quo", sub("w_", "", wn[-1])), solar = round(100 * ifelse(wn == "w_sq", 1, -1) * eU[wn]),
    wind_minus_solar = round(100 * ifelse(wn == "w_sq", 1, -1) * eU[xs]), t = round(eU[xs] / sU[xs], 1), row.names = NULL), row.names = FALSE)
})


################################################################################
# 13  Summary of the pooled mixed logit tests (step 11): likelihood-ratio tests, Holm correction, and the differences in WTA with delta-method s.e.
################################################################################
section("13", {
  # Round 2 (restarts from the best first-round solution, jobs ending in 2); 'g' of round 1 is the starting solution and is also used as it is.
  fl <- c(U = "U2", scale = "scale2", cov = "cov2", attr = "attr2", sq = "sq2", g = "g", all = "all2")
  jobs <- names(fl)
  # The unrestricted model is U3 (step 15: same model and draws as U2, optimiser with parameter scaling); it reached a higher likelihood than U2 (-4552.17 vs -4552.42),
  # and U2 is below the restricted 'scale2' (-4552.39), which a nested model cannot be. U3 is therefore the reference for all tests.
  M <- lapply(setNames(jobs, jobs), function(j) {
    e <- new.env()
    f <- if (j == "U") "state_15_U3.RData" else sprintf("state_11_%s.RData", fl[[j]])
    load(file.path(V2, f), envir = e)
    e$mP
  })
  cat("Unrestricted model must have the highest likelihood; check:", all(sapply(M, function(m) m$maximum) <= M$U$maximum + 1e-6), "\n")
  lab <- c(U = "U  all 34 parameters technology-specific", scale = "scale: mean, spread, correlations of the log-scale equal (6)", cov = "cov: covariance of the 4 taste dimensions equal (10)",
    attr = "attr: means of 10 attribute WTAs equal (10)", sq = "sq: mean status quo constant equal (1)", g = "g: 7 shifters of the status quo equal (7)", all = "all: one preference distribution for both technologies (34)")
  print(data.frame(model = lab[jobs], LL = round(sapply(M, function(m) m$maximum), 1), k = sapply(M, nfree), BIC = round(sapply(M, bic_of), 1),
    converged = sapply(M, function(m) m$message), row.names = NULL), row.names = FALSE)
  LLU <- M$U$maximum
  kU <- nfree(M$U)
  tst <- do.call(rbind, lapply(jobs[-1], function(j) {
    x <- 2 * (LLU - M[[j]]$maximum)
    d <- kU - nfree(M[[j]])
    data.frame(restriction = lab[[j]], LR = round(x, 1), df = d, p = signif(pchisq(x, d, lower.tail = FALSE), 2))
  }))
  ps <- tst$p[tst$restriction != lab[["all"]]]
  tst$p_Holm <- NA
  tst$p_Holm[tst$restriction != lab[["all"]]] <- signif(p.adjust(ps, "holm"), 2)
  cat("\nLikelihood-ratio tests against the unrestricted pooled model (= solar and wind separate, same draws)\n")
  print(tst, row.names = FALSE)
  cat("Holm adjustment covers the five single-block tests; 'all' is the joint test.\n")

  # ---- differences in WTA with delta-method s.e. (robust covariance of the unrestricted pooled model)
  e <- M$U$estimate
  V <- M$U$robvarcov
  nv <- rownames(V)
  dm <- function(f) {
    J <- sapply(nv, function(n) {
      h <- 1e-5 * max(1, abs(e[[n]]))
      tp <- e
      tp[n] <- e[[n]] + h
      tn <- e
      tn[n] <- e[[n]] - h
      (f(tp) - f(tn)) / (2 * h)
    })
    c(est = f(e), se = sqrt(drop(J %*% V %*% J)))
  }
  at <- list(`status quo (reference respondent)` = c("m_sq", +1), `no protection` = c("m_nop", -1), `< 1 km` = c("m_d1", -1), `visible` = c("m_view", -1),
    `production 38k` = c("w_el38", -1), `production 57k` = c("w_el57", -1), `plants moderate` = c("w_pl3", -1), `plants many` = c("w_pl6", -1),
    `wilderness only` = c("w_wil", -1), `1-3 km` = c("w_d2", -1), `3-6 km` = c("w_d3", -1))
  tb <- do.call(rbind, lapply(names(at), function(n) {
    p <- at[[n]][1]
    s <- as.numeric(at[[n]][2])
    xp <- paste0("x_", p)
    sol <- dm(function(th) 100 * s * th[[p]])
    dif <- dm(function(th) 100 * s * th[[xp]])
    data.frame(attribute = n, solar = round(sol["est"]), wind = round(sol["est"] + dif["est"]), difference = round(dif["est"]), se = round(dif["se"]), z = round(dif["est"] / dif["se"], 2), row.names = NULL)
  }))
  cat("\nWTA in NOK per household-year, solar vs wind from the unrestricted pooled model (difference = wind - solar; delta-method s.e.)\n")
  print(tb, row.names = FALSE)
  r <- dm(function(th) (th[["m_sq"]] + th[["x_m_sq"]]) / th[["m_sq"]])
  r35 <- dm(function(th) (th[["m_sq"]] + th[["x_m_sq"]]) - 3.5 * th[["m_sq"]])
  cat(sprintf("\nStatus quo constant, wind / solar (reference respondent): %.2f (delta s.e. %.2f); test of the article's '3.5 times': wind - 3.5 x solar = %.0f NOK (s.e. %.0f, z %.2f)\n",
    r["est"], r["se"], 100 * r35["est"], 100 * r35["se"], r35["est"] / r35["se"]))
  # ---- scale
  sc <- dm(function(th) th[["x_m_cost"]])
  sd5 <- dm(function(th) th[["l55"]] + th[["x_l55"]])
  sd5s <- dm(function(th) th[["l55"]])
  cat(sprintf("\nScale (log of the cost coefficient): solar mean %.2f, wind mean %.2f; wind - solar = %.2f (delta s.e. %.2f, z %.2f) i.e. a factor %.2f on the cost coefficient / scale\n",
    e[["m_cost"]], e[["m_cost"]] + e[["x_m_cost"]], sc["est"], sc["se"], sc["est"] / sc["se"], exp(sc["est"])))
  cat(sprintf("Spread of the log-scale (Cholesky l55): solar %.2f, wind %.2f\n", abs(sd5s["est"]), abs(sd5["est"])))
  cat("Cost-scale correlations (solar | wind) are in the output of step 04; here: the 'scale' LR test above covers mean, spread and correlations jointly.\n")
})


################################################################################
# 14  Figure: acceptance curves (share of households accepting a scenario) against annual bill reduction
################################################################################
section("14", {
  # Same model and delta-method standard errors as step 09 part C (LC2-D, population-weighted).  The code of step 09 is sourced silently so
  # that the function 'accept' and the delta-method function 'dm' are identical to the ones used for the tables.
  invisible(capture.output(withAutoprint(POLICY, evaluated = TRUE, local = FALSE, echo = FALSE)))
  cvg <- seq(0, 50000, by = 2500)
  sc <- list(P6 = "only production 38k, no visibility (cheapest scenario)", P1 = "38k, < 1 km from home", P4 = "38k, many plants, no protection")
  cols <- c(P6 = "#1f6f8b", P1 = "#d9822b", P4 = "#7a7a7a")
  png(file.path(V2, "output", "fig_acceptance.png"), width = 2000, height = 900, res = 220)
  par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 0.8), mgp = c(2.6, 0.7, 0), las = 1)
  for (s in c("solar", "wind")) {
    m <- lc[[s]]$lD
    plot(NA, xlim = c(0, 50), ylim = c(0, 60), xlab = "Annual bill reduction (1,000 NOK per household)", ylab = "Share accepting the project (%)", main = toupper(s), cex.main = 1)
    abline(v = 9, lty = 3, col = "grey60")
    text(9.6, 58, "highest amount on the cards (9,000 NOK)", adj = 0, cex = 0.6, col = "grey40")
    for (p in names(sc)) {
      res <- sapply(cvg, function(cc) dm(m, function(th) accept(th, s, scen[[p]], FALSE, cc, TRUE)))
      est <- 100 * res["est", ]
      se <- 100 * res["se", ]
      polygon(c(cvg / 1000, rev(cvg / 1000)), c(pmax(est - 1.96 * se, 0), rev(est + 1.96 * se)), col = adjustcolor(cols[[p]], 0.18), border = NA)
      lines(cvg / 1000, est, col = cols[[p]], lwd = 2.2)
    }
    if (s == "solar") legend("bottomright", legend = unname(sc), col = cols, lwd = 2.2, bty = "n", cex = 0.6)
  }
  dev.off()
  cat("figure written to", file.path(V2, "output", "fig_acceptance.png"), "\n")
})


################################################################################
# 16  Collinearity check on the covariates used as status quo / class shifters (step 01, 02, 06, 11)
################################################################################
section("16", {
  # Variables: age_c, female, univ, inc_c, inc_missing, org, exposed, att_c, conc_c, natv_c, natp_c, conseq
  # (the same set as the glm in 01_explore.R and the g_ / gd_ shifters in 02, 06, 11).
  # Pairwise correlations (like a textbook logistic-regression collinearity check) plus VIFs from a linear
  # probe of each variable on all the others (1 / (1 - R^2)); no external package needed.
  vars <- c("age_c", "female", "univ", "inc_c", "inc_missing", "org", "exposed", "att_c", "conc_c", "natv_c", "natp_c", "conseq")
  lab  <- c("age (centred)", "women", "university degree", "income (centred)", "income missing", "env. organisation", "exposed (near/owns land)",
    "attitude Q7", "concern Q10", "visit nature Q31", "preserve nature Q32", "believes results used")

  resp <- lapply(dat[c("solar", "wind")], function(d) d[!duplicated(d$id), vars])
  resp$pooled <- rbind(resp$solar, resp$wind)

  vif <- function(d) {
    sapply(vars, function(v) {
      f <- as.formula(paste(v, "~", paste(setdiff(vars, v), collapse = " + ")))
      r2 <- summary(lm(f, d))$r.squared
      1 / (1 - r2)
    })
  }

  for (s in c("solar", "wind", "pooled")) {
    d <- resp[[s]]
    cat("\n===== ", toupper(s), " (n =", nrow(d), ") =====\n")
    cm <- cor(d, use = "pairwise.complete.obs")
    dimnames(cm) <- list(lab, lab)
    cat("\nPairwise correlation matrix\n")
    print(round(cm, 2))
    off <- cm
    diag(off) <- NA
    hi <- which(abs(off) >= 0.5, arr.ind = TRUE)
    if (nrow(hi)) {
      cat("\nPairs with |r| >= 0.5:\n")
      print(data.frame(var1 = lab[hi[, 1]], var2 = lab[hi[, 2]], r = round(cm[hi], 2))[hi[, 1] < hi[, 2], ], row.names = FALSE)
    } else cat("\nNo pair reaches |r| >= 0.5\n")
    cat("\nVariance inflation factors (rule of thumb: concern from VIF > 5, serious from VIF > 10)\n")
    print(data.frame(variable = lab, VIF = round(vif(d), 2)), row.names = FALSE)
  }
})


################################################################################
# SAVE THE RESULTS OF THE STEP (later steps load them)
################################################################################
.new <- setdiff(ls(all.names = TRUE), c(.before, ".before", ".new", ".out", ".t0"))
.out <- file.path(V2, paste0(if (QUICK) "quick_" else "", if (STEP == "00") "state_base.RData" else sprintf("state_%s_%s.RData", STEP, SAMPLE)))
save(list = .new, file = .out)
cat(sprintf("\n[step %s (%s) finished in %.1f min]\n", STEP, SAMPLE, as.numeric(difftime(Sys.time(), .t0, units = "mins"))))
