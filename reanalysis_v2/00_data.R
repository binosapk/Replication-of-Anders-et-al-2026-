# =====================================================================
# 00  Data preparation: choice data, respondent covariates, population weights
# =====================================================================
suppressPackageStartupMessages({ library(apollo); library(dplyr); library(tidyr) })
# Paths come from run_step.R (.V2 = this folder, .ROOT = its parent with the data csv files); if sourced directly, the working directory is taken as this folder.
V2     <- if (exists(".V2")) .V2 else normalizePath(getwd(), winslash = "/")
ROOT   <- dirname(V2)
setwd(ROOT); set.seed(2026)
outdir <- file.path(V2, "apollo_output"); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------- helpers
quiet <- function(expr) { out <- NULL; invisible(capture.output(out <- suppressMessages(expr))); out }
num1  <- function(x) suppressWarnings(as.numeric(substr(x, 1, 1)))        # "4 Neutral" -> 4
nfree <- function(m) length(m$estimate) - length(m$apollo_fixed)
bic_of <- function(m) -2 * m$maximum + nfree(m) * log(m$nObs)
aic_of <- function(m) -2 * m$maximum + 2 * nfree(m)
rse   <- function(m) { v <- sqrt(diag(m$robvarcov)); names(v) <- rownames(m$robvarcov); v }   # robust s.e. with names
cse   <- function(m) { v <- sqrt(diag(m$varcov));    names(v) <- rownames(m$varcov);    v }   # classical s.e.
lrp   <- function(big, small) signif(pchisq(2 * (big$maximum - small$maximum), nfree(big) - nfree(small), lower.tail = FALSE), 2)

# one wrapper that estimates any Apollo model silently
fit_model <- function(name, data, probs, beta, fixed = c(), draws = NULL, randcoeff = NULL, lcpars = NULL,
                      ncores = 1, estimate_settings = list()) {
  database <<- data
  quiet(apollo_initialise())
  apollo_control <<- list(modelName = name, indivID = "id", mixing = !is.null(draws), nCores = ncores,
                          noDiagnostics = TRUE, outputDirectory = outdir)
  apollo_beta <<- beta; apollo_fixed <<- fixed
  if (!is.null(draws))     apollo_draws     <<- draws
  if (!is.null(randcoeff)) apollo_randCoeff <<- randcoeff
  if (!is.null(lcpars))    apollo_lcPars    <<- lcpars
  on.exit(suppressWarnings(rm(apollo_draws, apollo_randCoeff, apollo_lcPars, envir = globalenv())))
  apollo_inputs <<- apollo_validateInputs(silent = TRUE)
  es <- modifyList(list(silent = TRUE, writeIter = FALSE), estimate_settings)
  quiet(apollo_estimate(beta, fixed, probs, apollo_inputs, estimate_settings = es))
}

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
    if (grepl("^Mer", s))    return(4.5)
    if (length(nums) >= 2)   return(mean(nums[1:2]) / 1e6)
    NA_real_ })
}

# ---------------------------------------------------------------- one survey file -> one row per answered card
prepare <- function(file, tech, id_start) {
  r <- read.csv(file, encoding = "UTF-8", stringsAsFactors = FALSE)
  r$id <- id_start + seq_len(nrow(r)); r$tech <- tech; r$wind <- as.numeric(tech == "wind")
  # --- demographics (exogenous)
  r$female <- as.numeric(r$gender == "Kvinne")
  r$age_c  <- (r$age - 50) / 10
  r$agegrp <- cut(r$age, c(17, 29, 39, 49, 59, 120), labels = c("18-29", "30-39", "40-49", "50-59", "60+"))
  r$educ4  <- ifelse(grepl("^Grunnskole", r$Q38), 1,
              ifelse(grepl("^Videreg|^Fagbrev", r$Q38), 2,
              ifelse(grepl("^3-4", r$Q38), 3,
              ifelse(grepl("^5-|^PhD", r$Q38), 4, NA))))
  r$educ_missing <- as.numeric(is.na(r$educ4)); r$educ4[is.na(r$educ4)] <- 2     # 1% "unsure": modal category
  r$univ   <- as.numeric(r$educ4 >= 3)
  inc <- income_mid(r$Q40); r$inc_missing <- as.numeric(is.na(inc))
  r$inc_c  <- ifelse(is.na(inc), 0, inc - 0.8)                                    # household income, million NOK, centred at 0.8
  r$org    <- as.numeric(r$Q37 == "Ja")
  r$near   <- as.numeric(grepl("bor i n", r$Q39)); r$owner <- as.numeric(grepl("eier arealer", r$Q39))
  r$exposed <- as.numeric(r$near == 1 | r$owner == 1)                             # lives near / owns land
  r$sees   <- as.numeric(r$Q13 == "Ja")
  r$own_pv <- as.numeric(r$Q36 == "Ja")
  # --- attitudes and values (potentially endogenous; used as indicators or as diagnostics)
  r$att    <- num1(r$Q7);  r$conc <- num1(r$Q10);  r$nat_visit <- num1(r$Q31)
  r$nat_pres <- ifelse(grepl("ikke redusere", r$Q32), 3, ifelse(grepl("noe reduksjon", r$Q32), 2, ifelse(grepl("^Nei", r$Q32), 1, NA)))
  r$impact_own <- if (tech == "wind") num1(r$Q33r1) else num1(r$Q33r2)           # perceived landscape impact of own technology
  r$imp_el <- num1(r$Q25r1); r$imp_pl <- num1(r$Q25r2); r$imp_wild <- num1(r$Q25r3)
  r$imp_dist <- num1(r$Q25r4); r$imp_view <- num1(r$Q25r5); r$imp_cost <- num1(r$Q25r6)
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
pop_cell <- aggregate(n ~ agegrp + female, pop, sum); pop_cell$share <- pop_cell$n / sum(pop_cell$n)
edu_pop <- c(`1` = 83586, `2` = 127937 + 11217, `3` = 72522, `4` = 24563)          # basic; upper secondary (+ tertiary vocational); higher short; higher long
edu_share <- edu_pop / sum(edu_pop)

rake <- function(r, iter = 100) {            # r: one row per respondent
  w <- rep(1, nrow(r))
  cell <- paste(r$agegrp, r$female); pcell <- setNames(pop_cell$share, paste(pop_cell$agegrp, pop_cell$female))
  for (k in 1:iter) {
    s <- tapply(w, cell, sum) / sum(w);            w <- w * (pcell[cell] / s[cell])
    s <- tapply(w, r$educ4, sum) / sum(w);         w <- w * (edu_share[as.character(r$educ4)] / s[as.character(r$educ4)])
  }
  w / mean(w)
}
for (s in c("solar", "wind")) {
  d <- dat[[s]]; r <- d[!duplicated(d$id), ]
  w <- setNames(rake(r), r$id); d$wt <- unname(w[as.character(d$id)]); dat[[s]] <- d
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
  a <- prop.table(table(r$agegrp)); e <- prop.table(table(r$educ4))
  f <- mean(r$female)
  data.frame(variable = c(paste("age", levels(r$agegrp)), "women", "basic school", "upper secondary", "higher short (<=4y)", "higher long (5y+)"),
             sample = round(100 * c(a, f, e), 1))
}
pa <- c(tapply(pop$n, pop$agegrp, sum) / sum(pop$n), women = sum(pop$n[pop$female == 1]) / sum(pop$n), edu_share)
tab <- cbind(cmp("solar"), wind = cmp("wind")$sample, population = round(100 * pa, 1))
names(tab)[2] <- "solar"; print(tab, row.names = FALSE)
cat("\nWeights after raking on age x sex and education: summary and effective sample size\n")
for (s in c("solar", "wind")) { w <- dat[[s]]$wt[!duplicated(dat[[s]]$id)]
  cat(sprintf("%-6s min %.2f  max %.2f  effective n = %.0f of %d\n", s, min(w), max(w), sum(w)^2 / sum(w^2), length(w))) }
cat("\nShare of respondents with missing information: income", round(mean(dat$pooled$inc_missing[!duplicated(dat$pooled$id)]), 3),
    " education", round(mean(dat$pooled$educ_missing[!duplicated(dat$pooled$id)]), 3), "\n")
