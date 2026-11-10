# =====================================================================
# 01  Study the data first (Leeds [A] 1B: "start with the data, even before a simple MNL, analyse them non-parametrically")
# =====================================================================
# (a) acceptance by compensation level, with confidence intervals and population weights
# (b) who never trades? binary logit of "chose the status quo on all 8 cards"
# (c) does belief that the results will be used (consequentiality, Q41) matter?
wmean <- function(x, w) sum(x * w) / sum(w)
cat("(a) Share choosing the development option by bill reduction (%), unweighted [95% CI] and weighted\n")
acc <- do.call(rbind, lapply(c("solar", "wind"), function(s) {
  d <- dat[[s]]
  do.call(rbind, lapply(sort(unique(d$cost_alt1)), function(cst) {
    z <- d[d$cost_alt1 == cst, ]; p <- mean(z$choice == 1); n <- nrow(z)
    # cluster-robust s.e. (respondents answer several cards): use respondent-level means
    zi <- tapply(z$choice == 1, z$id, mean); se <- sd(zi) / sqrt(length(zi))
    data.frame(sample = s, cost = cst, n = n, share = round(100 * p, 1), lo = round(100 * (p - 1.96 * se), 1), hi = round(100 * (p + 1.96 * se), 1),
               weighted = round(100 * wmean(z$choice == 1, z$wt), 1)) }))}))
print(reshape(acc[, c("sample", "cost", "share", "lo", "hi", "weighted")], idvar = "cost", timevar = "sample", direction = "wide"), row.names = FALSE)

# slope test: logit of development choice on cost, cluster-robust by respondent (cards differ in other attributes, so this is only descriptive)
cat("\nDescriptive slope: change in P(development) per +1,000 NOK (linear probability model, respondent-clustered s.e.)\n")
for (s in c("solar", "wind")) {
  d <- dat[[s]]; d$y <- as.numeric(d$choice == 1); d$c1000 <- d$cost_alt1 / 1000
  f <- lm(y ~ c1000, d); X <- model.matrix(f); e <- resid(f)
  meat <- crossprod(rowsum(X * e, d$id)); bread <- solve(crossprod(X)); V <- bread %*% meat %*% bread
  cat(sprintf("%-6s slope %.4f (cluster s.e. %.4f), t = %.2f\n", s, coef(f)[2], sqrt(V[2, 2]), coef(f)[2] / sqrt(V[2, 2])))
}

cat("\n(b) Who NEVER trades?  Binary logit, odds ratios (z in brackets), respondent level\n")
resp <- lapply(c("solar", "wind"), function(s) { r <- dat[[s]][!duplicated(dat[[s]]$id), ]; r })
names(resp) <- c("solar", "wind")
form <- always_sq ~ I(age_c) + female + univ + inc_c + inc_missing + org + exposed + att_c + conc_c + natv_c + natp_c + conseq
lab <- c("age (+10 y)", "women", "university degree", "income (+1 M NOK)", "income missing", "env. organisation", "exposed (lives near / owns land)",
         "attitude Q7 (+1)", "concern Q10 (+1)", "visit nature Q31 (+1)", "preserve nature Q32 (+1)", "believes results are used")
fits <- lapply(resp, function(r) glm(form, binomial, r))
tab <- sapply(fits, function(f) { s <- summary(f)$coefficients[-1, ]; sprintf("%5.2f (%5.1f)", exp(s[, 1]), s[, 3]) })
print(data.frame(variable = lab, solar = tab[, 1], wind = tab[, 2]), row.names = FALSE)
cat(sprintf("McFadden pseudo-R2: solar %.3f, wind %.3f\n", 1 - fits$solar$deviance / fits$solar$null.deviance, 1 - fits$wind$deviance / fits$wind$null.deviance))

cat("\n(c) Consequentiality (Q41): 'very/fairly sure the authorities will use the results'\n")
print(do.call(rbind, lapply(c("solar", "wind"), function(s) { r <- resp[[s]]
  data.frame(sample = s, share_sure = round(mean(r$conseq), 3),
             always_SQ_if_sure = round(mean(r$always_sq[r$conseq == 1]), 3), always_SQ_if_unsure = round(mean(r$always_sq[r$conseq == 0]), 3)) })), row.names = FALSE)

cat("\nWeighted vs unweighted share of respondents who ALWAYS chose the status quo\n")
print(round(sapply(resp, function(r) c(unweighted = mean(r$always_sq), weighted = wmean(r$always_sq, r$wt))), 3))
