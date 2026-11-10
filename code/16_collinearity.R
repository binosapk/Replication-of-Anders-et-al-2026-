# =====================================================================
# 16  Collinearity check on the covariates used as status quo / class shifters (step 01, 02, 06, 11)
# =====================================================================
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
  off <- cm; diag(off) <- NA
  hi <- which(abs(off) >= 0.5, arr.ind = TRUE)
  if (nrow(hi)) {
    cat("\nPairs with |r| >= 0.5:\n")
    print(data.frame(var1 = lab[hi[, 1]], var2 = lab[hi[, 2]], r = round(cm[hi], 2))[hi[, 1] < hi[, 2], ], row.names = FALSE)
  } else cat("\nNo pair reaches |r| >= 0.5\n")
  cat("\nVariance inflation factors (rule of thumb: concern from VIF > 5, serious from VIF > 10)\n")
  print(data.frame(variable = lab, VIF = round(vif(d), 2)), row.names = FALSE)
}
