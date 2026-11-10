# =====================================================================
# 12  Pooled MNL in WTA space: classic scale test (Swait-Louviere style) before the mixed logit of step 11
# =====================================================================
# Leeds [W] 4B.  Wind utility is multiplied by a scale lambda (solar = 1).  In WTA space lambda is the ratio of the cost coefficients.
#   U  all WTAs and the status quo constant technology-specific, lambda free      (24 parameters)
#   A  one set of WTAs for both technologies, lambda free                          -> U vs A: do the WTAs differ, once scale is allowed to differ?
#   B  like A, lambda = 1                                                          -> A vs B: is the scale different?
#   C  only the status quo constant technology-specific, lambda free               -> A vs C: status quo differs?  C vs U: do the ten other WTAs differ?
#   D  like C, lambda = 1                                                          -> C vs D: scale difference once the status quo may differ
wn <- c("w_sq", "w_el38", "w_el57", "w_pl3", "w_pl6", "w_nop", "w_wil", "w_d1", "w_d2", "w_d3", "w_view")
ld1 <- function(f) { e <- new.env(); load(file.path(V2, f), envir = e); e }
l2 <- ld1("state_02_all.RData")$lad
bs <- l2$solar$m0$estimate; bw <- l2$wind$m0$estimate
prob_S <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  sc <- 1 + (lambda - 1) * wind
  fx <- function(n) get(n) + get(paste0("x_", n)) * wind
  V <- list()
  V[["dev"]] <- sc * b_cost * (fx("w_el38") * (el_alt1 == 38) + fx("w_el57") * (el_alt1 == 57) + fx("w_pl3") * (plants_alt1 == 3) + fx("w_pl6") * (plants_alt1 == 6) +
                               fx("w_nop") * (wild_alt1 == 0) + fx("w_wil") * (wild_alt1 == 1) + fx("w_d1") * (dist_alt1 == 1) + fx("w_d2") * (dist_alt1 == 2) +
                               fx("w_d3") * (dist_alt1 == 3) + fx("w_view") * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- sc * b_cost * fx("w_sq")
  P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
sol <- c(bs[c(wn, "b_cost")])
st <- c(sol, setNames(bw[wn] - bs[wn], paste0("x_", wn)), lambda = unname(bw[["b_cost"]] / bs[["b_cost"]]))
xs <- paste0("x_", wn)
runS <- function(nm, fixed, start) fit_model(nm, dat$pooled, prob_S, start, fixed = fixed)
st0 <- function(zero) { s <- st; s[zero] <- 0; s }
mU <- runS("mnl_U", c(), st)
mA <- runS("mnl_A", xs, st0(xs))
mB <- runS("mnl_B", c(xs, "lambda"), replace(st0(xs), "lambda", 1))
mC <- runS("mnl_C", xs[-1], st0(xs[-1]))
mD <- runS("mnl_D", c(xs[-1], "lambda"), replace(st0(xs[-1]), "lambda", 1))
mods <- list(U = mU, A = mA, B = mB, C = mC, D = mD)
nm <- c("U  everything technology-specific, scale free", "A  common WTAs, scale free", "B  common WTAs, scale = 1", "C  only status quo differs, scale free", "D  only status quo differs, scale = 1")
print(data.frame(model = nm, LL = round(sapply(mods, function(m) m$maximum), 1), k = sapply(mods, nfree), BIC = round(sapply(mods, bic_of), 1)), row.names = FALSE)
lr <- function(big, small) { x <- 2 * (big$maximum - small$maximum); d <- nfree(big) - nfree(small); data.frame(LR = round(x, 1), df = d, p = signif(pchisq(x, d, lower.tail = FALSE), 2)) }
cat("\nLikelihood-ratio tests\n")
print(rbind("U vs A  Do the WTAs differ (scale allowed to differ)?" = lr(mU, mA), "A vs B  Is the scale different (WTAs common)?" = lr(mA, mB),
            "A vs C  Does the status quo constant differ (scale free)?" = lr(mC, mA), "C vs U  Do the ten other WTAs differ (SQ and scale free)?" = lr(mU, mC),
            "C vs D  Is the scale different (SQ may differ)?" = lr(mC, mD)))
e <- mC$estimate; se <- rse(mC)
cat(sprintf("\nModel C: wind scale lambda = %.2f (robust s.e. %.2f; test of lambda = 1: z = %.2f);  status quo WTA difference wind - solar = %.0f NOK (z %.1f)\n",
            e[["lambda"]], se[["lambda"]], (e[["lambda"]] - 1) / se[["lambda"]], 100 * e[["x_w_sq"]], e[["x_w_sq"]] / se[["x_w_sq"]]))
eU <- mU$estimate; sU <- rse(mU)
cat("\nModel U: WTA differences wind - solar (NOK; WTA = -100 x coefficient, status quo +100 x), robust t\n")
print(data.frame(attribute = c("status quo", sub("w_", "", wn[-1])), solar = round(100 * ifelse(wn == "w_sq", 1, -1) * eU[wn]),
                 wind_minus_solar = round(100 * ifelse(wn == "w_sq", 1, -1) * eU[xs]), t = round(eU[xs] / sU[xs], 1), row.names = NULL), row.names = FALSE)
