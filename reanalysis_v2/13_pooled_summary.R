# =====================================================================
# 13  Summary of the pooled mixed logit tests (step 11): likelihood-ratio tests, Holm correction, and the differences in WTA with delta-method s.e.
# =====================================================================
# Round 2 (restarts from the best first-round solution, jobs ending in 2); 'g' of round 1 is the starting solution and is also used as it is.
fl <- c(U = "U2", scale = "scale2", cov = "cov2", attr = "attr2", sq = "sq2", g = "g", all = "all2")
jobs <- names(fl)
# The unrestricted model is U3 (step 15: same model and draws as U2, optimiser with parameter scaling); it reached a higher likelihood than U2 (-4552.17 vs -4552.42),
# and U2 is below the restricted 'scale2' (-4552.39), which a nested model cannot be. U3 is therefore the reference for all tests.
M <- lapply(setNames(jobs, jobs), function(j) { e <- new.env(); f <- if (j == "U") "state_15_U3.RData" else sprintf("state_11_%s.RData", fl[[j]]); load(file.path(V2, f), envir = e); e$mP })
cat("Unrestricted model must have the highest likelihood; check:", all(sapply(M, function(m) m$maximum) <= M$U$maximum + 1e-6), "\n")
lab <- c(U = "U  all 34 parameters technology-specific", scale = "scale: mean, spread, correlations of the log-scale equal (6)", cov = "cov: covariance of the 4 taste dimensions equal (10)",
         attr = "attr: means of 10 attribute WTAs equal (10)", sq = "sq: mean status quo constant equal (1)", g = "g: 7 shifters of the status quo equal (7)", all = "all: one preference distribution for both technologies (34)")
print(data.frame(model = lab[jobs], LL = round(sapply(M, function(m) m$maximum), 1), k = sapply(M, nfree), BIC = round(sapply(M, bic_of), 1),
                 converged = sapply(M, function(m) m$message), row.names = NULL), row.names = FALSE)
LLU <- M$U$maximum; kU <- nfree(M$U)
tst <- do.call(rbind, lapply(jobs[-1], function(j) { x <- 2 * (LLU - M[[j]]$maximum); d <- kU - nfree(M[[j]]); data.frame(restriction = lab[[j]], LR = round(x, 1), df = d, p = signif(pchisq(x, d, lower.tail = FALSE), 2)) }))
ps <- tst$p[tst$restriction != lab[["all"]]]
tst$p_Holm <- NA; tst$p_Holm[tst$restriction != lab[["all"]]] <- signif(p.adjust(ps, "holm"), 2)
cat("\nLikelihood-ratio tests against the unrestricted pooled model (= solar and wind separate, same draws)\n"); print(tst, row.names = FALSE)
cat("Holm adjustment covers the five single-block tests; 'all' is the joint test.\n")

# ---- differences in WTA with delta-method s.e. (robust covariance of the unrestricted pooled model)
e <- M$U$estimate; V <- M$U$robvarcov; nv <- rownames(V)
dm <- function(f) { J <- sapply(nv, function(n) { h <- 1e-5 * max(1, abs(e[[n]])); tp <- e; tp[n] <- e[[n]] + h; tn <- e; tn[n] <- e[[n]] - h; (f(tp) - f(tn)) / (2 * h) }); c(est = f(e), se = sqrt(drop(J %*% V %*% J))) }
at <- list(`status quo (reference respondent)` = c("m_sq", +1), `no protection` = c("m_nop", -1), `< 1 km` = c("m_d1", -1), `visible` = c("m_view", -1),
           `production 38k` = c("w_el38", -1), `production 57k` = c("w_el57", -1), `plants moderate` = c("w_pl3", -1), `plants many` = c("w_pl6", -1),
           `wilderness only` = c("w_wil", -1), `1-3 km` = c("w_d2", -1), `3-6 km` = c("w_d3", -1))
tb <- do.call(rbind, lapply(names(at), function(n) { p <- at[[n]][1]; s <- as.numeric(at[[n]][2]); xp <- paste0("x_", p)
  sol <- dm(function(th) 100 * s * th[[p]]); dif <- dm(function(th) 100 * s * th[[xp]])
  data.frame(attribute = n, solar = round(sol["est"]), wind = round(sol["est"] + dif["est"]), difference = round(dif["est"]), se = round(dif["se"]), z = round(dif["est"] / dif["se"], 2), row.names = NULL) }))
cat("\nWTA in NOK per household-year, solar vs wind from the unrestricted pooled model (difference = wind - solar; delta-method s.e.)\n"); print(tb, row.names = FALSE)
r <- dm(function(th) (th[["m_sq"]] + th[["x_m_sq"]]) / th[["m_sq"]])
r35 <- dm(function(th) (th[["m_sq"]] + th[["x_m_sq"]]) - 3.5 * th[["m_sq"]])
cat(sprintf("\nStatus quo constant, wind / solar (reference respondent): %.2f (delta s.e. %.2f); test of the article's '3.5 times': wind - 3.5 x solar = %.0f NOK (s.e. %.0f, z %.2f)\n",
            r["est"], r["se"], 100 * r35["est"], 100 * r35["se"], r35["est"] / r35["se"]))
# ---- scale
sc <- dm(function(th) th[["x_m_cost"]]); sd5 <- dm(function(th) th[["l55"]] + th[["x_l55"]]); sd5s <- dm(function(th) th[["l55"]])
cat(sprintf("\nScale (log of the cost coefficient): solar mean %.2f, wind mean %.2f; wind - solar = %.2f (delta s.e. %.2f, z %.2f) i.e. a factor %.2f on the cost coefficient / scale\n",
            e[["m_cost"]], e[["m_cost"]] + e[["x_m_cost"]], sc["est"], sc["se"], sc["est"] / sc["se"], exp(sc["est"])))
cat(sprintf("Spread of the log-scale (Cholesky l55): solar %.2f, wind %.2f\n", abs(sd5s["est"]), abs(sd5["est"])))
cat("Cost-scale correlations (solar | wind) are in the output of step 04; here: the 'scale' LR test above covers mean, spread and correlations jointly.\n")
