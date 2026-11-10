# =====================================================================
# 14  Figure: acceptance curves (share of households accepting a scenario) against annual bill reduction
# =====================================================================
# Same model and delta-method standard errors as step 09 part C (LC2-D, population-weighted).  The code of step 09 is sourced silently so
# that the function 'accept' and the delta-method function 'dm' are identical to the ones used for the tables.
invisible(capture.output(source(file.path(V2, "09_policy.R"), echo = FALSE)))
cvg <- seq(0, 50000, by = 2500)
sc <- list(P6 = "only production 38k, no visibility (cheapest scenario)", P1 = "38k, < 1 km from home", P4 = "38k, many plants, no protection")
cols <- c(P6 = "#1f6f8b", P1 = "#d9822b", P4 = "#7a7a7a")
png(file.path(V2, "output", "fig_acceptance.png"), width = 2000, height = 900, res = 220)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.6, 0.8), mgp = c(2.6, 0.7, 0), las = 1)
for (s in c("solar", "wind")) {
  m <- lc[[s]]$lD
  plot(NA, xlim = c(0, 50), ylim = c(0, 60), xlab = "Annual bill reduction (1,000 NOK per household)", ylab = "Share accepting the project (%)", main = toupper(s), cex.main = 1)
  abline(v = 9, lty = 3, col = "grey60"); text(9.6, 58, "highest amount on the cards (9,000 NOK)", adj = 0, cex = 0.6, col = "grey40")
  for (p in names(sc)) {
    res <- sapply(cvg, function(cc) dm(m, function(th) accept(th, s, scen[[p]], FALSE, cc, TRUE)))
    est <- 100 * res["est", ]; se <- 100 * res["se", ]
    polygon(c(cvg / 1000, rev(cvg / 1000)), c(pmax(est - 1.96 * se, 0), rev(est + 1.96 * se)), col = adjustcolor(cols[[p]], 0.18), border = NA)
    lines(cvg / 1000, est, col = cols[[p]], lwd = 2.2)
  }
  if (s == "solar") legend("bottomright", legend = unname(sc), col = cols, lwd = 2.2, bty = "n", cex = 0.6)
}
dev.off()
cat("figure written to", file.path(V2, "output", "fig_acceptance.png"), "\n")
