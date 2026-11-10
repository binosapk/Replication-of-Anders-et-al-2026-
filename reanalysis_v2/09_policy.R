# =====================================================================
# 09  Cost per kWh with the delta method (as in the article), across models, and a population-weighted acceptance curve
# =====================================================================
# Anders: cost/kWh = 100 * ( w_sq - sum of the w's of the scenario ) NOK * 4,100 households / 600,000,000 kWh, in ore; delta method on the robust covariance.
# (A) reproduce his Figure 6 from the published Table A.3                              (checks the formula and the units)
# (B) the same formula, delta method, for every model of this re-analysis
# (C) acceptance curves, post-stratified to the population (Leeds [A] 1B slides 54-56: weighting AFTER estimation, by sample enumeration)
ld <- function(f) { e <- new.env(); load(file.path(V2, f), envir = e); e }
source(file.path(V2, "a3_published.R"))
e02 <- ld("state_02_all.RData"); lad <- e02$lad
e06 <- ld("state_06_both.RData"); lc <- e06$lc
rep_ <- lapply(c(solar = "solar", wind = "wind"), function(s) ld(sprintf("state_03_%s.RData", s))$mA)
mvn  <- lapply(c(solar = "solar", wind = "wind"), function(s) ld(sprintf("state_04_%s.RData", s))$mM)

H <- 4100; KWH <- 600e6
scen <- list(P1 = c("el38", "d1"), P2 = c("el38", "pl6", "d1"), P3 = c("el38", "nop"), P4 = c("el38", "pl6", "nop"), P5 = c("el38", "pl3", "wil", "d3"), P6 = "el38")
fig6 <- list(wind  = rbind(vis = c(17.6, 18.2, 14.7, 15.3, 14.2, 11.1), novis = c(15.7, 16.3, 12.8, 13.3, 12.3, 9.2)),
             solar = rbind(vis = c(9.5, 10.7, 11.6, 12.8, 8.6, 6.5),   novis = c(6.3, 7.5, 8.4, 9.6, 5.4, 3.3)))
# generic cost: z = list(sq = status quo constant, w = named vector of attribute w's, view)
ore <- function(z, lev, visible) 100 * (z$sq - sum(z$w[lev]) - if (visible) z$view else 0) * H / KWH * 100

# ------------------------------------------------------------------ (A) Figure 6 from Table A.3
cat("(A) Figure 6 recomputed from the published Table A.3 (ore/kWh) vs the figure in the article\n")
chk <- do.call(rbind, lapply(c("wind", "solar"), function(s) { a <- A3[[s]]; z <- list(sq = a$mean[1], w = setNames(a$mean[2:11], a$par[2:11]), view = a$mean[11])
  z$w <- setNames(a$mean[2:10], c("el38", "el57", "pl3", "pl6", "nop", "wil", "d1", "d2", "d3")); z$view <- a$mean[11]
  do.call(rbind, lapply(c(TRUE, FALSE), function(v) data.frame(sample = s, visible = v, scenario = names(scen),
         recomputed = round(sapply(scen, function(l) ore(z, l, v)), 1), article = fig6[[s]][ifelse(v, "vis", "novis"), ])))}))
chk$diff <- round(chk$recomputed - chk$article, 1)
print(reshape(chk[, c("sample", "visible", "scenario", "recomputed", "article")], idvar = c("scenario", "visible"), timevar = "sample", direction = "wide"), row.names = FALSE)
cat("max absolute difference:", max(abs(chk$diff)), "ore/kWh\n")
cat("\nWorked example, Table A.6 of the article (alternative 1 with visibility, i.e. < 1 km, all protected, few plants):\n")
for (s in c("solar", "wind")) { a <- A3[[s]]; el <- a$mean[2]; el57 <- a$mean[3]; d1 <- a$mean[8]; vw <- a$mean[11]     # index 8 = < 1 km, 11 = visible
  s38 <- 100 * (a$mean[1] - el - d1 - vw); s57 <- 100 * (a$mean[1] - el57 - d1 - vw)
  cat(sprintf("  %-5s with the el38 coefficient (%.0f NOK): sum %.0f NOK -> %.2f ore/kWh;  with the el57 coefficient (%.0f NOK): sum %.0f NOK -> %.2f ore/kWh\n", s,
      -100 * el, s38, s38 * H / KWH * 100, -100 * el57, s57, s57 * H / KWH * 100)) }

# ------------------------------------------------------------------ delta method (numerical Jacobian; identical to car::deltaMethod for smooth functions)
dm <- function(m, f) {
  th <- m$estimate; V <- m$robvarcov; nmv <- rownames(V)
  J <- sapply(nmv, function(n) { h <- 1e-5 * max(1, abs(th[[n]])); tp <- th; tp[n] <- th[[n]] + h; tn <- th; tn[n] <- th[[n]] - h; (f(tp) - f(tn)) / (2 * h) })
  c(est = f(th), se = sqrt(drop(J %*% V %*% J)))
}
wbar <- function(s, v) { r <- dat[[s]][!duplicated(dat[[s]]$id), ]; sum(r[[v]] * r$wt) / sum(r$wt) }
zb <- function(s) sapply(c("age_c", "female", "univ", "inc_c", "inc_missing", "exposed", "org"), function(v) wbar(s, v))
mk_z <- list(
  MNL    = function(th, s) list(sq = th[["w_sq"]], w = c(el38 = th[["w_el38"]], el57 = th[["w_el57"]], pl3 = th[["w_pl3"]], pl6 = th[["w_pl6"]], nop = th[["w_nop"]], wil = th[["w_wil"]], d1 = th[["w_d1"]], d2 = th[["w_d2"]], d3 = th[["w_d3"]]), view = th[["w_view"]]),
  LCtr   = function(th, s) list(sq = th[["w_sq"]], w = c(el38 = th[["w_el38"]], el57 = th[["w_el57"]], pl3 = th[["w_pl3"]], pl6 = th[["w_pl6"]], nop = th[["w_nop"]], wil = th[["w_wil"]], d1 = th[["w_d1"]], d2 = th[["w_d2"]], d3 = th[["w_d3"]]), view = th[["w_view"]]),
  REP    = function(th, s) list(sq = th[["m_sq"]], w = c(el38 = th[["m_el38"]], el57 = th[["m_el57"]], pl3 = th[["m_pl3"]], pl6 = th[["m_pl6"]], nop = th[["m_nop"]], wil = th[["m_wil"]], d1 = th[["m_d1"]], d2 = th[["m_d2"]], d3 = th[["m_d3"]]), view = th[["m_view"]]),
  MVN    = function(th, s) { z <- zb(s); g <- c(th[["g_age"]], th[["g_fem"]], th[["g_univ"]], th[["g_inc"]], th[["g_incm"]], th[["g_exp"]], th[["g_org"]])
             list(sq = th[["m_sq"]] + sum(g * z), w = c(el38 = th[["w_el38"]], el57 = th[["w_el57"]], pl3 = th[["w_pl3"]], pl6 = th[["w_pl6"]], nop = th[["m_nop"]], wil = th[["w_wil"]], d1 = th[["m_d1"]], d2 = th[["w_d2"]], d3 = th[["w_d3"]]), view = th[["m_view"]]) })
models <- function(s) list(MNL = list(lad[[s]]$m0, "MNL"), `LC traders (LC2-DA)` = list(lc[[s]]$lDA, "LCtr"), `MMNL replication of Anders` = list(rep_[[s]], "REP"), `MMNL full covariance` = list(mvn[[s]], "MVN"))

cat("\n(B) Cost per kWh (ore), delta-method s.e. in brackets, for every model; scenario x visibility\n")
tabB <- list()
for (s in c("solar", "wind")) { rows <- list()
  for (mn in names(models(s))) { m <- models(s)[[mn]][[1]]; ty <- models(s)[[mn]][[2]]
    for (v in c(FALSE, TRUE)) rows[[length(rows) + 1]] <- c(model = mn, visible = ifelse(v, "yes", "no"),
      sapply(names(scen), function(p) { r <- dm(m, function(th) ore(mk_z[[ty]](th, s), scen[[p]], v)); sprintf("%5.1f (%.1f)", r["est"], r["se"]) })) }
  tabB[[s]] <- as.data.frame(do.call(rbind, rows)); cat("\n=====", toupper(s), "=====\n"); print(tabB[[s]], row.names = FALSE) }

cat("\nRatio of the status-quo term, wind / solar ('3.5 times higher' in the article), delta method (independent samples):\n")
for (mn in c("MNL", "MMNL replication of Anders")) { ty <- models("solar")[[mn]][[2]]
  a <- dm(models("solar")[[mn]][[1]], function(th) mk_z[[ty]](th, "solar")$sq); b <- dm(models("wind")[[mn]][[1]], function(th) mk_z[[ty]](th, "wind")$sq)
  r <- b["est"] / a["est"]; se <- r * sqrt((b["se"] / b["est"])^2 + (a["se"] / a["est"])^2)
  cat(sprintf("  %-28s solar %.0f NOK (se %.0f), wind %.0f NOK (se %.0f); ratio %.2f (delta se %.2f)\n", mn, 100 * a["est"], 100 * a["se"], 100 * b["est"], 100 * b["se"], r, se)) }

# ------------------------------------------------------------------ (C) acceptance curves from the exogenous-only latent class model, population-weighted
cat("\n(C) Share of households accepting a scenario at a given annual bill reduction, LC2-D (class shares by exogenous characteristics only)\n")
accept <- function(th, s, lev, visible, cvals, weighted) {
  r <- dat[[s]][!duplicated(dat[[s]]$id), ]
  v2 <- th[["d_cls"]] + th[["gd_age"]] * r$age_c + th[["gd_fem"]] * r$female + th[["gd_univ"]] * r$univ + th[["gd_inc"]] * r$inc_c + th[["gd_incm"]] * r$inc_missing + th[["gd_exp"]] * r$exposed + th[["gd_org"]] * r$org
  pi2 <- plogis(v2); wts <- if (weighted) r$wt else rep(1, nrow(r))
  z <- mk_z$LCtr(th, s); Sw <- sum(z$w[lev]) + if (visible) z$view else 0
  sapply(cvals, function(cc) sum(wts * ((1 - pi2) * plogis(th[["b_cost"]] * (Sw + cc / 100 - z$sq)) + pi2 * (1 - plogis(th[["a_sq2"]])))) / sum(wts))
}
cv <- c(0, 2000, 5000, 9000, 15000, 25000, 50000)
for (s in c("solar", "wind")) { m <- lc[[s]]$lD
  for (p in c("P6", "P1", "P4")) {
    out <- do.call(cbind, lapply(c(FALSE, TRUE), function(w) sapply(seq_along(cv), function(i) { r <- dm(m, function(th) accept(th, s, scen[[p]], FALSE, cv[i], w)); c(r[["est"]], r[["se"]]) })))
    cat(sprintf("\n[%s, %s, not visible]  compensation (NOK) | accepted share, population-weighted (delta s.e.) | unweighted\n", toupper(s), p))
    print(data.frame(NOK = cv, weighted = sprintf("%4.1f%% (%.1f)", 100 * out[1, seq_along(cv)], 100 * out[2, seq_along(cv)]),
                     unweighted = sprintf("%4.1f%% (%.1f)", 100 * out[1, length(cv) + seq_along(cv)], 100 * out[2, length(cv) + seq_along(cv)])), row.names = FALSE)
  } }
cat("\nShare of households in the status-quo-only class (LC2-D), weighted vs unweighted:\n")
print(round(sapply(c(solar = "solar", wind = "wind"), function(s) { th <- lc[[s]]$lD$estimate; r <- dat[[s]][!duplicated(dat[[s]]$id), ]
  p <- plogis(th[["d_cls"]] + th[["gd_age"]] * r$age_c + th[["gd_fem"]] * r$female + th[["gd_univ"]] * r$univ + th[["gd_inc"]] * r$inc_c + th[["gd_incm"]] * r$inc_missing + th[["gd_exp"]] * r$exposed + th[["gd_org"]] * r$org)
  c(unweighted = mean(p), weighted = sum(p * r$wt) / sum(r$wt)) }), 3))
