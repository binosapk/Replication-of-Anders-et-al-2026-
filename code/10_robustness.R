# =====================================================================
# 10  Robustness: who is removed, and does belief in the survey's consequences matter?
# =====================================================================
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
  apollo_attach(apollo_beta, apollo_inputs); on.exit(apollo_detach(apollo_beta, apollo_inputs))
  V <- list()
  V[["dev"]] <- b_cost * (w_el38 * (el_alt1 == 38) + w_el57 * (el_alt1 == 57) + w_pl3 * (plants_alt1 == 3) + w_pl6 * (plants_alt1 == 6) +
                          w_nop * (wild_alt1 == 0) + w_wil * (wild_alt1 == 1) + w_d1 * (dist_alt1 == 1) + w_d2 * (dist_alt1 == 2) +
                          w_d3 * (dist_alt1 == 3) + w_view * (view_alt1 == 2) + cost_alt1 / 100)
  V[["sq"]] <- b_cost * w_sq
  P <- list(model = apollo_mnl(list(alternatives = c(dev = 1, sq = 2), avail = list(dev = 1, sq = 1), choiceVar = choice, utilities = V), functionality))
  P <- apollo_panelProd(P, apollo_inputs, functionality); P <- apollo_prepareProb(P, apollo_inputs, functionality); return(P)
}
st0 <- setNames(c(rep(0, 11), 0.005), wn)
key <- c("w_sq", "w_nop", "w_d1", "w_view")
rob <- list()
for (s in c("solar", "wind")) {
  d <- dat[[s]]; cut10 <- quantile(d$qtime[!duplicated(d$id)], 0.10)
  vv <- list(`all respondents` = d,
             `Anders's protest rule removed` = d[!grepl(pat_anders, d$Q27), ],
             `protest rule without protected-value reasons` = d[!grepl(pat_noval, d$Q27), ],
             `no always-status-quo respondents` = d[d$always_sq == 0, ],
             `only believers in consequences (Q41)` = d[d$conseq == 1, ],
             `without 10% fastest` = d[d$qtime > cut10, ],
             `without 'difficult' (6-7)` = d[is.na(d$difficult) | d$difficult < 6, ])
  out <- lapply(names(vv), function(v) {
    m <- fit_model(paste0("rb_", s, "_", abs(sum(utf8ToInt(v)))), vv[[v]], prob_b, st0)
    e <- m$estimate; se <- rse(m); sgn <- c(100, -100, -100, -100); names(sgn) <- key
    data.frame(variant = v, respondents = length(unique(vv[[v]]$id)),
               status_quo = sprintf("%6.0f (%4.0f)", sgn[["w_sq"]] * e[["w_sq"]], 100 * se[["w_sq"]]),
               no_protection = sprintf("%6.0f (%4.0f)", sgn[["w_nop"]] * e[["w_nop"]], 100 * se[["w_nop"]]),
               lt_1km = sprintf("%6.0f (%4.0f)", sgn[["w_d1"]] * e[["w_d1"]], 100 * se[["w_d1"]]),
               visible = sprintf("%6.0f (%4.0f)", sgn[["w_view"]] * e[["w_view"]], 100 * se[["w_view"]]))
  })
  rob[[s]] <- do.call(rbind, out)
  cat("\n=====", toupper(s), " WTA in NOK per household-year (robust s.e.) =====\n"); print(rob[[s]], row.names = FALSE)
}
cat("\nShare of respondents removed by the two protest rules (Anders's / without protected-value reasons):\n")
print(round(sapply(c("solar", "wind"), function(s) { q <- dat[[s]][!duplicated(dat[[s]]$id), "Q27"]; c(grepl(pat_anders, q) |> mean(), grepl(pat_noval, q) |> mean()) }), 3))
cat("Share of the always-status-quo respondents whose reason is the protected-value statement 'Supplying Europe ... at the expense of Norwegian nature':\n")
print(round(sapply(c("solar", "wind"), function(s) { r <- dat[[s]][!duplicated(dat[[s]]$id), ]; r <- r[r$always_sq == 1, ]; mean(grepl("bekostning av norsk", r$Q27)) }), 3))
