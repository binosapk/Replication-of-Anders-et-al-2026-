# Re-analysis of Dugstad, Hoel and Lindhjem (2026): solar versus wind power

Re-analysis of the discrete choice experiment in Dugstad, A., Hoel, J. and Lindhjem, H. (2026), *Kraft med mindre konflikt? Preferanser for sol- vs. vindkraft*, Samfunnsøkonomen 138(2), with the same data and my own R code.

**Contributor:** Binod

**Data:** the survey and design files of the original study (<https://github.com/dugstad/Solar-vs-Wind-Experiment>) are copied into this repository, so the code runs as it is.

## What is in this repository

| Path | Content |
|---|---|
| `sol_complete.csv`, `vind_complete.csv`, `CE_design_final.csv` | Survey responses (solar, wind) and choice-experiment design, from the original study |
| `reanalysis_v2/Solar_vs_Wind_Reanalysis.R` | The R code: one script with one section per analysis step |
| `reanalysis_v2/ssb_pop.csv` | Population margins (Statistics Norway, Innlandet, age 18+) used for the raking weights |
| `reanalysis_v2/output/` | Text output of each step, and the acceptance-curve figure `fig_acceptance.png` (the only figure the scripts produce) |

Not included: the saved model objects (`state_*.RData`, Apollo output), which the scripts create when they run, and the first-round pooled mixed-logit outputs. The first-round pooled runs stopped early and are invalid (see "Convergence of the pooled mixed logit").

## How to run

Software: R 4.6.1 with the packages `apollo`, `dplyr`, `tidyr` and `car`. The script finds its own folder; there are no paths to edit.

From the repository folder:

1. `Rscript reanalysis_v2/Solar_vs_Wind_Reanalysis.R 00 all` builds the data, covariates and raking weights.
2. `Rscript reanalysis_v2/Solar_vs_Wind_Reanalysis.R 02 all` (or any other step and sample) runs one step. The output goes to the console; redirect it to a file, for example `Rscript reanalysis_v2/Solar_vs_Wind_Reanalysis.R 13 all > reanalysis_v2/output/13.txt`. The header of the script lists all steps and samples.
3. Each step saves its results in `state_<step>_<sample>.RData`, which later steps load (for example 04 loads the results of 02, 05 and 06b load those of 06, 13 loads those of 11 and 15, and 09 loads results of earlier model steps).
4. Checks after the original step scripts were merged into this one file: steps 00, 01, 02, 09, 10, 12 and 13 give exactly the same output as before, and step 14 gives a byte-identical figure. Every slow step (03 to 08, 11, 15 and 17) was run once in quick test mode (`SW_QUICK=1`: 10 draws and 2 iterations, writing only `quick_*` files), which shows that the code runs but gives no results. The slow estimations were not repeated.
5. The random-coefficient models are slow on one core: about 8 hours for the replication (step 03) and about 39 to 69 hours per job for the pooled mixed logit (step 11).

## Steps in the script

Each step is one section of `Solar_vs_Wind_Reanalysis.R`; the published Table A.3 is defined once at the top.

| Step | Content | Sample or job |
|---|---|---|
| 00 | Data, covariates, raking weights | any |
| 01 | Descriptive statistics | all |
| 02 | MNL ladder: shifters, income, attitudes | all |
| 03 | Benchmark mixed-logit replication, 500 MLHS draws (follows the specification and Apollo structure of the published script) | solar, wind |
| 04 | Full-covariance mixed logit with random cost, staged estimation, 250 MLHS draws | solar, wind |
| 05 | Pooled two-class latent class | both |
| 06, 06b | Latent classes LC2 and LC3 (06b: trader cost coefficient fixed) | both |
| 07 | Hybrid choice model with five indicators | solar, wind |
| 08 | Bayesian full-covariance model (did not converge; exploratory) | solar, wind |
| 09 | Cost per kWh, acceptance curves; check of Figure 6 | all |
| 10 | Respondent exclusion rules | all |
| 11 | Pooled mixed logit with all restrictions | job: U2, scale2, cov2, attr2, sq2, all2 |
| 12 | Pooled MNL scale test | all |
| 13 | Likelihood-ratio tests and Holm correction for the pooled mixed logit | all |
| 14 | Acceptance curves (figure) | all |
| 15 | Pooled mixed logit with a scaled optimiser (the unrestricted reference model U3) | job: U3, scale3, cov3, attr3, sq3, all3 |
| 16 | Collinearity check of the covariates | all |
| 17 | Draw-seed check of the pooled mixed logit (seed 2027). Started but not completed: there is no result, so the tests above rest on one set of draws | job: U4, all4 |

## Main results

**Data.** 775 respondents (solar) and 788 (wind), eight choice tasks each: 6,200 and 6,304 observations. Utility is written in willingness-to-accept (WTA) space; units are NOK per household and year. Respondents who chose "no development" in all eight tasks (35.9% solar, 52.0% wind) are kept in the main sample. Raking weights (age by sex, education) give an effective sample of about 398 (solar) and 365 (wind).

**Published figures and replication.** Figure 6 of the article is reproduced exactly from Table A.3 for all twelve scenario-and-visibility combinations (maximum difference 0 øre/kWh).

| | Solar: published | Solar: replication | Wind: published | Wind: replication |
|---|---:|---:|---:|---:|
| Log-likelihood | −2536.4 | −2531.3 | −2065.2 | −2054.7 |
| BIC | 5,335 | 5,324.6 | 4,392 | 4,371.9 |
| Parameters | 30 | 30 | 30 | 30 |

Points raised about the article: the worked example for wind in Table A.6 uses the 57,000-household coefficient (17.16 øre/kWh) while Figure 6 uses the 38,000-household coefficient (17.62 øre/kWh); the solar status quo value is NOK 3,822 in the text and NOK 3,872 in Table A.3; the published standard errors for wind are much smaller than robust ones from the replication (status quo mean: 0.703 against 5.70, in units of 100 NOK). The replication script also contains a block labelled "Figure 8" whose delta-method expression is divided by an unexplained constant (13.88); the published article has no Figure 8 and no such constant, and that block was not used.

**Latent classes (LC3).** Shares: cost-responsive traders 25.4% (solar) and 15.3% (wind); cost-ignorers 36.9% and 27.9%; status-quo-only 37.8% and 56.8%. The trader cost coefficient is fixed at its two-class value because the unrestricted wind model collapsed toward a second cost-ignorer class. Trader WTAs are imprecise.

**Policy.** For scenario P6 without visibility, cost per kWh is 3.3 (solar) and 9.2 (wind) øre/kWh in the article. Acceptance curves show that a large status-quo-only class remains unlikely to accept even at high compensation (at NOK 50,000 per year: 50.0% solar, 39.2% wind, scenario P6).

### Tests of differences between solar and wind

Pooled MNL (step 12) and pooled two-class latent class (step 05):

| Test | LR | df | p |
|---|---:|---:|---:|
| MNL: all eleven WTAs equal (scale free) | 43.4 | 11 | 9.2e-06 |
| MNL: status quo term equal (scale free) | 29.7 | 1 | 5.1e-08 |
| MNL: other ten WTAs equal | 13.7 | 10 | 0.19 |
| MNL: scale equal | 0.7 | 1 | 0.41 |
| Latent class: technology-specific constants | 45.9 | 2 | 1.1e-10 |
| Latent class: ten trader attribute WTAs equal | 21.7 | 10 | 0.016 |

Pooled mixed logit (steps 11, 13, 15; 12,504 observations, 250 MLHS draws). The unrestricted model has 68 parameters and is exactly the two separate models on the same draws; each restriction sets one block of wind-minus-solar differences to zero.

| Restriction | Log-lik. | Parameters | LR | df | p | p (Holm) |
|---|---:|---:|---:|---:|---:|---:|
| Unrestricted (U3) | −4552.2 | 68 | | | | |
| Cost scale equal | −4552.4 | 62 | 0.4 | 6 | 1.0 | 1.0 |
| Taste covariance equal | −4565.6 | 58 | 26.8 | 10 | 0.0028 | 0.0084 |
| Ten attribute WTA means equal | −4572.7 | 58 | 41.1 | 10 | 1.1e-05 | 4.4e-05 |
| Status quo constant equal | −4572.7 | 67 | 41.0 | 1 | 1.5e-10 | 7.5e-10 |
| Seven status quo shifters equal | −4557.1 | 61 | 9.9 | 7 | 0.19 | 0.38 |
| Everything equal (joint test) | −4621.8 | 34 | 139.3 | 34 | 1.1e-14 | |

The technologies differ in the status quo constant, the means of the attribute valuations and the taste covariance, but not in the cost scale or the status quo shifters. In the pooled MNL only the status quo constant differs.

BIC is lowest for the model with everything equal (9,564.4 against 9,745.8 for the unrestricted model) because the penalty for 34 extra parameters exceeds the likelihood gain; the likelihood-ratio test nevertheless rejects that restriction.

### Convergence of the pooled mixed logit

- In the first round the model with equal shifters (−4557.14) fitted better than the unrestricted model that contains it (−4568.42), so the first-round unrestricted model had stopped early. Tests against it were discarded.
- All models were restarted from the best solution (jobs ending in 2). The unscaled unrestricted restart (−4552.42) was still below the restricted scale model (−4552.39).
- The unrestricted model was therefore re-estimated with parameter scaling in the optimiser (step 15; same model, draws and likelihood). It reached −4552.17, above every restricted model, and is the reference in step 13.

## Limits

- The covariance matrix of the pooled mixed logit is not positive definite (smallest eigenvalue −1.1, five undefined standard errors) and the robust standard errors of the main coefficients are implausibly small. Delta-method z-values for the solar-wind differences are therefore not reported; the likelihood-ratio tests do not use them.
- The tests rest on one set of 250 MLHS draws. The scaled and unscaled estimates of the unrestricted model differ by 0.25 log-likelihood points. No re-estimation with another draw seed was done.
- No bootstrap of standard errors was done. Standard errors of random-coefficient models, like the published ones, are probably too small.
- The point estimates of the spread of the log-scale are 0.41 (solar) and 2.11 (wind), but the joint scale test does not reject equality (LR 0.4): this spread is weakly identified and is not interpreted.
- Random-coefficient models use 250 to 500 draws (the article 2,000). The Bayesian model did not converge. Attitude variables are measured after the choice tasks and may partly rationalise the choices; BIC comparisons apply to choice models only (the hybrid likelihood includes the indicators).
