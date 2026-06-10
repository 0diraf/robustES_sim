# Robust Effect Size Estimation Under Heavy-Tailed Distributions

Code accompanying the paper 'Robust Effect Size Estimation in Heavy-Tailed Data: Evaluating Alternatives to Cohen's d', containing simulation code for comparing the performance of Cohen's d, AKP robust d, and Blaine's d under heavy-tailed treatment data.

## Overview

Standard effect size estimators like Cohen's d rely on the mean and standard deviation, which are sensitive to heavy tails. This simulation evaluates whether median-based alternatives (Blaine's d) remain stable where Cohen's d does not.

The treatment group is generated as a contaminated mixture: `(1 - f) * N(μ, 1) + f * D_tail`, where `D_tail` is either a Pareto or Lognormal distribution. The contamination fraction `f` is varied from 1% to 10%. Estimator output is compared against the bulk location μ (the center of the non-contaminated component).

## Estimators

- **Cohen's d**: standardized mean difference using pooled SD
- **AKP robust d**: 20% trimmed means with Winsorized SD, rescaled by 0.642
- **Blaine's d (Δ_MAD)**: absolute median difference standardized by pooled MAD

## Simulation design

- **Sample size**: n = 10,000 per group
- **Effect sizes**: small (μ = 0.2), medium (μ = 0.5), large (μ = 0.8)
- **Pareto contamination**: α ∈ {0.5, 1.5, 2.5} (infinite mean/variance → finite mean/variance)
- **Lognormal contamination**: σ_log ∈ {1, 1.8, 3.0} (moderate → extreme tail heaviness)
- **Iterations**: 10,000 per condition
- **Contamination**: right-sided, treatment group only


## Reference

Ricca, B. & Blaine, B. (2022). Notes on a nonparametric estimate of effect size. *Journal of Experimental Education*.

