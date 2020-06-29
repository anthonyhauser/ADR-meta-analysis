# ADR-meta-analysis
Acquired HIV drug resistance mutations on first-line antiretroviral therapy in Southern Africa: Bayesian evidence synthesis

We searched several bibliographic databases, including Embase and Medline, from inception to May 2019 to identify studies reporting NRTI DRMs observed among adult HIV-positive patients experiencing virological failure on first-line NNRTI-based regimens in countries of the Southern African IeDEA region. After screening titles and abstracts, two independent reviewers assessed full manuscripts of potentially eligible studies and extracted data. We developed a hierarchical logistic meta-regression model to synthesize the effect of different NRTI drug combinations on DRM emergence across studies, accounting for ART duration and study-specific effects. Analyses were performed in a Bayesian framework using the rstan package in R.

The model assessed both the baseline (before treatment) prevalence of eight NRTI DRMs and the prevalence after 3 years, by first-line regimen.
