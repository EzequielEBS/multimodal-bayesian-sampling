# Multimodal Bayesian Sampling

This repository contains code and experiments related to sampling from multimodal posterior distributions in Bayesian mixture models. We compare standard Gibbs sampling with advanced MCMC methods including **parallel tempering**, **simulated tempering**, and **tempered transitions**. The goal is to evaluate their effectiveness in handling the label switching problem and exploring multimodal posteriors.

## Contents

- `code/`: R scripts.
- `data/`: Simulated datasets for a 4-component Gaussian mixture.
- `figures/`: Generated plots including trace plots and posterior distributions.
- `main.pdf`: Report

## Methods Implemented

- **Gibbs Sampling** with latent variables and semi-conjugate priors.
- **Parallel Tempering** with multiple chains at different temperatures and swap proposals.
- **Tempered Transitions** following Neal (1996) with bidirectional transition kernels.
- Post-processing for **label switching**, including imposing identifiability constraints.

## Simulated Data

We simulate data from a 4-component Normal mixture and assess the performance of each method based on:
- Trace plots
- Posterior histograms
- Acceptance rates
- Effective sample sizes (ESS)
