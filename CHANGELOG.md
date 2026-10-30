# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

- Accelerate parallel tempering with exact integer contact-energy caches and
  tabulated Boltzmann factors for the bundled MJ interactions; retain the general
  floating-point kernel and an explicit cache-selection option.
- Add checkpointed continuation and convergence diagnostics for the full-model
  Tubiana figure, including independent-ladder R-hat, effective sample sizes,
  conditional-frequency uncertainties, and downloadable trajectories.
- Replace the full-structure Tubiana figure with the PT sample that passes all
  monitored convergence criteria: 316,000 draws, with 172,717 at `P_nat > 0.99`.
  Document how to redraw it from the saved alignment.
- Validate parallel tempering against an exactly enumerable model and direct
  energies; preserve continuous trajectories across checkpoint boundaries.

## 2.1.0

- Add full-model parallel tempering with reproducible RNGs, cached mutation
  energies, temperature traces, and separate burn-in/production exchange and
  round-trip diagnostics. Include an equal-budget Metropolis comparison script.

## 1.0.0 - 2023-08-03

- No changes.

## 0.1.0

- Register at https://github.com/cossio/CossioJuliaRegistry.
