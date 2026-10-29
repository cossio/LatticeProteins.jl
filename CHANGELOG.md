# Changelog

All notable changes to this project will be documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

- Parallelize local replica updates with temperature-specific RNG streams, producing
  identical results across thread counts and with `threaded=false`. Seeded trajectories
  differ from v2.1.0. Share the stable folding-loss calculation across samplers and
  probability evaluation, and add an exact-distribution sampling check.

## 2.1.0

- Add full-model parallel tempering with reproducible RNGs, cached mutation
  energies, temperature traces, and separate burn-in/production exchange and
  round-trip diagnostics. Include an equal-budget Metropolis comparison script.

## 1.0.0 - 2023-08-03

- No changes.

## 0.1.0

- Register at https://github.com/cossio/CossioJuliaRegistry.
