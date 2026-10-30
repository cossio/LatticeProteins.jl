"""Check independent PT ladders before interpreting a filtered alignment.

Usage: python diagnose_pt.py OUT [--threshold .99] [--discard 0]
Requires Python >=3.11, numpy, scipy, pandas, matplotlib and arviz.
The reported gates are evidence for the monitored observables, not a proof that
every sequence family has been explored. Rejected draws remain in time-series
diagnostics; filtering must not conceal drift or compress the sampling clock.
"""
import argparse
import hashlib
import json
from pathlib import Path
import tomllib

import arviz as az
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import eigh

AA = "CMFILVWYAGTSNQDEHRKP"
CODE = {a: i for i, a in enumerate(AA)}


def diagnostic(values):
    values = np.asarray(values, dtype=float)
    if np.ptp(values) == 0:
        return {"rhat": None, "bulk_ess": None, "mcse": None}
    return {"rhat": float(az.rhat(values, method="rank")),
            "bulk_ess": float(az.ess(values, method="bulk")),
            "mcse": float(az.mcse(values, method="mean"))}


def encode(sequences):
    return np.array([[CODE[a] for a in s] for s in sequences], dtype=np.int8)


def frequencies(sequences):
    return np.concatenate([np.bincount(sequences[:, i], minlength=20)
                           for i in range(sequences.shape[1])]) / len(sequences)


def sequence_covariance(sequences):
    """The same one-hot covariance, without allocating a 540 x sample array."""
    sequences = np.asarray(sequences, dtype=np.int16)
    sites = sequences.shape[1]
    center = frequencies(sequences)
    second = np.empty((sites, 20, sites, 20))
    for i in range(sites):
        for j in range(i, sites):
            joint = np.bincount(sequences[:, i] * 20 + sequences[:, j],
                                minlength=400).reshape(20, 20) / len(sequences)
            second[i, :, j, :] = joint
            second[j, :, i, :] = joint.T
    return center, second.reshape(sites * 20, sites * 20) - np.outer(center, center)


def project_sequences(sequences, center, vectors):
    sites = sequences.shape[-1]
    weights = vectors.reshape(sites, 20, -1)
    result = np.zeros((*sequences.shape[:-1], vectors.shape[-1]))
    for i in range(sites):
        result += weights[i][sequences[..., i]]
    return result - center @ vectors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out", type=Path)
    parser.add_argument("--threshold", type=float, default=.99)
    parser.add_argument("--discard", type=float, default=0,
                        help="explicit fraction of each saved trajectory to exclude")
    parser.add_argument("--rhat", type=float, default=1.01)
    parser.add_argument("--min-ess", type=float, default=400)
    parser.add_argument("--reference", type=Path,
                        help="reference align_A.msa path when inspecting results on another machine")
    args = parser.parse_args()
    if not 0 <= args.discard < 1 or not 0 < args.threshold < 1:
        parser.error("discard must be in [0,1), threshold in (0,1)")
    out = args.out
    parameters = tomllib.loads((out / "parameters.toml").read_text())
    runs = parameters["runs"]
    frames = [pd.read_csv(out / f"run_{i}_cold.csv") for i in range(1, runs + 1)]
    if len({len(f) for f in frames}) != 1:
        raise ValueError("runs must have equal numbers of chronological samples")
    discard = int(len(frames[0]) * args.discard)
    frames = [f.iloc[discard:].copy() for f in frames]
    if len(frames[0]) < 4:
        raise ValueError("need at least four retained draws per run")
    p = np.array([f.pnat.values for f in frames])
    seq = np.array([encode(f.sequence) for f in frames])
    keep = p > args.threshold
    if keep.sum() < 6 or np.any(keep.sum(axis=1) < 2):
        raise ValueError("too few sequences above the threshold for PCA/comparison")
    center, covariance = sequence_covariance(seq[keep])
    eigenvalues, eigenvectors = eigh(covariance, subset_by_index=[535, 539])
    vectors = eigenvectors[:, ::-1]
    projections = project_sequences(seq, center, vectors)

    # A fixed coordinate from the published reference alignment supplements PCA
    # fitted to the new sample, which could otherwise miss an unexplored family.
    ref = pd.read_csv(args.reference or parameters["reference_alignment"], sep=r"\s+", header=None,
                      names=["beta", "pnat", "sequence"])
    ref_center, ref_covariance = sequence_covariance(encode(ref.loc[ref.pnat > .995, "sequence"]))
    _, ref_vec = eigh(ref_covariance, subset_by_index=[539, 539])
    ref_projection = project_sequences(seq, ref_center, ref_vec)[..., 0]

    observables = {"pnat": diagnostic(p), "above_threshold": diagnostic(keep),
                   "reference_pc1": diagnostic(ref_projection)}
    for j in range(5):
        observables[f"pc{j+1}"] = diagnostic(projections[:, :, j])
        observables[f"selected_pc{j+1}"] = diagnostic(keep * projections[:, :, j])
    raw_freq = frequencies(seq.reshape(-1, seq.shape[-1]))
    common = np.flatnonzero((raw_freq > .01) & (raw_freq < .99))
    selected_common = np.flatnonzero((center > .01) & (center < .99))
    for j in common:
        observables[f"site_{j//20+1}_{AA[j%20]}"] = diagnostic(seq[:, :, j//20] == j%20)
    conditional_mcse = []
    # Estimate uncertainty of conditional frequencies using ratio-estimator
    # influence functions on the ORIGINAL time grid, not shortened traces.
    for j in selected_common:
        indicator = seq[:, :, j//20] == j%20
        observables[f"selected_site_{j//20+1}_{AA[j%20]}"] = diagnostic(keep * indicator)
        influence = keep * (indicator - center[j]) / keep.mean()
        conditional_mcse.append(float(az.mcse(influence, method="mean")))
    finite = {name: d for name, d in observables.items() if d["rhat"] is not None}
    failures = {name: d for name, d in finite.items()
                if not np.isfinite(d["rhat"]) or d["rhat"] >= args.rhat
                or not np.isfinite(d["bulk_ess"]) or d["bulk_ess"] < args.min_ess}
    summaries = [tomllib.loads((out / f"run_{i}_summary.toml").read_text())
                 for i in range(1, runs + 1)]
    trips = [sum(s["production_round_trips"]) for s in summaries]
    # These trip counts cover the full recorded production phase. With --discard
    # they are still reported, but cannot certify trips within the retained tail.
    passed = not failures and all(t > 0 for t in trips) and discard == 0
    report = {
        "threshold": args.threshold, "beta": parameters["betas"][-1],
        "runs": runs, "raw_samples_after_discard": int(p.size),
        "discarded_samples_per_run": discard, "retained_sequences": int(keep.sum()),
        "retained_by_run": keep.sum(axis=1).tolist(), "raw_mean_pnat": float(p.mean()),
        "retained_mean_pnat": float(p[keep].mean()), "maximum_pnat": float(p.max()),
        "raw_mean_by_run": p.mean(axis=1).tolist(),
        "first_half_mean_by_run": p[:, :p.shape[1]//2].mean(axis=1).tolist(),
        "last_half_mean_by_run": p[:, p.shape[1]//2:].mean(axis=1).tolist(),
        "production_round_trips": trips,
        "maximum_probability_error": max(s["maximum_probability_error"] for s in summaries),
        "rhat_limit": args.rhat, "minimum_bulk_ess_required": args.min_ess,
        "max_rhat": max(d["rhat"] for d in finite.values()),
        "min_bulk_ess": min(d["bulk_ess"] for d in finite.values()),
        "max_conditional_frequency_mcse": max(conditional_mcse),
        "observables_checked": len(finite), "failed_observables": failures,
        "observables": observables, "diagnostic_gates_passed": passed,
        "filtered_pca_eigenvalues": eigenvalues[::-1].tolist(),
        "arviz_version": az.__version__,
        "numpy_version": np.__version__, "pandas_version": pd.__version__,
        "diagnostic_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "interpretation": "Checks are necessary evidence, not a proof of global equilibrium. "
                          "A high folding probability alone does not establish mixing.",
    }
    # Allow infinity in a human-readable CSV; JSON uses null with a failure label.
    def clean(value):
        if isinstance(value, dict):
            return {k: clean(v) for k, v in value.items()}
        if isinstance(value, list):
            return [clean(v) for v in value]
        return None if isinstance(value, float) and not np.isfinite(value) else value
    (out / "diagnostics.json").write_text(json.dumps(clean(report), indent=2) + "\n")
    pd.DataFrame(observables).T.to_csv(out / "observables.csv", index_label="observable")
    # Export exactly the analyzed time window, still unfiltered. The figure's
    # analysis script applies and records the explicitly requested threshold.
    with (out / "analysis_alignment.msa").open("w") as stream:
        for frame in frames:
            for row in frame.itertuples():
                stream.write(f"{parameters['betas'][-1]} {row.pnat:.17g} {row.sequence}\n")
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    stride = parameters["thin_rounds"] * parameters["local_steps"]
    steps = frames[0]["sample"].values * stride
    for i in range(runs):
        axes[0].plot(steps, p[i], lw=.5, alpha=.6, label=f"Run {i+1}")
        axes[1].plot(steps, ref_projection[i], lw=.5, alpha=.6)
        axes[2].plot(steps, projections[i, :, 0], lw=.5, alpha=.6)
    axes[0].axhline(args.threshold, color="black", ls="--", lw=1)
    axes[0].legend(ncol=runs, fontsize=8)
    axes[0].set_ylabel("P_nat (unfiltered)")
    axes[1].set_ylabel("Reference PC1")
    axes[2].set_ylabel("Full-model PC1")
    axes[2].set_xlabel("Mutation proposals per replica after burn-in")
    fig.suptitle(f"Full-model parallel tempering, β={parameters['betas'][-1]:g}\n"
                 f"Monitored diagnostic gates: {'PASS' if passed else 'NOT MET'}")
    fig.savefig(out / "diagnostics.png", dpi=150)
    print(json.dumps({k: report[k] for k in (
        "raw_samples_after_discard", "retained_sequences", "raw_mean_pnat",
        "retained_mean_pnat", "production_round_trips", "max_rhat", "min_bulk_ess",
        "max_conditional_frequency_mcse", "diagnostic_gates_passed")}, indent=2))
    print(f"Failed observables: {len(failures)} of {len(finite)}; full results in {out}")


if __name__ == "__main__":
    main()
