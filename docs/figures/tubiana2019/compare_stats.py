"""Quantitative comparison of the Fig. 4 statistics between two MSAs.
Usage: python compare_stats.py <stats dir 1> <stats dir 2> <output png>"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
d1, d2, out = sys.argv[1], sys.argv[2], sys.argv[3]
load = lambda d, n: np.loadtxt(f"{d}/{n}.csv", delimiter=",")
AA = "CMFILVWYAGTSNQDEHRKP"
r = lambda a, b: np.corrcoef(a.ravel(), b.ravel())[0, 1]
f1, f2 = load(d1, "freq"), load(d2, "freq")
c1, c2 = load(d1, "conservation"), load(d2, "conservation")
w1, w2 = load(d1, "pc1"), load(d2, "pc1")
C1, C2 = load(d1, "contact_cov"), load(d2, "contact_cov")
p1, p2 = load(d1, "pnat"), load(d2, "pnat")
pr1, pr2 = load(d1, "proj"), load(d2, "proj")
cos = abs(np.dot(w1.ravel(), w2.ravel())) / np.linalg.norm(w1) / np.linalg.norm(w2)
print(f"sequences (pnat>0.995): {len(p1)} vs {len(p2)}; mean pnat {p1.mean():.4f} vs {p2.mean():.4f}")
print(f"site frequencies: Pearson r = {r(f1, f2):.4f}, max abs difference = {np.abs(f1 - f2).max():.4f}")
print(f"conservation profile: Pearson r = {r(c1, c2):.4f}, mean abs difference = {np.abs(c1 - c2).mean():.4f}")
print(f"PC1 weights: |cosine| = {cos:.4f}")
print(f"contact covariance matrix: Pearson r = {r(C1, C2):.4f}")
fig, axs = plt.subplots(2, 2, figsize=(11, 9))
axs[0, 0].plot(range(1, 28), c1, "o-", label="Jacquin MSA"); axs[0, 0].plot(range(1, 28), c2, "s--", label="fresh sample")
axs[0, 0].set_title("conservation per site"); axs[0, 0].set_xlabel("site"); axs[0, 0].legend()
axs[0, 1].scatter(f1.ravel(), f2.ravel(), s=8); axs[0, 1].plot([0, 1], [0, 1], "k:")
axs[0, 1].set_title(f"site frequencies f_i(a), r = {r(f1, f2):.3f}"); axs[0, 1].set_xlabel("Jacquin MSA"); axs[0, 1].set_ylabel("fresh sample")
axs[1, 0].scatter(w1.ravel(), w2.ravel() * np.sign(np.dot(w1.ravel(), w2.ravel())), s=8)
axs[1, 0].set_title(f"PC1 coefficients, |cos| = {cos:.3f}"); axs[1, 0].set_xlabel("Jacquin MSA"); axs[1, 0].set_ylabel("fresh sample")
axs[1, 1].scatter(C1.ravel(), C2.ravel(), s=8); lim = max(np.abs(C1).max(), np.abs(C2).max()); axs[1, 1].plot([-lim, lim], [-lim, lim], "k:")
axs[1, 1].set_title(f"contact covariance C(a,b), r = {r(C1, C2):.3f}"); axs[1, 1].set_xlabel("Jacquin MSA"); axs[1, 1].set_ylabel("fresh sample")
fig.tight_layout(); fig.savefig(out, dpi=130); print("saved", out)
