"""Draw the panels of Fig. 4 of Tubiana, Cocco & Monasson (2019) from the statistics
written by analyze_msa.jl.  Usage: python plot_fig4.py <stats dir> <output png> <title> [pnat threshold shown in labels]"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, TwoSlopeNorm
import logomaker

statsdir, out, title = sys.argv[1], sys.argv[2], sys.argv[3]
thr = sys.argv[4] if len(sys.argv) > 4 else "0.995"  # pnat threshold used by analyze_msa.jl, for the title
AA = "CMFILVWYAGTSNQDEHRKP"
classes = {"ED": ("negative", "#d62728"), "HKR": ("positive", "#1f77b4"), "NTSQ": ("polar", "#9467bd"),
           "FWY": ("aromatic", "#e6c700"), "ILMV": ("aliphatic", "black"), "C": ("cysteine", "#2ca02c"), "AGP": ("small", "#999999")}
color = {a: c for letters, (_, c) in classes.items() for a in letters}
cls_index = {a: k for k, (letters, _) in enumerate(classes.items()) for a in letters}

load = lambda name: np.loadtxt(f"{statsdir}/{name}.csv", delimiter=",")
freq, cons, pc1, proj, C, pnat = load("freq"), load("conservation"), load("pc1"), load("proj"), load("contact_cov"), load("pnat")
walks = {n: load(f"walk_{n}") for n in "AG"}
contacts = {n: load(f"contacts_{n}").astype(int) for n in "AG"}
seqs = [l.strip() for l in open(f"{statsdir}/sample_seqs.txt") if l.strip()]

fig = plt.figure(figsize=(16, 15))
gs = fig.add_gridspec(4, 3, height_ratios=[1.3, 0.8, 0.8, 1.2], hspace=0.45, wspace=0.3)

def draw_walk(ax, xyz, cm, label):
    ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], "-", color="#444444", lw=2.5)
    for i, j in cm:
        p, q = xyz[i - 1], xyz[j - 1]
        ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], ":", color="#bbbbbb", lw=1)
    ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], s=140, c=np.arange(27), cmap="viridis", depthshade=False, edgecolor="k")
    for k, (x, y, z) in enumerate(xyz, 1):
        ax.text(x, y, z + 0.12, str(k), fontsize=7, ha="center")
    ax.set_title(label); ax.set_xticks([0, 1, 2]); ax.set_yticks([0, 1, 2]); ax.set_zticks([0, 1, 2])
    ax.set_box_aspect((1, 1, 1))

draw_walk(fig.add_subplot(gs[0, 0], projection="3d"), walks["A"], contacts["A"], "(a) structure $S_A$ (chain solid, contacts dotted)")
draw_walk(fig.add_subplot(gs[0, 1], projection="3d"), walks["G"], contacts["G"], "(b) competing structure $S_G$ (ref. set #6863)")

ax = fig.add_subplot(gs[0, 2])
block = np.array([[cls_index[a] for a in s] for s in seqs])
ax.imshow(block, cmap=ListedColormap([c for _, (_, c) in classes.items()]), aspect="auto", interpolation="nearest")
ax.set_title(f"(c) {len(seqs)} sampled sequences (pnat > {thr})"); ax.set_xlabel("site"); ax.set_ylabel("sequence")
ax.set_xticks(range(0, 27, 5)); ax.set_xticklabels(range(1, 28, 5))
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=c, label=f"{n} ({l})") for l, (n, c) in classes.items()], fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1))

def logo(ax, mat, ttl, ylabel):
    df = pd.DataFrame(mat.T, columns=list(AA), index=range(1, 28))
    lg = logomaker.Logo(df, ax=ax, color_scheme=color, flip_below=False)
    ax.set_title(ttl); ax.set_ylabel(ylabel); ax.set_xticks(range(1, 28)); ax.set_xlabel("site")

logo(fig.add_subplot(gs[1, :]), freq * cons[None, :], "(d) position weight matrix (letter height ∝ frequency, total = conservation score)", "conservation")
logo(fig.add_subplot(gs[2, :]), pc1, "(e) first principal component PC1 (coefficients $w_i^{PC1}(v)$)", "PC1 weight")

ax = fig.add_subplot(gs[3, 0:2])
v = np.abs(C).max()
im = ax.imshow(C, cmap="RdBu_r", norm=TwoSlopeNorm(vcenter=0, vmin=-v, vmax=v))
ax.set_xticks(range(20)); ax.set_xticklabels(list(AA)); ax.set_yticks(range(20)); ax.set_yticklabels(list(AA))
for lab in ax.get_xticklabels() + ax.get_yticklabels():
    lab.set_color(color[lab.get_text()]); lab.set_fontweight("bold")
ax.set_title("(f) average covariance between amino acids on sites in contact in $S_A$")
fig.colorbar(im, ax=ax, fraction=0.03)

ax = fig.add_subplot(gs[3, 2])
ax.hist(proj, bins=120, color="#1f77b4")
ax.set_title("(g) projections of sequences on PC1"); ax.set_xlabel("projection on PC1"); ax.set_ylabel("count")

fig.suptitle(f"{title}\n{len(pnat)} sequences with pnat > {thr}; mean pnat = {pnat.mean():.4f}", fontsize=13)
fig.savefig(out, dpi=130, bbox_inches="tight")
print("saved", out)
