"""Entropy curve with the area criterion of Ochoa et al. (2021), linear axes.
Usage: python3 entropy_curve.py H.npy 70 600 out.png
  H.npy : entropies of all positions (any order)
  70    : threshold x (%)
  600   : last rank shown on the x-axis
"""
import sys, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

H = np.sort(np.load(sys.argv[1]))[::-1]          # non-increasing entropy
x = float(sys.argv[2]); xmax = int(sys.argv[3]); out = sys.argv[4]
r = np.arange(1, len(H) + 1)
tail = np.cumsum(H[::-1])[::-1]                   # tail[k] = area from rank k+1 onward
q = 1 + max(k for k in range(len(H)) if tail[k] >= x / 100 * H.sum())

DARK, LIGHT, LINE, INK, MUTED = "#1c5cab", "#cde2fb", "#2a78d6", "#2b2b2b", "#6b6b6b"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED})
fig, ax = plt.subplots(figsize=(3.5, 2.2))
ax.fill_between(r[:q], 0, H[:q], step="post", color=DARK, alpha=.85, lw=0)
ax.fill_between(r[q-1:], 0, H[q-1:], step="post", color=LIGHT, lw=0)
ax.step(r, H, where="post", color=LINE, lw=1.4)
ax.axvline(q, color=INK, lw=.9, ls="--")
ax.text(q + 0.015 * xmax, 0.9 * H[0] + 0.07, rf"$q = {q}$", color=INK, fontsize=8.5)
ax.annotate(rf"{x:.0f}% of the total area", xy=(3.2 * q, 0.25 * H[0] * 0.5),
            xytext=(4.8 * q, 0.45 * H[0]), fontsize=8, color=INK,
            arrowprops=dict(arrowstyle="->", color=INK, lw=.7))
ax.set_xlim(0, xmax); ax.set_ylim(0, 1.0)
ax.set_xlabel("Position rank"); ax.set_ylabel("Entropy")
ax.grid(axis="y", color="#e6e6e6", lw=.6); ax.set_axisbelow(True)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.tick_params(colors=MUTED, labelcolor=INK)
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
print("q =", q)
