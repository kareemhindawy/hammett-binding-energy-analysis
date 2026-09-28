#!/usr/bin/env python3
"""
Hammett correlation: binding energy (kcal/mol) vs. Hammett sigma
- Ordinary least-squares linear regression
- Publication-style scatter plot with labelled points, trendline, equation, R^2

Usage:   python hammett_binding_energy_analysis.py
Output:  hammett_binding_energy.png (300 dpi) + numeric results in the console
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy import stats

# ----------------------------------------------------------------------
# 1) Data - used exactly as provided (no rounding / no modification)
# ----------------------------------------------------------------------
names = ["H", "CH3", "OCH3", "Cl"]              # console output
tags = ["H", r"CH$_3$", r"OCH$_3$", "Cl"]       # plot labels
sigma = np.array([0.00, -0.17, -0.27, 0.23])    # Hammett sigma
be = np.array([-29.1, -43.2, -41.1, -35.5])     # Binding energy, kcal/mol

# ----------------------------------------------------------------------
# 2) Linear regression:  BE = slope * sigma + intercept
# ----------------------------------------------------------------------
fit = stats.linregress(sigma, be)
slope, intercept = fit.slope, fit.intercept
r2 = fit.rvalue ** 2
n = len(sigma)

t_crit = stats.t.ppf(0.975, df=n - 2)           # 95 % CI, n-2 degrees of freedom
slope_ci = (slope - t_crit * fit.stderr, slope + t_crit * fit.stderr)

predicted = slope * sigma + intercept
residuals = be - predicted

sgn = "-" if intercept < 0 else "+"
print("=" * 60)
print("Linear regression:  Binding energy = slope * sigma + intercept")
print("=" * 60)
print(f"n                = {n}  (df = {n - 2})")
print(f"slope            = {slope:.6f} kcal/mol per sigma unit  (SE = {fit.stderr:.4f})")
print(f"slope 95% CI     = [{slope_ci[0]:.2f}, {slope_ci[1]:.2f}]")
print(f"intercept        = {intercept:.6f} kcal/mol  (SE = {fit.intercept_stderr:.4f})")
print(f"R^2              = {r2:.6f}")
print(f"r (Pearson)      = {fit.rvalue:.6f}")
print(f"p-value (slope)  = {fit.pvalue:.4f}")
print(f"Equation         : y = {slope:.4f}x {sgn} {abs(intercept):.4f}")
print("-" * 60)
print(f"{'Subst.':<8}{'sigma':>8}{'BE_obs':>10}{'BE_fit':>10}{'Resid.':>10}")
for nm, s, y, f, r in zip(names, sigma, be, predicted, residuals):
    print(f"{nm:<8}{s:>8.2f}{y:>10.1f}{f:>10.3f}{r:>10.3f}")

# ----------------------------------------------------------------------
# 3) Figure
# ----------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "mathtext.default": "regular",
    "font.size": 11,
    "axes.linewidth": 1.1,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 5, "ytick.major.size": 5,
    "xtick.minor.size": 2.5, "ytick.minor.size": 2.5,
})

fig, ax = plt.subplots(figsize=(6.5, 4.8))

# trendline (extended slightly beyond the data range)
x_line = np.linspace(sigma.min() - 0.05, sigma.max() + 0.05, 200)
ax.plot(x_line, slope * x_line + intercept, color="#C0392B", lw=1.8,
        zorder=2, label="Linear fit")

# data points
ax.scatter(sigma, be, s=80, color="#1F4E79", edgecolor="black",
           linewidth=0.8, zorder=3, label="Computed data")

# point labels: (dx, dy) offset in points, horizontal / vertical alignment
offsets = {
    "H":    (0,  10, "center", "bottom"),
    "CH3":  (0, -10, "center", "top"),
    "OCH3": (-7,  8, "right",  "bottom"),
    "Cl":   (0, -10, "center", "top"),
}
for nm, tag, x, y in zip(names, tags, sigma, be):
    dx, dy, ha, va = offsets[nm]
    ax.annotate(tag, (x, y), xytext=(dx, dy), textcoords="offset points",
                ha=ha, va=va, fontsize=12,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85),
                zorder=4)

# equation / R^2 box
minus = "\u2212" if intercept < 0 else "+"
eq_text = (rf"$\mathit{{y}}$ = {slope:.4f}$\mathit{{x}}$ {minus} {abs(intercept):.4f}"
           "\n" rf"$\mathit{{R}}^{{2}}$ = {r2:.4f}"
           "\n" r"$\mathit{n}$ = 4")
ax.text(0.03, 0.97, eq_text, transform=ax.transAxes, ha="left", va="top",
        fontsize=11.5, linespacing=1.5,
        bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.5", lw=0.8))

# axes
ax.set_xlabel(r"Hammett substituent constant, $\sigma$", fontsize=12)
ax.set_ylabel("Binding energy (kcal/mol)", fontsize=12)
ax.set_xlim(-0.35, 0.32)
ax.set_ylim(-46, -26)
ax.xaxis.set_major_locator(MultipleLocator(0.1))
ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.5)
ax.legend(loc="lower right", frameon=True, fontsize=10, edgecolor="0.6")

fig.tight_layout()
fig.savefig("hammett_binding_energy.png", dpi=300, facecolor="white")
# fig.savefig("hammett_binding_energy.pdf")   # vector version (LaTeX / journals)
print("\nSaved: hammett_binding_energy.png")
