import numpy as np
import matplotlib.pyplot as plt
from math import log

# --- Parameters ---
N = 340_000_000  # U.S. population
I1 = 100         # initial infections
r = 0.016        # weekly exponential growth rate (HIV-like pathogen)

# --- Time and outbreak dynamics ---
weeks = np.arange(1, 10 * 52 + 1)  # 10 years in weeks
years = weeks / 52.0
weekly_incidence = I1 * np.exp(r * (weeks - 1))
cumulative = np.cumsum(weekly_incidence)

# --- Cumulative incidence milestones ---
ci_targets = [0.0001, 0.001, 0.01, 0.1]  # 0.01%, 0.1%, 1%, 10%
ci_labels = ["0.01%", "0.1%", "1%", "10%"]

milestones = []
for ci in ci_targets:
    threshold = N * ci
    idx = np.argmax(cumulative >= threshold)
    if idx > 0 or cumulative[0] >= threshold:
        milestones.append((years[idx], threshold))
    else:
        milestones.append((None, threshold))

# --- 1 million infection point ---
one_million = 1_000_000
idx_million = np.argmax(cumulative >= one_million)
year_million = years[idx_million]
ci_million = (cumulative[idx_million] / N) * 100

# --- Plot ---
plt.figure(figsize=(10, 5))
plt.plot(years, cumulative, linewidth=2, color='navy')
plt.xlabel("Years since start")
plt.ylabel("Cumulative infections")
plt.title("Modeled cumulative infections (HIV-like pathogen; r=0.016/week)", pad=30)
plt.xlim(0, 10)
plt.xticks(np.arange(0, 11, 1))
plt.grid(True, linestyle="--", linewidth=0.5)

# --- Red milestone lines ---
for (x, y), label in zip(milestones, ci_labels):
    if x is not None:
        infections = int(y)
        plt.axvline(x=x, color='red', linestyle='--', linewidth=1)
        plt.text(
            x, cumulative.max() * 1.03,
            f"{label}\n{infections:,} infections\n({x:.2f} yrs)",
            color='red', fontsize=7, ha='center', va='bottom'
        )

# --- Gray 1M infections line ---
plt.axvline(x=year_million, color='gray', linestyle='--', linewidth=1.5)
plt.text(
    year_million, cumulative.max() * 1.03,
    f"1M infections\n({ci_million:.3f}% CI)\n({year_million:.2f} yrs)",
    color='gray', fontsize=7, ha='center', va='bottom'
)

plt.tight_layout()
plt.show()

