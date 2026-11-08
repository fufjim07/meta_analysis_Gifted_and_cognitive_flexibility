import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
import pandas as pd


def create_forest_plot(study_results, summary, output_path):
    effect_sizes = study_results["g_python"].to_numpy(dtype=float)
    lower_limits = study_results["ci_lower"].to_numpy(dtype=float)
    upper_limits = study_results["ci_upper"].to_numpy(dtype=float)
    study_names = study_results["study"].astype(str).tolist()

    pooled_g = float(summary["pooled_g"].iloc[0])
    pooled_lower = float(summary["ci_lower"].iloc[0])
    pooled_upper = float(summary["ci_upper"].iloc[0])
    tau_squared = float(summary["tau_squared"].iloc[0])

    # Random-effects weights control square sizes
    variances = study_results["variance"].to_numpy(dtype=float)
    weights = 1 / (variances + tau_squared)
    weight_percentages = 100 * weights / weights.sum()

    number_of_studies = len(study_names)
    positions = np.arange(number_of_studies, 0, -1)

    figure, axis = plt.subplots(
        figsize=(11, max(5, 0.45 * number_of_studies + 2))
    )

    # Individual study confidence intervals
    axis.hlines(
        positions,
        lower_limits,
        upper_limits,
        color="steelblue",
        linewidth=1.5,
    )

    axis.scatter(
        effect_sizes,
        positions,
        s=35 + 160 * weights / weights.max(),
        marker="s",
        color="steelblue",
        zorder=3,
    )

    # Pooled effect and its confidence interval
    diamond = Polygon(
        [
            (pooled_lower, 0),
            (pooled_g, 0.22),
            (pooled_upper, 0),
            (pooled_g, -0.22),
        ],
        closed=True,
        facecolor="darkred",
        edgecolor="darkred",
        zorder=4,
    )
    axis.add_patch(diamond)

    # Zero indicates no difference between groups
    axis.axvline(0, color="gray", linestyle="--", linewidth=1)
    axis.axhline(0.5, color="lightgray", linewidth=1)

    axis.set_yticks(list(positions) + [0])
    axis.set_yticklabels(study_names + ["Pooled effect"])
    axis.set_ylim(-0.7, number_of_studies + 1)

    minimum = min(lower_limits.min(), pooled_lower, 0)
    maximum = max(upper_limits.max(), pooled_upper, 0)
    padding = max(0.15, 0.1 * (maximum - minimum))
    axis.set_xlim(minimum - padding, maximum + padding)

    # Numeric results to the right of the plot
    for position, effect, lower, upper, weight in zip(
        positions,
        effect_sizes,
        lower_limits,
        upper_limits,
        weight_percentages,
    ):
        axis.text(
            1.03,
            position,
            f"{effect:.2f} [{lower:.2f}, {upper:.2f}]",
            transform=axis.get_yaxis_transform(),
            va="center",
            fontsize=9,
        )
        axis.text(
            1.57,
            position,
            f"{weight:.1f}%",
            transform=axis.get_yaxis_transform(),
            va="center",
            fontsize=9,
        )

    axis.text(
        1.03,
        0,
        f"{pooled_g:.2f} [{pooled_lower:.2f}, {pooled_upper:.2f}]",
        transform=axis.get_yaxis_transform(),
        va="center",
        color="darkred",
        fontweight="bold",
        fontsize=9,
    )

    axis.text(
        1.03,
        number_of_studies + 0.65,
        "Hedges' g [95% CI]",
        transform=axis.get_yaxis_transform(),
        fontweight="bold",
        fontsize=9,
    )
    axis.text(
        1.57,
        number_of_studies + 0.65,
        "Weight",
        transform=axis.get_yaxis_transform(),
        fontweight="bold",
        fontsize=9,
    )

    axis.set_xlabel(
        "Hedges' g\n"
        "Negative: comparison advantage | Positive: gifted advantage"
    )
    axis.set_title("Cognitive flexibility: random-effects meta-analysis")

    for border in ["top", "right", "left"]:
        axis.spines[border].set_visible(False)

    axis.tick_params(axis="y", length=0)
    figure.subplots_adjust(left=0.36, right=0.73, bottom=0.18, top=0.90)

    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)

    return output_path