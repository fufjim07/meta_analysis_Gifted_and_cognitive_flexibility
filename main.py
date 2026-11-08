from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from analysis_functions import (
    hedges_g,
    run_meta_analysis,
    cochran_q_test,
    calculate_heterogeneity,
    leave_one_out_analysis,
)
from visualization import create_forest_plot


def calculate_study(row):
    gifted_size = float(row["n מחוננים"])
    comparison_size = float(row["n השוואה"])
    direction = float(row["כיוון (+1 גבוה=טוב, -1 נמוך=טוב)"])

    if (
        not np.isfinite(gifted_size)
        or not np.isfinite(comparison_size)
        or gifted_size <= 1
        or comparison_size <= 1
        or direction not in (1, -1)
    ):
        raise ValueError(f"Invalid inputs for study: {row['מחקר']}")

    means_and_sds = [
        row["M מחוננים"],
        row["SD מחוננים"],
        row["M השוואה"],
        row["SD השוואה"],
    ]

    if all(pd.notna(value) for value in means_and_sds):
        if not np.all(np.isfinite(np.asarray(means_and_sds, dtype=float))):
            raise ValueError(f"Invalid means or SDs: {row['מחקר']}")

        result = hedges_g(
            n1=gifted_size,
            m1=float(row["M מחוננים"]),
            sd1=float(row["SD מחוננים"]),
            n2=comparison_size,
            m2=float(row["M השוואה"]),
            sd2=float(row["SD השוואה"]),
            direction=direction,
        )

        method = "Means and SDs"

    elif pd.notna(row["t (כשאין M/SD)"]):
        # The t sign must represent gifted minus comparison.
        t_value = float(row["t (כשאין M/SD)"])

        if not np.isfinite(t_value):
            raise ValueError(f"Invalid t statistic: {row['מחקר']}")

        degrees_of_freedom = gifted_size + comparison_size - 2
        correction_factor = 1 - 3 / (4 * degrees_of_freedom - 1)

        d_value = direction * t_value * np.sqrt(
            1 / gifted_size + 1 / comparison_size
        )

        result = {
            "d": d_value,
            "J": correction_factor,
            "g": correction_factor * d_value,
        }

        method = "t statistic"

    else:
        raise ValueError(f"Missing data for study: {row['מחקר']}")

    return pd.Series({
        "d_python": result["d"],
        "J_python": result["J"],
        "g_python": result["g"],
        "calculation_method": method,
    })


def main():
    # Locate files in the same folder as main.py
    base_folder = Path(__file__).resolve().parent
    input_path = base_folder / "1נתוני_חילוץ_מטא_אנליזה.xlsx"
    output_path = base_folder / "meta_analysis_results.xlsx"

    # Read input data
    study_data = pd.read_excel(
        input_path,
        sheet_name="M ו-SD (מדד ראשי)",
    )

    # Keep study rows and reset the index
    study_data = study_data[
        study_data["מחקר"].notna()
    ].copy().reset_index(drop=True)

    # Calculate an effect size for each study
    calculations = study_data.apply(calculate_study, axis=1)

    results = pd.concat(
        [
            study_data[["מחקר"]].rename(columns={"מחקר": "study"}),
            calculations,
        ],
        axis=1,
    )

    # Run all analyses
    study_results, summary = run_meta_analysis(study_data, results)

    q_results = cochran_q_test(study_results)

    heterogeneity_results = calculate_heterogeneity(
        q_results, summary
    )

    loo_results = leave_one_out_analysis(study_results)

    # Save all data and results to one Excel workbook
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        study_data.to_excel(
            writer, sheet_name="Input data", index=False
        )
        study_results.to_excel(
            writer, sheet_name="Study results", index=False
        )
        summary.to_excel(
            writer, sheet_name="Meta-analysis", index=False
        )
        q_results.to_excel(
            writer, sheet_name="Cochran Q", index=False
        )
        heterogeneity_results.to_excel(
            writer, sheet_name="Heterogeneity", index=False
        )
        loo_results.to_excel(
            writer, sheet_name="Leave-one-out", index=False
        )

    print("\nMeta-analysis summary:")
    print(summary.to_string(index=False))
    print(f"\nAll results saved to: {output_path}")
    forest_path = base_folder / "forest_plot.png"

    create_forest_plot(
        study_results=study_results,
        summary=summary,
        output_path=forest_path,
    )

    print(f"Forest plot saved to: {forest_path}")
   

    print(f"Forest plot saved to: {forest_path}")

if __name__ == "__main__":
    main()