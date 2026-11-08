import numpy as np
import pandas as pd
from statsmodels.stats.meta_analysis import combine_effects
from scipy.stats import t as t_distribution,chi2



def hedges_g(n1, m1, sd1, n2, m2, sd2, direction=1):
    """this function calculates hedges g effect size for two independent groups, 
    given their sample sizes, means, and standard deviations.
    It also includes a correction factor for small sample sizes."""
    # input validation
    if n1 <= 1 or n2 <= 1:
        raise ValueError("in each group there must be at least two participants") 
    if sd1 < 0 or sd2 < 0:
        raise ValueError("standard deviations cannot be negative") 
    if direction not in (1, -1):
        raise ValueError("direction must be 1 or -1")

    # freedom degrees
    df = n1 + n2 - 2

    # סטיית תקן משותפת
    pooled_sd = np.sqrt(
        ((n1 - 1) * sd1**2 + (n2 - 1) * sd2**2) / df
    )

    if pooled_sd == 0:
        raise ValueError("pooled standard deviation is zero")

    # mean difference standardized
    d = direction * (m1 - m2) / pooled_sd

    # fix for small sample sizes
    J = 1 - 3 / (4 * df - 1)

    # effect size fixed for small sample sizes
    g = J * d

    return {
        "pooled_sd": pooled_sd,
        "d": d,
        "J": J,
        "g": g
    }

def run_meta_analysis(study_data, results):
    """this function calculates the meta-analysis results for the given study data and effect sizes.
    It computes the variances, standard errors, confidence intervals, and performs a random-effects model"""
    # Copy to avoid changing the original table
    study_results = results.copy()

    gifted_sizes = study_data["n מחוננים"].to_numpy(dtype=float)
    comparison_sizes = study_data["n השוואה"].to_numpy(dtype=float)

    degrees_of_freedom = gifted_sizes + comparison_sizes - 2

    d_values = study_results["d_python"].to_numpy(dtype=float)
    j_values = study_results["J_python"].to_numpy(dtype=float)
    g_values = study_results["g_python"].to_numpy(dtype=float)

    # Approximate sampling variance for independent groups
    variances = j_values**2 * (
        1 / gifted_sizes
        + 1 / comparison_sizes
        + d_values**2 / (2 * degrees_of_freedom)
    )

    if (
        len(g_values) < 2
        or not np.all(np.isfinite(g_values))
        or not np.all(np.isfinite(variances))
        or np.any(variances <= 0)
    ):
        raise ValueError(
            "At least two valid effects with positive variances are required"
        )

    study_results["variance"] = variances
    study_results["standard_error"] = np.sqrt(variances)
    study_results["ci_lower"] = (
        g_values - 1.96 * study_results["standard_error"]
    )
    study_results["ci_upper"] = (
        g_values + 1.96 * study_results["standard_error"]
    )

    # Random-effects model
    meta_result = combine_effects(
        effect=g_values,
        variance=variances,
        method_re="iterated",
        use_t=True
    )

    pooled_g = meta_result.mean_effect_re
    pooled_se = np.sqrt(meta_result.var_hksj_re)

    meta_df = len(g_values) - 1
    critical_t = t_distribution.ppf(0.975, df=meta_df)

    if pooled_se <= 0 or not np.isfinite(pooled_se):
        raise ValueError("The pooled standard error is invalid")

    t_value = pooled_g / pooled_se
    p_value = 2 * t_distribution.sf(abs(t_value), df=meta_df)

    summary = pd.DataFrame([{
        "method": "Paule-Mandel + Hartung-Knapp",
        "number_of_studies": len(g_values),
        "pooled_g": pooled_g,
        "standard_error": pooled_se,
        "ci_lower": pooled_g - critical_t * pooled_se,
        "ci_upper": pooled_g + critical_t * pooled_se,
        "degrees_of_freedom": meta_df,
        "t_value": t_value,
        "p_value": p_value,
        "tau_squared": meta_result.tau2
    }])

    return study_results, summary

def cochran_q_test(study_results):
    """this function calculates the Cochran's Q tes
      for heterogeneity among the effect sizes."""
    effect_sizes = study_results["g_python"].to_numpy(dtype=float)
    variances = study_results["variance"].to_numpy(dtype=float)

    if (
        len(effect_sizes) < 2
        or not np.all(np.isfinite(effect_sizes))
        or not np.all(np.isfinite(variances))
        or np.any(variances <= 0)
    ):
        raise ValueError(
            "At least two valid effects with positive variances are required"
        )

    # Inverse-variance weights, without between-study variance
    weights = 1 / variances

    weighted_mean = np.sum(weights * effect_sizes) / np.sum(weights)

    q_value = np.sum(weights * (effect_sizes - weighted_mean)**2)

    degrees_of_freedom = len(effect_sizes) - 1
    p_value = chi2.sf(q_value, df=degrees_of_freedom)

    return pd.DataFrame([{
        "Q": q_value,
        "degrees_of_freedom": degrees_of_freedom,
        "p_value": p_value
    }])

def calculate_heterogeneity(q_results, summary):
    """this function calculates the heterogeneity statistics (I squared and tau squared)"""
    q_value = float(q_results["Q"].iloc[0]) 
    degrees_of_freedom = int(
        q_results["degrees_of_freedom"].iloc[0]
    )
    tau_squared = float(summary["tau_squared"].iloc[0])

    if (
        not np.isfinite(q_value)
        or q_value < 0
        or degrees_of_freedom < 1
        or not np.isfinite(tau_squared)
        or tau_squared < 0
    ):
        raise ValueError("Invalid heterogeneity inputs")

    # Q-based I squared, expressed as a percentage
    if q_value == 0:
        i_squared = 0.0
    else:
        i_squared = max(
            0.0,
            (q_value - degrees_of_freedom) / q_value
        ) * 100

    return pd.DataFrame([{
        "I_squared_percent": i_squared,
        "tau_squared": tau_squared
    }])

def leave_one_out_analysis(study_results):
    """this function performs a leave-one-out analysis for the meta-analysis results.
    It recalculates the pooled effect size and confidence intervals after omitting each study in turn"""
    if len(study_results) < 3:
        raise ValueError("At least three studies are required")

    effect_sizes = study_results["g_python"].to_numpy(dtype=float)
    variances = study_results["variance"].to_numpy(dtype=float)
    study_names = study_results["study"].to_numpy()

    if (
        not np.all(np.isfinite(effect_sizes))
        or not np.all(np.isfinite(variances))
        or np.any(variances <= 0)
    ):
        raise ValueError("Effects must be finite and variances positive")

    # Full model for comparison
    full_model = combine_effects(
        effect=effect_sizes,
        variance=variances,
        method_re="iterated",
        use_t=True
    )

    full_pooled_g = full_model.mean_effect_re
    output_rows = []

    for omitted_index, study_name in enumerate(study_names):
        remaining_mask = np.arange(len(effect_sizes)) != omitted_index

        remaining_effects = effect_sizes[remaining_mask]
        remaining_variances = variances[remaining_mask]

        # Re-estimate between-study variance after each omission
        reduced_model = combine_effects(
            effect=remaining_effects,
            variance=remaining_variances,
            method_re="iterated",
            use_t=True
        )

        pooled_g = reduced_model.mean_effect_re
        pooled_se = np.sqrt(reduced_model.var_hksj_re)

        if not np.isfinite(pooled_se) or pooled_se <= 0:
            raise ValueError(
                f"Invalid pooled standard error after omitting: {study_name}"
            )

        degrees_of_freedom = len(remaining_effects) - 1
        critical_t = t_distribution.ppf(
            0.975, df=degrees_of_freedom
        )

        t_value = pooled_g / pooled_se
        p_value = 2 * t_distribution.sf(
            abs(t_value), df=degrees_of_freedom
        )

        output_rows.append({
            "omitted_study": study_name,
            "remaining_studies": len(remaining_effects),
            "full_pooled_g": full_pooled_g,
            "pooled_g": pooled_g,
            "change_in_g": pooled_g - full_pooled_g,
            "ci_lower": pooled_g - critical_t * pooled_se,
            "ci_upper": pooled_g + critical_t * pooled_se,
            "p_value": p_value,
            "tau_squared": reduced_model.tau2
        })

    return pd.DataFrame(output_rows)