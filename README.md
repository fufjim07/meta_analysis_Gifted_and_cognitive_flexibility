# Cognitive Flexibility Meta-Analysis

A Python project for a quantitative meta-analysis comparing cognitive flexibility in gifted and comparison groups. The analysis calculates individual study effect sizes, pools them using a random-effects model, assesses heterogeneity, and evaluates sensitivity to individual studies.

## Features

- Hedges' g for two independent groups, using sample sizes, means, and standard deviations.
- Effect-size calculation from a signed independent-groups t statistic when means and standard deviations are unavailable.
- Direction alignment so positive effects indicate an advantage for the gifted group.
- Approximate sampling variances and 95% confidence intervals for individual effects.
- Random-effects pooling using Paule-Mandel estimation with Hartung-Knapp inference.
- Cochran's Q, Q-based I squared, and tau squared.
- Leave-one-out sensitivity analysis with the model re-estimated after each omission.
- Excel export and a forest plot with study effects, confidence intervals, and the pooled effect.

## Project files

| File | Purpose |
| --- | --- |
| `main.py` | Reads the input workbook, runs the analyses, and saves results. |
| `analysis_functions.py` | Contains effect-size, meta-analysis, heterogeneity, sensitivity, and forest-plot functions. |
| `1נתוני_חילוץ_מטא_אנליזה.xlsx` | Input workbook. |
| `meta_analysis_results.xlsx` | Generated results workbook. |
| `forest_plot.png` | Generated forest plot. |

Keep both Python files and the input workbook in the same folder. The forest plot requires `create_forest_plot` to be defined in `analysis_functions.py`, imported into `main.py`, and called inside `main()`.

## Installation

Install the required packages in the Python environment used to run the project:

```bash
python -m pip install numpy pandas scipy statsmodels openpyxl matplotlib
```

## Run

From the project folder:

```bash
python main.py
```

If your Python executable has a different name or path, use that same executable for installation and execution.

The script saves its outputs next to `main.py`. Close an existing results workbook in Excel before running the script again, as Windows may prevent overwriting an open file. An existing output with the same name will be replaced.

## Input format

The script reads the worksheet named `M ו-SD (מדד ראשי)`. Preserve these column names:

| Column | Description |
| --- | --- |
| `מחקר` | Study name. Rows without a study name are excluded. |
| `n מחוננים` | Gifted-group sample size. |
| `M מחוננים` | Gifted-group mean. |
| `SD מחוננים` | Gifted-group standard deviation. |
| `n השוואה` | Comparison-group sample size. |
| `M השוואה` | Comparison-group mean. |
| `SD השוואה` | Comparison-group standard deviation. |
| `כיוון (+1 גבוה=טוב, -1 נמוך=טוב)` | Use `1` when higher scores are better and `-1` when lower scores are better. |
| `t (כשאין M/SD)` | Signed t statistic used when group means or standard deviations are missing. |

Use numeric cells for numeric inputs. Standard errors must be converted to standard deviations before entering them as SDs. For the t-based calculation, the t statistic must represent gifted minus comparison before applying the direction multiplier, and must come from an ordinary independent-groups t test with pooled variance. A reported absolute t value requires checking the original direction of the difference.

Each analysis row should represent one independent study effect. Multiple outcomes from the same participants require an approach that accounts for their dependence.

## Analysis workflow

1. Read and filter the study data.
2. Calculate Cohen's d, the small-sample correction J, and Hedges' g for each study.
3. Calculate approximate sampling variances and individual confidence intervals.
4. Fit the random-effects model and calculate the pooled effect and its confidence interval.
5. Calculate Cochran's Q, I squared, and tau squared.
6. Repeat the meta-analysis while omitting each study in turn.
7. Save the result tables and forest plot.

The current model uses **Paule-Mandel with Hartung-Knapp**, not REML. Results may differ from analyses using another estimator or another sampling-variance formula. Individual study intervals use a normal approximation; the pooled interval uses Hartung-Knapp inference with a t distribution.

## Excel output

| Worksheet | Contents |
| --- | --- |
| `Input data` | Study rows used in the analysis. |
| `Study results` | Individual effect sizes, correction factors, variances, standard errors, and confidence intervals. |
| `Meta-analysis` | Pooled Hedges' g, confidence interval, t statistic, p value, and tau squared. |
| `Cochran Q` | Q statistic, degrees of freedom, and p value. |
| `Heterogeneity` | I squared as a percentage and tau squared. |
| `Leave-one-out` | Pooled results after omitting each study and changes relative to the full model. |

## Reading the results

- Positive Hedges' g indicates an advantage for the gifted group after direction alignment; negative values indicate an advantage for the comparison group.
- A 95% confidence interval excluding zero corresponds to a two-sided test significant at the 0.05 level using the same inference method.
- Cochran's Q tests heterogeneity. A nonsignificant result does not establish that all true effects are identical.
- I squared estimates the proportion of observed variability attributed to heterogeneity rather than sampling error. It is not the percentage of studies that disagree.
- Tau squared estimates the variance of the true effects on the squared effect-size scale.
- Leave-one-out results indicate how much the pooled estimate changes when a study is omitted. Influence alone is not a reason to exclude a study.

In the forest plot, squares show study effects, horizontal lines show 95% confidence intervals, and the diamond shows the pooled effect and its confidence interval. Square areas reflect random-effects weights. The vertical line at zero represents no group difference.

## Scope

The project implements a study-level, independent-groups meta-analysis. It does not currently implement dependent effect sizes, paired pre/post effects, subgroup comparisons, meta-regression, or publication-bias analyses. Statistical results should be interpreted alongside study selection, outcome definitions, and study quality.
