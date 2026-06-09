# Statistical Analysis Code

This file documents the statistical analysis code included in `mednlp_data`.

## Archived Code Files

- `statistical_analysis_code/generate_answer_statistics_round_consistency.py`

## Input Files

- `../benchmark_queries/AnswerStatistics.xlsx`: adjudicated binary scoring table used for the main benchmark analysis.
- `../benchmark_queries/benchmark_queries_18.xlsx`: 18 benchmark query scenarios.
- `../benchmark_queries/candidates_results.xlsx`: returned candidate lists or model responses.
- `../curated_dataset/final_database.xlsx`: frozen curated MedNLP-Hub dataset.

## Output Files

- `../statistical_outputs/answer_statistics_round_consistency_stats.xlsx`

## Runtime Environment

Recommended environment:

- Python 3.10 or later
- pandas
- numpy
- openpyxl
- scipy
- statsmodels

Install dependencies:

```bash
pip install pandas numpy openpyxl scipy statsmodels
```

## Reproduction Commands

Run from `mednlp_data/statistical_analysis_code`:

```bash
python generate_answer_statistics_round_consistency.py
```

## Statistical Methods

The benchmark evaluation used 18 predefined query scenarios. Binary outcomes were coded as successful or unsuccessful. In the archived scoring workbook, `1` indicates success and `0` indicates failure.

For each system and round, the scripts calculate:

- Hit rate.
- Unconditional recency rate.
- Recency rate conditional on a successful hit.
- 95% Wilson confidence intervals for proportions.

Between-system paired comparisons use exact McNemar tests. Multiple pairwise comparisons are adjusted using the Benjamini-Hochberg false discovery rate procedure.

Round-to-round consistency is summarized using:

- Agreement rate.
- Cohen kappa for binary outcomes.
- Exact McNemar shift test.
- Benjamini-Hochberg adjusted q values for shift tests.

## Notes

The archived statistical code is intentionally limited to the main round-level and consistency analysis. It is intended to reproduce manuscript-level benchmark summaries from the frozen scoring file and is not dependent on the live MedNLP-Hub website.
