import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.proportion import proportion_confint

BASE_DIR = Path(__file__).resolve().parents[1]
INPUT_XLSX = BASE_DIR / "benchmark_queries" / "AnswerStatistics.xlsx"
OUTPUT_XLSX = BASE_DIR / "statistical_outputs" / "answer_statistics_round_consistency_stats.xlsx"


def binarize(value):
    if pd.isna(value):
        return np.nan
    text = str(value).strip().lower()
    if text in {"yes", "y", "1", "true"}:
        return 1
    if text in {"x", "no", "n", "0", "false", "no."}:
        return 0
    return np.nan


def wilson_ci(k, n):
    if n <= 0:
        return np.nan, np.nan
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return float(lo), float(hi)


def pct(x):
    if pd.isna(x):
        return np.nan
    return float(x * 100.0)


def format_pct(x):
    if pd.isna(x):
        return "NA"
    return f"{x:.2f}%"


def bh_fdr(p_values):
    p = np.asarray(p_values, dtype=float)
    m = len(p)
    if m == 0:
        return np.array([])
    order = np.argsort(p)
    q = np.empty(m, dtype=float)
    for rank, idx in enumerate(order, start=1):
        q[idx] = p[idx] * m / rank
    for i in range(m - 2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.minimum(q, 1.0)


def cohen_kappa_binary(x, y):
    valid = x.notna() & y.notna()
    xv = x[valid].astype(int)
    yv = y[valid].astype(int)
    n = int(valid.sum())
    if n == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    n11 = int(((xv == 1) & (yv == 1)).sum())
    n10 = int(((xv == 1) & (yv == 0)).sum())
    n01 = int(((xv == 0) & (yv == 1)).sum())
    n00 = int(((xv == 0) & (yv == 0)).sum())

    po = (n11 + n00) / n
    px1 = (n11 + n10) / n
    py1 = (n11 + n01) / n
    px0 = 1.0 - px1
    py0 = 1.0 - py1
    pe = px1 * py1 + px0 * py0

    if math.isclose(1.0 - pe, 0.0):
        kappa = np.nan
    else:
        kappa = (po - pe) / (1.0 - pe)

    return po, pe, kappa, n10, n01


def add_unique_headers(df):
    counts = {}
    unique = []
    for col in df.columns:
        counts[col] = counts.get(col, 0) + 1
        unique.append(col if counts[col] == 1 else f"{col}__{counts[col]}")
    df.columns = unique
    return df


def build_model_map(df):
    columns = list(df.columns)

    def recency_after(hit_col):
        idx = columns.index(hit_col)
        if idx + 1 >= len(columns):
            raise ValueError(f"No recency column found after {hit_col}")
        return columns[idx + 1]

    model_map = {
        "GPT": {
            "hit_r1": "ChatGPT 5.2 Hit rate - round 1",
            "hit_r2": "ChatGPT 5.2 Hit rate - round 2",
        },
        "Claude": {
            "hit_r1": "claude-sonnet-4-6 round 1",
            "hit_r2": "claude-sonnet-4-6 round 2",
        },
        "Gemini": {
            "hit_r1": "gemini-2.5-flash round 1",
            "hit_r2": "gemini-2.5-flash round 2",
        },
        "PubMed": {
            "hit_r1": "PubMed round 1",
            "hit_r2": "PubMed round 2",
        },
    }

    for model in model_map:
        model_map[model]["rec_r1"] = recency_after(model_map[model]["hit_r1"])
        model_map[model]["rec_r2"] = recency_after(model_map[model]["hit_r2"])

    return model_map


def descriptive_rows(df, model_map):
    rows = []
    n_tasks = len(df)
    for model, cols in model_map.items():
        for round_id in [1, 2]:
            hit_col = cols[f"hit_r{round_id}"]
            rec_col = cols[f"rec_r{round_id}"]

            hit = df[hit_col]
            rec = df[rec_col]

            n_hit = int(hit.notna().sum())
            k_hit = int(hit.sum())
            hit_rate = (k_hit / n_hit) if n_hit > 0 else np.nan
            hit_lo, hit_hi = wilson_ci(k_hit, n_hit)

            n_rec = int(rec.notna().sum())
            k_rec = int(rec.sum())
            rec_rate = (k_rec / n_rec) if n_rec > 0 else np.nan
            rec_lo, rec_hi = wilson_ci(k_rec, n_rec)

            mask = (hit == 1) & rec.notna()
            n_cond = int(mask.sum())
            k_cond = int((rec[mask] == 1).sum()) if n_cond > 0 else 0
            cond_rate = (k_cond / n_cond) if n_cond > 0 else np.nan
            cond_lo, cond_hi = wilson_ci(k_cond, n_cond) if n_cond > 0 else (np.nan, np.nan)

            rows.append(
                {
                    "Model": model,
                    "Round": round_id,
                    "Tasks_N": n_tasks,
                    "Hit_Success": k_hit,
                    "Hit_N": n_hit,
                    "Hit_Rate": hit_rate,
                    "Hit_Rate_pct": pct(hit_rate),
                    "Hit_95CI_low_pct": pct(hit_lo),
                    "Hit_95CI_high_pct": pct(hit_hi),
                    "Hit_SCI_fmt": f"{k_hit}/{n_hit}, {format_pct(pct(hit_rate))} ({format_pct(pct(hit_lo))}-{format_pct(pct(hit_hi))})",
                    "Recency_Success_unconditional": k_rec,
                    "Recency_N_unconditional": n_rec,
                    "Recency_Rate_unconditional": rec_rate,
                    "Recency_Rate_unconditional_pct": pct(rec_rate),
                    "Recency_uncond_95CI_low_pct": pct(rec_lo),
                    "Recency_uncond_95CI_high_pct": pct(rec_hi),
                    "Recency_uncond_SCI_fmt": f"{k_rec}/{n_rec}, {format_pct(pct(rec_rate))} ({format_pct(pct(rec_lo))}-{format_pct(pct(rec_hi))})",
                    "Recency_Success_conditional": k_cond,
                    "Recency_N_conditional": n_cond,
                    "Recency_Rate_conditional": cond_rate,
                    "Recency_Rate_conditional_pct": pct(cond_rate),
                    "Recency_cond_95CI_low_pct": pct(cond_lo),
                    "Recency_cond_95CI_high_pct": pct(cond_hi),
                    "Recency_cond_SCI_fmt": f"{k_cond}/{n_cond}, {format_pct(pct(cond_rate))} ({format_pct(pct(cond_lo))}-{format_pct(pct(cond_hi))})",
                    "Hit_Column": hit_col,
                    "Recency_Column": rec_col,
                }
            )
    return pd.DataFrame(rows)


def pairwise_model_tests(df, model_map):
    rows = []
    p_hits = []
    p_recs = []

    for round_id in [1, 2]:
        models = list(model_map.keys())
        for model_a, model_b in combinations(models, 2):
            hit_a = df[model_map[model_a][f"hit_r{round_id}"]].astype(int)
            hit_b = df[model_map[model_b][f"hit_r{round_id}"]].astype(int)
            rec_a = df[model_map[model_a][f"rec_r{round_id}"]].astype(int)
            rec_b = df[model_map[model_b][f"rec_r{round_id}"]].astype(int)

            n00_h = int(((hit_a == 0) & (hit_b == 0)).sum())
            n01_h = int(((hit_a == 0) & (hit_b == 1)).sum())
            n10_h = int(((hit_a == 1) & (hit_b == 0)).sum())
            n11_h = int(((hit_a == 1) & (hit_b == 1)).sum())
            p_h = 1.0 if (n01_h + n10_h) == 0 else float(mcnemar([[n00_h, n01_h], [n10_h, n11_h]], exact=True).pvalue)

            rows.append(
                {
                    "Round": round_id,
                    "Metric": "Hit",
                    "Model_A": model_a,
                    "Model_B": model_b,
                    "Paired_N": len(df),
                    "n00_A0_B0": n00_h,
                    "n01_A0_B1": n01_h,
                    "n10_A1_B0": n10_h,
                    "n11_A1_B1": n11_h,
                    "McNemar_p_exact": p_h,
                }
            )
            p_hits.append(p_h)

            n00_r = int(((rec_a == 0) & (rec_b == 0)).sum())
            n01_r = int(((rec_a == 0) & (rec_b == 1)).sum())
            n10_r = int(((rec_a == 1) & (rec_b == 0)).sum())
            n11_r = int(((rec_a == 1) & (rec_b == 1)).sum())
            p_r = 1.0 if (n01_r + n10_r) == 0 else float(mcnemar([[n00_r, n01_r], [n10_r, n11_r]], exact=True).pvalue)

            rows.append(
                {
                    "Round": round_id,
                    "Metric": "Recency_unconditional",
                    "Model_A": model_a,
                    "Model_B": model_b,
                    "Paired_N": len(df),
                    "n00_A0_B0": n00_r,
                    "n01_A0_B1": n01_r,
                    "n10_A1_B0": n10_r,
                    "n11_A1_B1": n11_r,
                    "McNemar_p_exact": p_r,
                }
            )
            p_recs.append(p_r)

    out = pd.DataFrame(rows)
    mask_hit = out["Metric"] == "Hit"
    mask_rec = out["Metric"] == "Recency_unconditional"

    q_hit = bh_fdr(out.loc[mask_hit, "McNemar_p_exact"].tolist())
    q_rec = bh_fdr(out.loc[mask_rec, "McNemar_p_exact"].tolist())

    out.loc[mask_hit, "FDR_q"] = q_hit
    out.loc[mask_rec, "FDR_q"] = q_rec
    out["Significant_q<0.05"] = np.where(out["FDR_q"] < 0.05, "Yes", "No")
    return out.sort_values(["Metric", "FDR_q", "McNemar_p_exact", "Round", "Model_A", "Model_B"]).reset_index(drop=True)


def consistency_tests(df, model_map):
    rows = []
    for model, cols in model_map.items():
        for metric, c1, c2 in [
            ("Hit", cols["hit_r1"], cols["hit_r2"]),
            ("Recency_unconditional", cols["rec_r1"], cols["rec_r2"]),
        ]:
            x = df[c1]
            y = df[c2]
            po, pe, kappa, n10, n01 = cohen_kappa_binary(x, y)
            valid = x.notna() & y.notna()
            xv = x[valid].astype(int)
            yv = y[valid].astype(int)
            n = int(valid.sum())
            n11 = int(((xv == 1) & (yv == 1)).sum())
            n00 = int(((xv == 0) & (yv == 0)).sum())

            p_shift = 1.0 if (n10 + n01) == 0 else float(mcnemar([[n00, n01], [n10, n11]], exact=True).pvalue)

            rows.append(
                {
                    "Model": model,
                    "Metric": metric,
                    "N_paired": n,
                    "Agree_n": n11 + n00,
                    "Agreement_rate": po,
                    "Agreement_rate_pct": pct(po),
                    "Expected_agreement": pe,
                    "Cohen_kappa": kappa,
                    "Round1_yes_Round2_no": n10,
                    "Round1_no_Round2_yes": n01,
                    "McNemar_p_shift": p_shift,
                    "Consistency_Interpretation": "Good/Very good" if (pd.notna(kappa) and kappa >= 0.60) else "Moderate/Lower",
                    "Round1_Column": c1,
                    "Round2_Column": c2,
                }
            )
    out = pd.DataFrame(rows)
    out["FDR_q_shift"] = np.nan
    for metric in out["Metric"].unique():
        idx = out["Metric"] == metric
        q = bh_fdr(out.loc[idx, "McNemar_p_shift"].tolist())
        out.loc[idx, "FDR_q_shift"] = q
    out["Significant_shift_q<0.05"] = np.where(out["FDR_q_shift"] < 0.05, "Yes", "No")
    return out.sort_values(["Metric", "Model"]).reset_index(drop=True)


def build_sci_main_table(descriptive_df):
    table = descriptive_df.copy()
    table["Row"] = table["Model"] + " (Round " + table["Round"].astype(str) + ")"
    table = table[
        [
            "Row",
            "Hit_SCI_fmt",
            "Recency_uncond_SCI_fmt",
            "Recency_cond_SCI_fmt",
        ]
    ].rename(
        columns={
            "Row": "Model_Round",
            "Hit_SCI_fmt": "Hit rate (n/N, %, 95% CI)",
            "Recency_uncond_SCI_fmt": "Recency unconditional (n/N, %, 95% CI)",
            "Recency_cond_SCI_fmt": "Recency among hits (n/N, %, 95% CI)",
        }
    )
    return table


def main():
    data = pd.read_excel(INPUT_XLSX, sheet_name="05_Answer_Data")
    data = data.rename(columns={col: str(col).replace("Pubemd", "PubMed") for col in data.columns})
    data = add_unique_headers(data)

    metadata_cols = ["Prompt", "PMID", "Task type", "Text type", "Variable", "Language", "Modeling type"]
    for col in data.columns:
        if col not in metadata_cols and not pd.api.types.is_numeric_dtype(data[col]):
            data[col] = data[col].map(binarize)

    model_map = build_model_map(data)
    desc = descriptive_rows(data, model_map)
    pairwise = pairwise_model_tests(data, model_map)
    consistency = consistency_tests(data, model_map)
    sci_table = build_sci_main_table(desc)

    clean_cols = metadata_cols.copy()
    for model in ["GPT", "Claude", "Gemini", "PubMed"]:
        clean_cols += [model_map[model]["hit_r1"], model_map[model]["rec_r1"], model_map[model]["hit_r2"], model_map[model]["rec_r2"]]
    clean_cols = [c for c in clean_cols if c in data.columns]
    clean_data = data[clean_cols].copy()

    methods = pd.DataFrame(
        [
            {"Item": "Input", "Value": INPUT_XLSX},
            {"Item": "Output", "Value": OUTPUT_XLSX},
            {"Item": "Task count", "Value": len(data)},
            {"Item": "Systems", "Value": "PubMed, GPT, Gemini, Claude"},
            {"Item": "Experiments", "Value": "Two rounds (round 1 and round 2)"},
            {"Item": "Binary coding", "Value": "1/Yes/True=True; 0/No/False=False"},
            {"Item": "Proportion CI", "Value": "95% Wilson CI"},
            {"Item": "Between-model paired test", "Value": "Exact McNemar"},
            {"Item": "Consistency tests", "Value": "Agreement rate, Cohen kappa, McNemar shift test"},
            {"Item": "Multiple testing", "Value": "Benjamini-Hochberg FDR"},
        ]
    )

    OUTPUT_XLSX.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        methods.to_excel(writer, sheet_name="00_Methods", index=False)
        sci_table.to_excel(writer, sheet_name="01_SCI_Main_Table", index=False)
        desc.to_excel(writer, sheet_name="02_Descriptive_By_Round", index=False)
        consistency.to_excel(writer, sheet_name="03_Consistency_Tests", index=False)
        pairwise.to_excel(writer, sheet_name="04_Pairwise_McNemar", index=False)
        clean_data.to_excel(writer, sheet_name="05_Cleaned_Data", index=False)

    print(f"Wrote {OUTPUT_XLSX}")


if __name__ == "__main__":
    main()
