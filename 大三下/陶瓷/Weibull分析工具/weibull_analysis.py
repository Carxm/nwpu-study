#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize


def parse_args():
    parser = argparse.ArgumentParser(
        description="Two-parameter Weibull analysis for specimen-level failure data."
    )
    parser.add_argument("input_csv", help="CSV file with at least one positive numeric value column.")
    parser.add_argument(
        "--value-column",
        default="value",
        help="Column containing specimen failure values such as strength or life. Default: value",
    )
    parser.add_argument(
        "--group-column",
        default="group",
        help="Optional grouping column for comparing batches or test conditions. Default: group",
    )
    parser.add_argument(
        "--facet-column",
        default=None,
        help="Optional column used to split the Weibull plot into subplots, such as scenario.",
    )
    parser.add_argument(
        "--label",
        default="Failure value",
        help="Axis label used in plots and summaries. Example: Flexural strength (MPa)",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=1000,
        help="Number of bootstrap resamples for confidence intervals. Default: 1000",
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=0.95,
        help="Confidence level for bootstrap intervals. Default: 0.95",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for bootstrap resampling. Default: 42",
    )
    parser.add_argument(
        "--outdir",
        default="weibull_output",
        help="Directory for CSV summary and plot. Default: weibull_output",
    )
    parser.add_argument(
        "--x-reference",
        type=float,
        default=None,
        help="Optional vertical reference line in x units, e.g. 1000 for Type X threshold.",
    )
    parser.add_argument(
        "--x-reference-label",
        default=None,
        help="Optional label shown next to the vertical reference line.",
    )
    args = parser.parse_args()
    if args.bootstrap < 50:
        parser.error("--bootstrap must be at least 50; use 1000 or more for final analysis.")
    if not 0 < args.confidence < 1:
        parser.error("--confidence must lie between 0 and 1.")
    return args


def _initial_guess(values):
    log_values = np.log(values)
    std_log = np.std(log_values, ddof=1) if len(values) > 1 else 0.2
    shape = max(0.2, min(20.0, 1.2 / max(std_log, 1e-6)))
    scale = float(np.median(values))
    return np.array([math.log(shape), math.log(scale)])


def fit_weibull_2p(values):
    values = np.asarray(values, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("Failure values must be finite numbers.")
    if np.any(values <= 0):
        raise ValueError("Weibull analysis requires strictly positive values.")
    if len(values) < 3:
        raise ValueError("The software requires at least 3 observations; this is not a statistical sample-size recommendation.")
    if np.ptp(values) == 0:
        raise ValueError("All observations are identical; Weibull scatter cannot be estimated.")

    def nll(theta):
        log_m, log_eta = theta
        m = math.exp(log_m)
        eta = math.exp(log_eta)
        scaled = values / eta
        power_term = np.exp(np.clip(m * np.log(scaled), -700, 700)).sum()
        return -(len(values) * math.log(m) - len(values) * m * math.log(eta) + (m - 1) * np.log(values).sum() - power_term)

    bounds = [
        (math.log(0.1), math.log(200.0)),
        (math.log(values.min() / 10.0), math.log(values.max() * 10.0)),
    ]
    result = minimize(nll, _initial_guess(values), method="L-BFGS-B", bounds=bounds)
    if not result.success:
        raise RuntimeError(f"Weibull fit failed: {result.message}")

    m = math.exp(result.x[0])
    eta = math.exp(result.x[1])
    return m, eta


def weibull_quantile(m, eta, failure_probability):
    return eta * (-math.log(1.0 - failure_probability)) ** (1.0 / m)


def bootstrap_fit(values, n_boot, confidence, rng):
    values = np.asarray(values, dtype=float)
    fits = []
    for _ in range(n_boot):
        sample = rng.choice(values, size=len(values), replace=True)
        try:
            fits.append(fit_weibull_2p(sample))
        except Exception:
            continue

    if len(fits) < max(50, n_boot // 5):
        raise RuntimeError("Bootstrap failed too often; inspect the dataset for duplicates or too few specimens.")

    fits = np.asarray(fits, dtype=float)
    alpha = 1.0 - confidence
    lo = alpha / 2.0
    hi = 1.0 - alpha / 2.0

    b10 = np.array([weibull_quantile(m, eta, 0.10) for m, eta in fits])
    b50 = np.array([weibull_quantile(m, eta, 0.50) for m, eta in fits])

    return {
        "fits": fits,
        "m_ci_low": float(np.quantile(fits[:, 0], lo)),
        "m_ci_high": float(np.quantile(fits[:, 0], hi)),
        "eta_ci_low": float(np.quantile(fits[:, 1], lo)),
        "eta_ci_high": float(np.quantile(fits[:, 1], hi)),
        "b10_ci_low": float(np.quantile(b10, lo)),
        "b10_ci_high": float(np.quantile(b10, hi)),
        "b50_ci_low": float(np.quantile(b50, lo)),
        "b50_ci_high": float(np.quantile(b50, hi)),
    }


def median_rank(n):
    i = np.arange(1, n + 1, dtype=float)
    return (i - 0.3) / (n + 0.4)


def weibull_probability_y(failure_probability):
    failure_probability = np.clip(failure_probability, 1e-9, 1.0 - 1e-9)
    return np.log(-np.log(1.0 - failure_probability))


def analyze_group(values, n_boot, confidence, rng):
    m, eta = fit_weibull_2p(values)
    intervals = bootstrap_fit(values, n_boot, confidence, rng)
    fits = intervals.pop("fits")
    return {
        "n": int(len(values)),
        "mean": float(np.mean(values)),
        "std": float(np.std(values, ddof=1)),
        "shape_m": float(m),
        "scale_eta": float(eta),
        "b10": float(weibull_quantile(m, eta, 0.10)),
        "b50": float(weibull_quantile(m, eta, 0.50)),
        "bootstrap_fits": fits,
        **intervals,
    }


def _failure_probability(value_grid, shape_m, scale_eta):
    return 1.0 - np.exp(-np.power(value_grid / scale_eta, shape_m))


def _display_name(group_name, facet_name):
    group_str = str(group_name)
    facet_str = str(facet_name)
    suffix = f"_{facet_str}"
    if group_str.endswith(suffix):
        group_str = group_str[: -len(suffix)]
    return group_str.replace("_", " ")


def make_probability_plot(grouped_values, fit_rows, outpath, label, facets=None, x_reference=None, x_reference_label=None, confidence=0.95):
    if facets:
        facet_names = list(facets.keys())
        fig, axes = plt.subplots(
            1,
            len(facet_names),
            figsize=(9 * len(facet_names), 8.5),
            sharey=True,
            gridspec_kw={"wspace": 0.08},
        )
        if len(facet_names) == 1:
            axes = [axes]
    else:
        facet_names = ["All data"]
        facets = {"All data": list(grouped_values.keys())}
        fig, ax = plt.subplots(figsize=(11, 8.5))
        axes = [ax]

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    yticks_p = np.array([0.01, 0.05, 0.10, 0.20, 0.50, 0.80, 0.95, 0.99])
    y_min = weibull_probability_y(np.array([0.005]))[0]
    y_max = weibull_probability_y(np.array([0.995]))[0]

    for ax, facet_name in zip(axes, facet_names):
        groups_in_facet = facets[facet_name]
        x_min = min(grouped_values[g].min() for g in groups_in_facet) * 0.85
        x_max = max(grouped_values[g].max() for g in groups_in_facet) * 1.12
        x_grid = np.geomspace(x_min, x_max, 300)
        text_y = 0.97

        for idx, group_name in enumerate(groups_in_facet):
            values = np.sort(np.asarray(grouped_values[group_name], dtype=float))
            ranks = median_rank(len(values))
            fit = fit_rows[group_name]
            color = colors[idx % len(colors)]
            display_name = _display_name(group_name, facet_name)

            ax.scatter(
                values,
                weibull_probability_y(ranks),
                s=34,
                color=color,
                label=f"{display_name} data",
                zorder=3,
            )

            failure_probability = _failure_probability(x_grid, fit["shape_m"], fit["scale_eta"])
            ax.plot(
                x_grid,
                weibull_probability_y(failure_probability),
                color=color,
                linewidth=2.0,
                label=f"{display_name} fit",
                zorder=2,
            )

            bootstrap_fits = fit["bootstrap_fits"]
            band_probs = []
            for shape_m, scale_eta in bootstrap_fits:
                band_probs.append(_failure_probability(x_grid, shape_m, scale_eta))
            band_probs = np.asarray(band_probs)
            alpha = (1.0 - confidence) / 2.0
            lower = np.quantile(band_probs, alpha, axis=0)
            upper = np.quantile(band_probs, 1.0 - alpha, axis=0)
            ax.fill_between(
                x_grid,
                weibull_probability_y(lower),
                weibull_probability_y(upper),
                color=color,
                alpha=0.14,
                linewidth=0,
                zorder=1,
            )

            ax.text(
                0.03,
                text_y,
                f"{display_name}: n={fit['n']}, m={fit['shape_m']:.2f}, "
                f"eta={fit['scale_eta']:.1f}, B10={fit['b10']:.1f}",
                transform=ax.transAxes,
                fontsize=9,
                va="top",
                ha="left",
                color=color,
                bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": color, "pad": 3},
            )
            text_y -= 0.12

        ax.axhline(
            weibull_probability_y(np.array([0.632]))[0],
            color="gray",
            linestyle="--",
            linewidth=1.0,
            alpha=0.8,
        )
        ax.text(
            0.98,
            weibull_probability_y(np.array([0.632]))[0],
            "63.2%",
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            color="gray",
            fontsize=9,
        )

        ax.set_xscale("log")
        if x_reference is not None:
            ax.axvline(x_reference, color="black", linestyle=":", linewidth=1.2, alpha=0.9)
            if x_reference_label:
                ax.text(
                    x_reference,
                    0.02,
                    x_reference_label,
                    transform=ax.get_xaxis_transform(),
                    rotation=90,
                    ha="right",
                    va="bottom",
                    fontsize=9,
                    color="black",
                )
        ax.set_yticks(weibull_probability_y(yticks_p))
        ax.set_yticklabels(["1%", "5%", "10%", "20%", "50%", "80%", "95%", "99%"])
        ax.set_ylim(y_min, y_max)
        ax.set_xlabel(label)
        ax.set_title(f"{facet_name} (2-parameter Weibull)")
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
        ax.tick_params(axis="y", labelsize=10)
        ax.tick_params(axis="x", labelsize=11)

    axes[0].set_ylabel("Failure probability")
    fig.suptitle(f"Weibull Probability Plot with {confidence:.0%} Pointwise Bootstrap Bands", fontsize=18, y=0.98)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.12, top=0.88, wspace=0.08)
    fig.savefig(outpath, dpi=220)
    plt.close(fig)


def main():
    args = parse_args()
    rng = np.random.default_rng(args.seed)

    df = pd.read_csv(args.input_csv)
    if args.value_column not in df.columns:
        raise KeyError(f"Missing value column: {args.value_column}")

    work = df.copy()
    if args.group_column not in work.columns:
        work[args.group_column] = "all"
    if args.facet_column and args.facet_column not in work.columns:
        raise KeyError(f"Missing facet column: {args.facet_column}")

    keep_cols = [args.group_column, args.value_column]
    if args.facet_column:
        keep_cols.append(args.facet_column)
    work = work[keep_cols].copy()
    numeric = pd.to_numeric(work[args.value_column], errors="coerce")
    invalid = work.isna().any(axis=1) | ~np.isfinite(numeric) | (numeric <= 0)
    if invalid.any():
        raise ValueError(f"Found {int(invalid.sum())} invalid row(s). Values must be positive, finite numbers and group/facet labels must be present. No rows were silently removed.")
    work[args.value_column] = numeric
    if work.empty:
        raise ValueError("The CSV contains no observations.")
    if args.facet_column:
        facet_counts = work.groupby(args.group_column)[args.facet_column].nunique()
        if (facet_counts > 1).any():
            raise ValueError("Each group must belong to one facet. Use unique group labels across facets.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    summaries = []
    grouped_values = {}
    fit_rows = {}
    facet_map = {}

    for group_name, group_df in work.groupby(args.group_column):
        values = group_df[args.value_column].to_numpy(dtype=float)
        summary = analyze_group(values, args.bootstrap, args.confidence, rng)
        summary["group"] = group_name
        summaries.append(summary)
        grouped_values[group_name] = values
        fit_rows[group_name] = summary
        if args.facet_column and args.facet_column in group_df.columns:
            facet_name = str(group_df[args.facet_column].iloc[0])
            facet_map.setdefault(facet_name, []).append(group_name)

    summary_df = pd.DataFrame(summaries)[
        [
            "group",
            "n",
            "mean",
            "std",
            "shape_m",
            "m_ci_low",
            "m_ci_high",
            "scale_eta",
            "eta_ci_low",
            "eta_ci_high",
            "b10",
            "b10_ci_low",
            "b10_ci_high",
            "b50",
            "b50_ci_low",
            "b50_ci_high",
        ]
    ].sort_values("group")

    summary_path = outdir / "weibull_summary.csv"
    plot_path = outdir / "weibull_probability_plot.png"
    summary_df.to_csv(summary_path, index=False)
    plot_fit_rows = {
        key: value for key, value in fit_rows.items()
    }
    make_probability_plot(
        grouped_values,
        plot_fit_rows,
        plot_path,
        args.label,
        facets=facet_map if facet_map else None,
        x_reference=args.x_reference,
        x_reference_label=args.x_reference_label,
        confidence=args.confidence,
    )

    print(f"Saved summary: {summary_path}")
    print(f"Saved plot: {plot_path}")
    print()
    print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
