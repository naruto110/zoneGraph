#!/usr/bin/env python3
"""Plot degree-aware cache threshold ablation."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["Skitter", "Pokec", "LiveJournal"]
MIN_DEGREES = ["0", "8", "16"]

LOAD_TOTAL = {
    "Skitter": np.array([33.8094, 11.2306, 5.35002]),
    "Pokec": np.array([78.1386, 38.3663, 28.3023]),
    "LiveJournal": np.array([np.nan, 110.267, 67.507]),
}

DUMP_COUNT = {
    "Skitter": np.array([13, 5, 2]),
    "Pokec": np.array([14, 10, 9]),
    "LiveJournal": np.array([np.nan, 23, 18]),
}

DUMP_TIME = {
    "Skitter": np.array([26.7087, 8.22219, 3.20483]),
    "Pokec": np.array([55.688, 27.7204, 20.2135]),
    "LiveJournal": np.array([np.nan, 69.5317, 44.4088]),
}


def set_paper_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 11,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "hatch.linewidth": 0.5,
        }
    )


def grouped_bars(ax, data, ylabel, log=False):
    x = np.arange(len(DATASETS))
    width = 0.23
    colors = ["#2F5D8C", "#D07A2D", "#5B8E3E"]
    hatches = ["", "\\\\", "//"]
    handles = []

    for i, degree in enumerate(MIN_DEGREES):
        values = np.array([data[d][i] for d in DATASETS], dtype=float)
        bars = ax.bar(
            x + (i - 1) * width,
            np.nan_to_num(values, nan=0.0),
            width,
            label=f"min_degree={degree}",
            color=colors[i],
            hatch=hatches[i],
            edgecolor="black",
            linewidth=0.45,
        )
        handles.append(bars[0])

    finite_values = []
    for dataset in DATASETS:
        finite_values.extend([v for v in data[dataset] if not np.isnan(v)])

    if log:
        ax.set_yscale("log")
        ax.set_ylim(min(finite_values) * 0.65, max(finite_values) * 1.9)
    else:
        ax.set_ylim(0, max(finite_values) * 1.22)

    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(DATASETS, rotation=0, ha="center")
    ax.grid(axis="y", which="major", linestyle="--", linewidth=0.45, alpha=0.45)
    ax.set_axisbelow(True)

    if np.isnan(data["LiveJournal"][0]):
        y = min(finite_values) * 0.85 if log else max(finite_values) * 0.05
        ax.text(
            x[2] - width,
            y,
            "Killed",
            ha="center",
            va="bottom",
            rotation=90,
            fontsize=8,
        )

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    return handles


def plot(output_dir: Path):
    set_paper_style()
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.85), sharex=False)
    axes = axes.ravel()

    handles = grouped_bars(axes[0], LOAD_TOTAL, "Load time (s)", log=False)
    axes[0].set_title("Load total", pad=3)

    grouped_bars(axes[1], DUMP_COUNT, "Dump count", log=False)
    axes[1].set_title("Dump count", pad=3)

    grouped_bars(axes[2], DUMP_TIME, "Dump time (s)", log=False)
    axes[2].set_title("Dump time", pad=3)

    fig.legend(
        handles,
        [f"min_degree={d}" for d in MIN_DEGREES],
        ncol=3,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.2,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.07, right=0.995, top=0.86, bottom=0.28, wspace=0.30)
    fig.savefig(output_dir / "ablation_degree_aware.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "ablation_degree_aware.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


if __name__ == "__main__":
    plot(Path("results/figures"))
