#!/usr/bin/env python3
"""Plot GC strategy evaluation results."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["Skitter", "Pokec", "LiveJournal"]
POLICIES = ["Basic", "Greedy", "ATZ-GC"]

LOAD_TOTAL = {
    "Skitter": np.array([np.nan, np.nan, np.nan]),
    "Pokec": np.array([82.6449, 98.3830, 79.0971]),
    "LiveJournal": np.array([178.059, 174.656, 150.781]),
}

LOAD_GC_TIME = {
    "Skitter": np.array([np.nan, np.nan, np.nan]),
    "Pokec": np.array([24.7555, 41.0163, 21.9039]),
    "LiveJournal": np.array([46.9262, 43.4211, 23.4232]),
}

LOAD_GC_COUNT = {
    "Skitter": np.array([np.nan, np.nan, np.nan]),
    "Pokec": np.array([3, 3, 1]),
    "LiveJournal": np.array([3, 2, 1]),
}

SKEW_TOTAL = {
    "Skitter": np.array([76.5065, 47.4261, 35.3489]),
    "Pokec": np.array([104.118, 5.02134, 4.67991]),
    "LiveJournal": np.array([6.78123, 6.83587, 5.77257]),
}

SKEW_GC_TIME = {
    "Skitter": np.array([70.6571, 42.2164, 31.4336]),
    "Pokec": np.array([98.9529, 0.0, 0.0]),
    "LiveJournal": np.array([0.0, 0.0, 0.0]),
}

SKEW_GC_COUNT = {
    "Skitter": np.array([6, 4, 2]),
    "Pokec": np.array([6, 0, 0]),
    "LiveJournal": np.array([0, 0, 0]),
}

MOVED_BLOCKS = {
    "Skitter": np.array([1_252_379, 745_266, 562_029]),
    "Pokec": np.array([1_722_148, 0, 0]),
    "LiveJournal": np.array([0, 0, 0]),
}

FOREGROUND = {
    "Skitter": np.array([3.88725, 3.24768, 2.03297]),
    "Pokec": np.array([3.09369, 2.91789, 2.72394]),
    "LiveJournal": np.array([3.91137, 3.95773, 2.78858]),
}


def set_paper_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 11,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "xtick.labelsize": 10,
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


def grouped_bars(ax, values_by_dataset, ylabel, ylim_pad=1.18, log=False):
    x = np.arange(len(DATASETS))
    width = 0.23
    colors = ["#6A7FDB", "#59A14F", "#D07A2D"]
    hatches = ["", "\\\\", "//"]

    handles = []
    for i, policy in enumerate(POLICIES):
        values = np.array([values_by_dataset[d][i] for d in DATASETS], dtype=float)
        bars = ax.bar(
            x + (i - 1) * width,
            values,
            width,
            label=policy,
            color=colors[i],
            hatch=hatches[i],
            edgecolor="black",
            linewidth=0.45,
        )
        handles.append(bars[0])

    finite_values = []
    for dataset in DATASETS:
        finite_values.extend([v for v in values_by_dataset[dataset] if not np.isnan(v)])

    if log:
        positive = [v for v in finite_values if v > 0]
        ax.set_yscale("log")
        ax.set_ylim(min(positive) * 0.6, max(positive) * 2.2)
    elif finite_values:
        ax.set_ylim(0, max(finite_values) * ylim_pad)

    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(DATASETS)
    ax.grid(axis="y", which="major", linestyle="--", linewidth=0.45, alpha=0.45)
    ax.set_axisbelow(True)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    return handles


def plot_total_time(output_dir: Path):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.75))

    handles = grouped_bars(axes[0], LOAD_TOTAL, "Load time (s)")
    axes[0].set_title("Load phase", pad=3)

    grouped_bars(axes[1], SKEW_TOTAL, "Skew update time (s)", log=True)
    axes[1].set_title("Skew update phase", pad=3)

    fig.legend(
        handles,
        POLICIES,
        ncol=3,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.1,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.075, right=0.995, top=0.86, bottom=0.28, wspace=0.32)
    fig.savefig(output_dir / "gc_total_time.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "gc_total_time.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


def plot_gc_details(output_dir: Path):
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 2.75))

    handles = grouped_bars(axes[0], SKEW_GC_COUNT, "Skew GC count")
    axes[0].set_title("GC invocations", pad=3)

    grouped_bars(axes[1], SKEW_GC_TIME, "Skew GC time (s)", log=True)
    axes[1].set_title("GC time", pad=3)

    moved_millions = {d: MOVED_BLOCKS[d] / 1_000_000.0 for d in DATASETS}
    grouped_bars(axes[2], moved_millions, "Moved blocks (M)")
    axes[2].set_title("Data movement", pad=3)

    fig.legend(
        handles,
        POLICIES,
        ncol=3,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.1,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.065, right=0.995, top=0.86, bottom=0.28, wspace=0.34)
    fig.savefig(output_dir / "gc_details.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "gc_details.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


def plot(output_dir: Path):
    set_paper_style()
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_total_time(output_dir)
    plot_gc_details(output_dir)


if __name__ == "__main__":
    plot(Path("results/figures"))
