#!/usr/bin/env python3
"""Plot GC strategy results as stacked end-to-end runtime bars."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["Skitter", "Pokec", "LiveJournal"]
POLICIES = ["Basic", "Greedy", "ATZ-GC"]

# For Skitter, the GC stress workload is the skew update phase.
# For Pokec and LiveJournal, we include both load and skew update phases.
TOTAL_TIME = {
    "Skitter": np.array([76.5065, 47.4261, 35.3489]),
    "Pokec": np.array([82.6449 + 104.118, 98.3830 + 5.02134, 79.0971 + 4.67991]),
    "LiveJournal": np.array([178.059 + 6.78123, 174.656 + 6.83587, 150.781 + 5.77257]),
}

GC_TIME = {
    "Skitter": np.array([70.6571, 42.2164, 31.4336]),
    "Pokec": np.array([24.7555 + 98.9529, 41.0163 + 0.0, 21.9039 + 0.0]),
    "LiveJournal": np.array([46.9262 + 0.0, 43.4211 + 0.0, 23.4232 + 0.0]),
}

GC_COUNT = {
    "Skitter": np.array([6 + 0, 4 + 0, 2 + 0]),
    "Pokec": np.array([3 + 6, 3 + 0, 1 + 0]),
    "LiveJournal": np.array([3 + 0, 2 + 0, 1 + 0]),
}

MOVED_BLOCKS = {
    "Skitter": np.array([1_252_379, 745_266, 562_029]),
    "Pokec": np.array([1_722_148, 0, 0]),
    "LiveJournal": np.array([0, 0, 0]),
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


def plot_runtime(output_dir: Path):
    colors = {
        "non_gc": "#8FB7D8",
        "gc": "#D07A2D",
    }

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.75), sharey=False)
    x = np.arange(len(POLICIES))
    width = 0.62
    legend_handles = None

    for ax, dataset in zip(axes, DATASETS):
        total = TOTAL_TIME[dataset]
        gc = GC_TIME[dataset]
        non_gc = total - gc

        bar_non_gc = ax.bar(
            x,
            non_gc,
            width,
            label="Non-GC time",
            color=colors["non_gc"],
            edgecolor="black",
            linewidth=0.45,
        )
        bar_gc = ax.bar(
            x,
            gc,
            width,
            bottom=non_gc,
            label="GC time",
            color=colors["gc"],
            edgecolor="black",
            linewidth=0.45,
            hatch="//",
        )

        ax.set_title(dataset, pad=3)
        ax.set_xticks(x)
        ax.set_xticklabels(POLICIES)
        ax.set_ylim(0, max(total) * 1.18)
        ax.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.45)
        ax.set_axisbelow(True)

        speedup = total[0] / total[2]
        ax.text(
            x[2],
            total[2] * 1.04,
            f"{speedup:.2f}x",
            ha="center",
            va="bottom",
            fontsize=9,
        )

        if ax is axes[0]:
            ax.set_ylabel("Runtime (s)")
            legend_handles = [bar_non_gc[0], bar_gc[0]]

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    fig.legend(
        legend_handles,
        ["Non-GC time", "GC time"],
        ncol=2,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.2,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.08, right=0.995, top=0.86, bottom=0.28, wspace=0.30)
    fig.savefig(output_dir / "gc_runtime_stacked.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "gc_runtime_stacked.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


def plot_reduction(output_dir: Path):
    """Plot normalized GC pressure metrics relative to Basic."""

    colors = ["#6A7FDB", "#59A14F", "#D07A2D"]
    hatches = ["", "\\\\", "//"]
    metrics = [
        ("Runtime", TOTAL_TIME),
        ("GC time", GC_TIME),
        ("GC count", GC_COUNT),
        ("Moved blocks", MOVED_BLOCKS),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(9.6, 2.65), sharey=True)
    x = np.arange(len(DATASETS))
    width = 0.23
    legend_handles = None

    for ax, (metric_name, metric_data) in zip(axes, metrics):
        handles = []
        for i, policy in enumerate(POLICIES):
            values = []
            for dataset in DATASETS:
                baseline = metric_data[dataset][0]
                value = metric_data[dataset][i]
                if baseline == 0:
                    values.append(np.nan)
                else:
                    values.append(value / baseline)
            values = np.array(values, dtype=float)
            bars = ax.bar(
                x + (i - 1) * width,
                np.nan_to_num(values, nan=0.0),
                width,
                label=policy,
                color=colors[i],
                hatch=hatches[i],
                edgecolor="black",
                linewidth=0.45,
            )
            handles.append(bars[0])

        ax.axhline(1.0, color="black", linewidth=0.6, linestyle=":")
        ax.set_title(metric_name, pad=3)
        ax.set_xticks(x)
        ax.set_xticklabels(DATASETS, rotation=18, ha="right")
        ax.set_ylim(0, 1.15)
        ax.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.45)
        ax.set_axisbelow(True)

        if ax is axes[0]:
            ax.set_ylabel("Normalized to Basic")
            legend_handles = handles

        if metric_name == "Moved blocks":
            ax.text(x[2], 0.05, "N/A", ha="center", va="bottom", rotation=90, fontsize=8)

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

    fig.legend(
        legend_handles,
        POLICIES,
        ncol=3,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.2,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.065, right=0.995, top=0.86, bottom=0.33, wspace=0.22)
    fig.savefig(output_dir / "gc_normalized_reduction.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "gc_normalized_reduction.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


def plot(output_dir: Path):
    set_paper_style()
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_runtime(output_dir)
    plot_reduction(output_dir)


if __name__ == "__main__":
    plot(Path("results/figures"))
