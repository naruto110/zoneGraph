#!/usr/bin/env python3
"""Plot cache block threshold ablation."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["Skitter", "Pokec", "LiveJournal"]
THRESHOLDS = ["50K", "100K", "200K", "1M"]

LOAD_TOTAL = {
    "Skitter": np.array([6.35757, 5.36125, 5.37699, 2.09253]),
    "Pokec": np.array([31.1605, 30.2686, 25.8533, 12.9108]),
    "LiveJournal": np.array([68.8927, 67.2249, 62.3618, 45.7485]),
}

DUMP_COUNT = {
    "Skitter": np.array([4, 2, 2, 0]),
    "Pokec": np.array([15, 9, 6, 1]),
    "LiveJournal": np.array([26, 18, 11, 4]),
}

DUMP_TIME = {
    "Skitter": np.array([4.15649, 3.21872, 3.23771, 0.0]),
    "Pokec": np.array([22.7535, 22.1251, 18.1303, 6.50174]),
    "LiveJournal": np.array([45.7643, 44.592, 39.5339, 28.1591]),
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


def grouped_bars(ax, data, ylabel):
    x = np.arange(len(DATASETS))
    width = 0.18
    colors = ["#2F5D8C", "#D07A2D", "#5B8E3E", "#8C6BB1"]
    hatches = ["", "\\\\", "//", "xx"]
    handles = []

    for i, threshold in enumerate(THRESHOLDS):
        values = np.array([data[d][i] for d in DATASETS], dtype=float)
        bars = ax.bar(
            x + (i - 1.5) * width,
            values,
            width,
            label=threshold,
            color=colors[i],
            hatch=hatches[i],
            edgecolor="black",
            linewidth=0.45,
        )
        handles.append(bars[0])

    max_value = max(max(data[d]) for d in DATASETS)
    ax.set_ylim(0, max_value * 1.20 if max_value > 0 else 1)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x)
    ax.set_xticklabels(DATASETS, rotation=0, ha="center")
    ax.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.45)
    ax.set_axisbelow(True)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    return handles


def plot(output_dir: Path):
    set_paper_style()
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.85), sharex=False)

    handles = grouped_bars(axes[0], LOAD_TOTAL, "Load time (s)")
    axes[0].set_title("Load total", pad=3)

    grouped_bars(axes[1], DUMP_COUNT, "Dump count")
    axes[1].set_title("Dump count", pad=3)

    grouped_bars(axes[2], DUMP_TIME, "Dump time (s)")
    axes[2].set_title("Dump time", pad=3)

    fig.legend(
        handles,
        THRESHOLDS,
        title="Cache blocks",
        ncol=4,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.0,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.07, right=0.995, top=0.86, bottom=0.29, wspace=0.30)
    fig.savefig(output_dir / "ablation_cache_block_threshold.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "ablation_cache_block_threshold.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


if __name__ == "__main__":
    plot(Path("results/figures"))
