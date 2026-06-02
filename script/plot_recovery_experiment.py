#!/usr/bin/env python3
"""Plot ZoneGraph recovery experiment results."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["Skitter", "Pokec", "LiveJournal"]
LOAD_TOTAL = np.array([3.40048, 22.8916, 44.4979])
RECOVERY_TOTAL = np.array([5.49613, 4.9163, 13.1166])
DUMP_COUNT = np.array([2, 7, 12])


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
        }
    )


def plot(output_dir: Path):
    set_paper_style()
    output_dir.mkdir(parents=True, exist_ok=True)

    x = np.arange(len(DATASETS))
    width = 0.30

    fig, ax = plt.subplots(figsize=(5.8, 2.75))

    load_bar = ax.bar(
        x - width / 2,
        LOAD_TOTAL,
        width,
        label="Load",
        color="#2F5D8C",
        edgecolor="black",
        linewidth=0.45,
    )
    recovery_bar = ax.bar(
        x + width / 2,
        RECOVERY_TOTAL,
        width,
        label="Recovery",
        color="#D07A2D",
        edgecolor="black",
        linewidth=0.45,
        hatch="//",
    )

    ax.set_ylabel("Time (s)")
    ax.set_xticks(x)
    ax.set_xticklabels(DATASETS, rotation=0, ha="center")
    ax.set_ylim(0, max(LOAD_TOTAL.max(), RECOVERY_TOTAL.max()) * 1.18)
    ax.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.45)
    ax.set_axisbelow(True)

    ax2 = ax.twinx()
    dump_line = ax2.plot(
        x,
        DUMP_COUNT,
        marker="o",
        markersize=4.5,
        linewidth=1.3,
        color="#5B8E3E",
        label="Dump count",
    )
    ax2.set_ylabel("Dump count")
    ax2.set_ylim(0, max(DUMP_COUNT) * 1.35)
    ax2.grid(False)

    for i, count in enumerate(DUMP_COUNT):
        ax2.text(
            x[i],
            count + max(DUMP_COUNT) * 0.06,
            str(count),
            ha="center",
            va="bottom",
            fontsize=9,
            color="#355E2B",
        )

    for spine in ["top"]:
        ax.spines[spine].set_visible(False)
        ax2.spines[spine].set_visible(False)

    handles = [load_bar[0], recovery_bar[0], dump_line[0]]
    labels = ["Load", "Recovery", "Dump count"]
    fig.legend(
        handles,
        labels,
        ncol=3,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.2,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.095, right=0.90, top=0.98, bottom=0.27)
    fig.savefig(output_dir / "recovery_experiment.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "recovery_experiment.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


if __name__ == "__main__":
    plot(Path("results/figures"))
