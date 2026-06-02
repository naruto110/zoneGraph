#!/usr/bin/env python3
"""Plot graph maintenance thread-count ablation."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


THREADS = np.array([1, 4, 8, 64])
DATA = {
    "Skitter": np.array([6.47689, 5.22297, 5.29004, 5.35957]),
    "Pokec": np.array([58.1474, 35.285, 31.4887, 28.4993]),
    "LiveJournal": np.array([128.227, 80.3506, 72.9466, 67.2383]),
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
        }
    )


def plot(output_dir: Path):
    set_paper_style()
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(5.4, 2.9))

    colors = {
        "Skitter": "#2F5D8C",
        "Pokec": "#D07A2D",
        "LiveJournal": "#5B8E3E",
    }
    markers = {
        "Skitter": "o",
        "Pokec": "s",
        "LiveJournal": "^",
    }

    x = np.arange(len(THREADS))
    for dataset, values in DATA.items():
        ax.plot(
            x,
            values,
            marker=markers[dataset],
            markersize=4.8,
            linewidth=1.6,
            label=dataset,
            color=colors[dataset],
        )

    ax.set_xticks(x)
    ax.set_xticklabels([str(t) for t in THREADS])
    ax.set_xlabel("Threads")
    ax.set_ylabel("Load time (s)")
    ax.set_ylim(0, max(v.max() for v in DATA.values()) * 1.12)
    ax.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.45)
    ax.set_axisbelow(True)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.legend(
        ncol=3,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        columnspacing=1.0,
        handlelength=1.5,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.105, right=0.995, top=0.98, bottom=0.28)
    fig.savefig(output_dir / "ablation_threads.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "ablation_threads.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


if __name__ == "__main__":
    plot(Path("results/figures"))
