#!/usr/bin/env python3
"""Plot edge insertion maintenance time for ZoneGraph and baselines."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["DBLP", "YouTube", "Twitch", "Skitter", "Pokec", "LiveJournal"]

ZONEGRAPH_TOTAL = np.array([0.134075, 0.458681, 1.99381, 3.40048, 22.8916, 44.4979])
ZONEGRAPH_DUMP = np.array([0.0, 0.0, 0.687911, 1.58822, 16.4243, 29.4439])
ZONEGRAPH_ALGORITHM = np.array([0.134075, 0.458681, 1.3059, 1.81226, 6.46722, 15.0539])

LIVEGRAPH = np.array([1.49208, 5.03401, 7.08399, 19.9943, 166.939, np.nan])
LSMGRAPH = np.array([1.6, 4.6, 12.0, 18.0, 40.0, 140.0])
MBFGRAPH = np.array([6.222188, 11.223344, 26.20806, 41.046239, 85.538215, 283.933646])


def set_paper_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "hatch.linewidth": 0.5,
        }
    )


def plot(output_dir: Path):
    set_paper_style()
    output_dir.mkdir(parents=True, exist_ok=True)

    x = np.arange(len(DATASETS))
    width = 0.18

    fig, ax = plt.subplots(figsize=(6.8, 2.7))

    colors = {
        "zone_algorithm": "#4C78A8",
        "zone_dump": "#9EC3E6",
        "live": "#F58518",
        "lsm": "#54A24B",
        "mbf": "#B279A2",
    }

    ax.bar(
        x - 1.5 * width,
        ZONEGRAPH_ALGORITHM,
        width,
        label="ZoneGraph algorithm",
        color=colors["zone_algorithm"],
        edgecolor="black",
        linewidth=0.45,
    )
    ax.bar(
        x - 1.5 * width,
        ZONEGRAPH_DUMP,
        width,
        bottom=ZONEGRAPH_ALGORITHM,
        label="ZoneGraph dump",
        color=colors["zone_dump"],
        edgecolor="black",
        linewidth=0.45,
        hatch="//",
    )
    ax.bar(
        x - 0.5 * width,
        LIVEGRAPH,
        width,
        label="LiveGraph",
        color=colors["live"],
        edgecolor="black",
        linewidth=0.45,
    )
    ax.bar(
        x + 0.5 * width,
        LSMGRAPH,
        width,
        label="LSMGraph",
        color=colors["lsm"],
        edgecolor="black",
        linewidth=0.45,
    )
    ax.bar(
        x + 1.5 * width,
        MBFGRAPH,
        width,
        label="MBFGraph",
        color=colors["mbf"],
        edgecolor="black",
        linewidth=0.45,
    )

    # LiveGraph did not finish on LiveJournal in this experiment.
    lj = len(DATASETS) - 1
    ax.text(
        x[lj] - 0.5 * width,
        1.1,
        "N/A",
        ha="center",
        va="bottom",
        rotation=90,
        fontsize=7,
    )

    ax.set_yscale("log")
    ax.set_ylabel("Maintenance time (s)")
    ax.set_xticks(x)
    ax.set_xticklabels(DATASETS, rotation=18, ha="right")
    ax.set_ylim(0.08, 500)
    ax.grid(axis="y", which="major", linestyle="--", linewidth=0.45, alpha=0.45)
    ax.set_axisbelow(True)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.legend(
        ncol=5,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.22),
        columnspacing=0.8,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.tight_layout(pad=0.35)
    fig.savefig(output_dir / "edge_insertion_maintenance.pdf", bbox_inches="tight")
    fig.savefig(output_dir / "edge_insertion_maintenance.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    plot(Path("results/figures"))
