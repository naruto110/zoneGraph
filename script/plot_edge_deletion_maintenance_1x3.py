#!/usr/bin/env python3
"""Plot edge deletion maintenance time as one-row small multiples."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["Skitter", "Pokec", "LiveJournal"]

ZONEGRAPH = np.array([0.398703, 0.278677, 0.218737])
LIVEGRAPH = np.array([11.1641, 22.0760, np.nan])
LSMGRAPH = np.array([0.946315, 1.98598, 1.63808])
MBFGRAPH = np.array([1.17, 1.21, 1.33])


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

    colors = {
        "zone": "#2F5D8C",
        "live": "#D07A2D",
        "lsm": "#5B8E3E",
        "mbf": "#8C6BB1",
    }

    fig, axes = plt.subplots(1, 3, figsize=(6.4, 2.55), sharey=False)
    x = np.arange(4)
    width = 0.62
    labels = ["ZG", "LG", "LSM", "MBF"]

    legend_handles = None
    for i, (ax, dataset) in enumerate(zip(axes, DATASETS)):
        bars = []
        bars.append(
            ax.bar(
                x[0],
                ZONEGRAPH[i],
                width,
                label="ZoneGraph",
                color=colors["zone"],
                edgecolor="black",
                linewidth=0.45,
            )
        )
        bars.append(
            ax.bar(
                x[1],
                LIVEGRAPH[i] if not np.isnan(LIVEGRAPH[i]) else 0,
                width,
                label="LiveGraph",
                color=colors["live"],
                edgecolor="black",
                linewidth=0.45,
            )
        )
        bars.append(
            ax.bar(
                x[2],
                LSMGRAPH[i],
                width,
                label="LSMGraph",
                color=colors["lsm"],
                edgecolor="black",
                linewidth=0.45,
            )
        )
        bars.append(
            ax.bar(
                x[3],
                MBFGRAPH[i],
                width,
                label="MBFGraph",
                color=colors["mbf"],
                edgecolor="black",
                linewidth=0.45,
            )
        )

        ymax = np.nanmax([ZONEGRAPH[i], LIVEGRAPH[i], LSMGRAPH[i], MBFGRAPH[i]])
        ax.set_ylim(0, ymax * 1.22)

        if np.isnan(LIVEGRAPH[i]):
            ax.text(x[1], ymax * 0.05, "N/A", ha="center", va="bottom", rotation=90, fontsize=9)

        ax.set_title(dataset, pad=3)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=0)
        ax.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.45)
        ax.set_axisbelow(True)

        if i == 0:
            ax.set_ylabel("Deletion time (s)")

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        if legend_handles is None:
            legend_handles = [b[0] for b in bars]

    fig.legend(
        legend_handles,
        ["ZoneGraph", "LiveGraph", "LSMGraph", "MBFGraph"],
        ncol=4,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=1.0,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.085, right=0.995, top=0.86, bottom=0.28, wspace=0.28)
    fig.savefig(output_dir / "edge_deletion_maintenance_1x3.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(
        output_dir / "edge_deletion_maintenance_1x3.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.01,
    )


if __name__ == "__main__":
    plot(Path("results/figures"))
