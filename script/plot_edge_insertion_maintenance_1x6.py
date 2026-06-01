#!/usr/bin/env python3
"""Plot edge insertion maintenance time as one-row small multiples."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["DBLP", "YouTube", "Twitch", "Skitter", "Pokec", "LiveJournal"]

ZONEGRAPH_TOTAL = np.array([0.134075, 0.458681, 1.99381, 3.40048, 22.8916, 44.4979])
ZONEGRAPH_FLUSH = np.array([0.0, 0.0, 0.687911, 1.58822, 16.4243, 29.4439])
ZONEGRAPH_UPDATE = np.array([0.134075, 0.458681, 1.3059, 1.81226, 6.46722, 15.0539])

LIVEGRAPH = np.array([1.49208, 5.03401, 7.08399, 19.9943, 166.939, np.nan])
LSMGRAPH = np.array([1.6, 4.6, 12.0, 18.0, 40.0, 140.0])
MBFGRAPH = np.array([6.222188, 11.223344, 26.20806, 41.046239, 85.538215, 283.933646])


def set_paper_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 11,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 10,
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

    colors = {
        "zone_update": "#2F5D8C",
        "zone_flush": "#8FB7D8",
        "live": "#D07A2D",
        "lsm": "#5B8E3E",
        "mbf": "#8C6BB1",
    }

    fig, axes = plt.subplots(1, 6, figsize=(10.2, 2.65), sharey=False)
    x = np.arange(4)
    width = 0.62

    legend_handles = None
    for i, (ax, dataset) in enumerate(zip(axes, DATASETS)):
        bars = []
        bars.append(
            ax.bar(
                x[0],
                ZONEGRAPH_UPDATE[i],
                width,
                label="ZoneGraph update",
                color=colors["zone_update"],
                edgecolor="black",
                linewidth=0.45,
            )
        )
        bars.append(
            ax.bar(
                x[0],
                ZONEGRAPH_FLUSH[i],
                width,
                bottom=ZONEGRAPH_UPDATE[i],
                label="ZoneGraph flush",
                color=colors["zone_flush"],
                edgecolor="black",
                linewidth=0.45,
                hatch="//",
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

        if np.isnan(LIVEGRAPH[i]):
            ax.text(x[1], max(LSMGRAPH[i], MBFGRAPH[i]) * 0.04, "N/A", ha="center", va="bottom", rotation=90, fontsize=8)

        ymax = np.nanmax([ZONEGRAPH_TOTAL[i], LIVEGRAPH[i], LSMGRAPH[i], MBFGRAPH[i]])
        ax.set_ylim(0, ymax * 1.22)
        ax.set_title(dataset, pad=3)
        ax.set_xticks(x)
        ax.set_xticklabels(["ZG", "LG", "LSM", "MBF"], rotation=0)
        ax.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.45)
        ax.set_axisbelow(True)

        if i == 0:
            ax.set_ylabel("Maintenance time (s)")
        else:
            ax.set_ylabel("")

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        if legend_handles is None:
            legend_handles = [b[0] for b in bars]

    fig.legend(
        legend_handles,
        ["ZoneGraph update", "ZoneGraph flush", "LiveGraph", "LSMGraph", "MBFGraph"],
        ncol=5,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.02),
        columnspacing=0.9,
        handlelength=1.3,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.055, right=0.995, top=0.88, bottom=0.28, wspace=0.34)
    fig.savefig(output_dir / "edge_insertion_maintenance_1x6.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(
        output_dir / "edge_insertion_maintenance_1x6.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.01,
    )


if __name__ == "__main__":
    plot(Path("results/figures"))
