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


def plot(output_dir: Path):
    set_paper_style()
    output_dir.mkdir(parents=True, exist_ok=True)

    x = np.arange(len(DATASETS))
    width = 0.18

    fig, ax = plt.subplots(figsize=(7.2, 2.65))

    colors = {
        "zone_algorithm": "#2F5D8C",
        "zone_dump": "#8FB7D8",
        "live": "#D07A2D",
        "lsm": "#5B8E3E",
        "mbf": "#8C6BB1",
    }

    ax.bar(
        x - 1.5 * width,
        ZONEGRAPH_ALGORITHM,
        width,
        label="ZoneGraph update",
        color=colors["zone_algorithm"],
        edgecolor="black",
        linewidth=0.45,
    )
    ax.bar(
        x - 1.5 * width,
        ZONEGRAPH_DUMP,
        width,
        bottom=ZONEGRAPH_ALGORITHM,
        label="ZoneGraph flush",
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
    ax.set_xticklabels(DATASETS, rotation=0, ha="center")
    ax.set_ylim(0.08, 500)
    ax.grid(axis="y", which="major", linestyle="--", linewidth=0.45, alpha=0.45)
    ax.set_axisbelow(True)

    # for spine in ["top", "right"]:
    #     ax.spines[spine].set_visible(False)

    ax.legend(
      ncol=5,
      frameon=True,          # 原来是 False
      edgecolor="black",     # 边框颜色
      fancybox=False,        # 直角边框
      framealpha=1.0,        # 不透明
      loc="upper center",
      bbox_to_anchor=(0.5, -0.12),
      columnspacing=0.8,
      handlelength=1.3,
      handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.075, right=0.995, top=0.995, bottom=0.245)
    fig.savefig(output_dir / "edge_insertion_maintenance.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(
        output_dir / "edge_insertion_maintenance.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.01,
    )
    plt.show()


if __name__ == "__main__":
    plot(Path("results/figures"))
