#!/usr/bin/env python3
"""Plot graph query time before/after ZoneGraph final dump and baselines."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATASETS = ["DBLP", "YouTube", "Twitch", "Skitter", "Pokec", "LiveJournal"]
QUERIES = ["BFS", "PageRank", "WCC"]

ZG_BEFORE = {
    "BFS": np.array([0.026914, 0.042919, 0.015887, 0.079379, 0.185206, 0.341869]),
    "PageRank": np.array([0.0786, 0.243987, 0.164126, 0.489104, 3.8072, 5.40071]),
    "WCC": np.array([0.028959, 0.085371, 0.069908, 0.085747, 0.420181, 0.994973]),
}

ZG_AFTER = {
    "BFS": np.array([0.053408, 0.166051, 0.403229, 0.600951, 0.250616, 1.22529]),
    "PageRank": np.array([0.151414, 0.577902, 1.63564, 2.00769, 4.77294, 8.76382]),
    "WCC": np.array([0.056992, 0.135559, 0.433172, 0.2982, 0.518758, 1.49308]),
}

LIVEGRAPH = {
    "BFS": np.array([0.0805436, 0.15158, 0.147265, 0.153366, 0.542994, np.nan]),
    "PageRank": np.array([0.33951, 0.898717, 1.20663, 2.37617, 6.33252, np.nan]),
    "WCC": np.array([0.0528316, 0.13502, 0.0952364, 0.637649, 1.40922, np.nan]),
}

LSMGRAPH = {
    "BFS": np.array([0.191396408, 0.30569013, 0.147846, 0.494324537, 0.411547188, 3.466617011]),
    "PageRank": np.array([0.814450924, 3.16108401, 2.534623026, 5.750116254, 8.874712505, 50.06178904]),
    "WCC": np.array([0.183235263, 0.440698265, 0.212878541, 0.531422759, 0.653071546, 3.812760809]),
}

LSMGRAPH_DEDUP = {
    "BFS": np.array([0.18125103, 0.345500576, 0.190377193, 0.614117868, 0.522896218, 4.016680903]),
    "PageRank": np.array([0.898948651, 4.630487858, 4.039286944, 10.81630277, 12.72835988, 58.42036032]),
    "WCC": np.array([0.20555358, 0.489561878, 0.31318945, 0.765834817, 0.830387633, 4.410231376]),
}

MBFGRAPH = {
    "BFS": np.array([13.0911, 13.0034, 6.88078, 19.5575, 9.99074, 22.9653]),
    "PageRank": np.array([9.63443, 12.9893, 4.70072, 17.6659, 7.82024, 18.3401]),
    "WCC": np.array([13.7446, 13.0392, 7.37041, 19.3516, 10.1674, 24.192]),
}

MBFGRAPH_DEDUP = {
    "BFS": np.array([19.2356, 29.4135, 34.3564, 103.51, 119.289, 379.842]),
    "PageRank": np.array([14.2213, 30.7117, 17.6768, 110.097, 118.724, 180.366]),
    "WCC": np.array([21.1603, 31.1209, 37.8207, 126.143, 147.046, 357.887]),
}


def set_paper_style():
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 9,
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

    series = [
        ("ZG-B", ZG_BEFORE, "#2F5D8C", ""),
        ("ZG-A", ZG_AFTER, "#8FB7D8", "//"),
        ("LG", LIVEGRAPH, "#D07A2D", ""),
        ("LSM", LSMGRAPH, "#5B8E3E", ""),
        ("LSM-D", LSMGRAPH_DEDUP, "#A4C97A", "//"),
        ("MBF", MBFGRAPH, "#8C6BB1", ""),
        ("MBF-D", MBFGRAPH_DEDUP, "#C4A3D8", "//"),
    ]

    fig, axes = plt.subplots(3, 6, figsize=(11.4, 6.0), sharey=False)
    x = np.arange(len(series))
    width = 0.72
    legend_handles = None

    for row, query in enumerate(QUERIES):
        for col, dataset in enumerate(DATASETS):
            ax = axes[row, col]
            handles = []
            values = []

            for idx, (label, table, color, hatch) in enumerate(series):
                value = table[query][col]
                values.append(value)
                bar = ax.bar(
                    x[idx],
                    0 if np.isnan(value) else value,
                    width,
                    label=label,
                    color=color,
                    edgecolor="black",
                    linewidth=0.4,
                    hatch=hatch,
                )
                handles.append(bar[0])

            finite_values = [v for v in values if not np.isnan(v) and v > 0]
            ax.set_yscale("log")
            ax.set_ylim(min(finite_values) * 0.55, max(finite_values) * 2.0)
            ax.grid(axis="y", which="major", linestyle="--", linewidth=0.4, alpha=0.45)
            ax.set_axisbelow(True)

            if row == 0:
                ax.set_title(dataset, pad=3)
            if col == 0:
                ax.set_ylabel(f"{query}\ntime (s)")

            ax.set_xticks(x)
            if row == 2:
                ax.set_xticklabels([s[0] for s in series], rotation=35, ha="right")
            else:
                ax.set_xticklabels([])

            if np.isnan(LIVEGRAPH[query][col]):
                ax.text(
                    x[2],
                    min(finite_values) * 0.8,
                    "N/A",
                    ha="center",
                    va="bottom",
                    rotation=90,
                    fontsize=7,
                )

            for spine in ["top", "right"]:
                ax.spines[spine].set_visible(False)

            if legend_handles is None:
                legend_handles = handles

    fig.legend(
        legend_handles,
        [s[0] for s in series],
        ncol=7,
        frameon=False,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.03),
        columnspacing=0.8,
        handlelength=1.2,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.055, right=0.995, top=0.955, bottom=0.16, wspace=0.32, hspace=0.22)
    fig.savefig(output_dir / "graph_query_3x6.pdf", bbox_inches="tight", pad_inches=0.01)
    fig.savefig(output_dir / "graph_query_3x6.png", dpi=600, bbox_inches="tight", pad_inches=0.01)


if __name__ == "__main__":
    plot(Path("results/figures"))
