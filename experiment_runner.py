#!/usr/bin/env python3
"""
Closeness Centrality — Experiment Runner
=========================================
Runs selected binaries N times, averages timings, and produces bar charts
with per-run scatter dots overlaid for three experiment types:

  Graph 1 — Time vs Nodes   : n = 100K … 1M (step 100K), edges = n·log₂n, threads=16
  Graph 2 — Time vs Edges   : n = 500K, edges = 1…10 × n·log₂n (step n·log₂n), threads=16
  Graph 3 — Time vs Threads : n = 500K, edges = 5·n·log₂n, threads = 1,2,4,8,16,32
"""

import os, re, csv, shutil, subprocess, math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from itertools import cycle

# ──────────────────────────────────────────────────────────────────────────────
#  Binary metadata:
#    id → (executable_name, short_label, thread_arg_mode, is_parallel)
#  thread_arg_mode: "none" | "optional" | "required"
# ──────────────────────────────────────────────────────────────────────────────
BINARIES = {
    "01": ("01_seq_naive_bfs",        "Seq Naive BFS",      "none",     False),
    "02": ("02_seq_multisource_bfs",  "Seq MS-BFS",         "none",     False),
    "03": ("03_seq_bcc_reduced",      "Seq BCC+R3/R4",      "none",     False),
    "04": ("04_par_simple_mimd",      "Par MIMD",           "required", True),
    "05": ("05_par_mimd_msbfs",       "Par MIMD+MS-BFS",    "required", True),
    "06": ("06_par_level_sync",       "Par Level-Sync",     "optional", True),
    "07": ("07_dynamic_shukla",       "Dyn Shukla",         "required", True),
    "08": ("08_novel_bcc_spmm",       "Novel BCC-SpMM",     "optional", True),
    "09": ("09_novel_vdbcc",          "Novel VD-BCC",       "none",     False),
}

GRAPH_GEN  = "./bi_connected_graph_gen1"   # generates connected_graph.csv
GRAPH_OUT  = "biconnected_graph.csv"
OUTPUT_DIR = "experiment_results"

# Timing patterns accepted from stdout
TIME_PATTERNS = [
    r'Time:\s*([\d.]+)\s*ms',
    r'Parallel time:\s*([\d.]+)\s*ms',
    r'LevelSync time:\s*([\d.]+)\s*ms',
    r'BCC-reduced time:\s*([\d.]+)\s*ms',
    r'BCC:\s*([\d.]+)\s*ms',
]

PALETTE = [
    "#2196F3", "#F44336", "#4CAF50", "#FF9800",
    "#9C27B0", "#00BCD4", "#FF5722", "#607D8B", "#E91E63",
]

# ──────────────────────────────────────────────────────────────────────────────
#  Helpers
# ──────────────────────────────────────────────────────────────────────────────

def nlogn(n: int) -> int:
    return int(n * math.log2(n))


def ensure_dirs():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs("a", exist_ok=True)


def generate_graph(n: int, m: int) -> bool:
    """Call connected_graph_gen and verify output file exists."""
    if not os.path.exists(GRAPH_GEN):
        print(f"  [ERROR] Generator not found: {GRAPH_GEN}")
        print("  Compile with:  g++ -O2 -std=c++17 connected_graph_gen.cpp -o connected_graph_gen")
        return False
    r = subprocess.run([GRAPH_GEN, str(n), str(m)],
                       capture_output=True, text=True, timeout=120)
    return os.path.exists(GRAPH_OUT)


def run_binary(sid: str, graph_file: str, threads: int | None) -> float | None:
    """Run one binary; return elapsed ms parsed from stdout, or None."""
    exe, _, thread_arg_mode, _ = BINARIES[sid]
    path = f"./{exe}"
    if not os.path.exists(path):
        print(f"  [SKIP] {path} not found")
        return None

    cmd = [path, graph_file]
    if thread_arg_mode in ("required", "optional") and threads is not None:
        cmd.append(str(threads))

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        for pat in TIME_PATTERNS:
            m = re.search(pat, result.stdout)
            if m:
                return float(m.group(1))
        print(f"  [WARN] No timing found in stdout of {exe}")
        if result.returncode != 0:
            print(f"         stderr: {result.stderr[:200]}")
        return None
    except subprocess.TimeoutExpired:
        print(f"  [TIMEOUT] {exe}")
        return None
    except Exception as e:
        print(f"  [ERROR] {exe}: {e}")
        return None


def save_csv(filepath: str, header: list, rows: list):
    with open(filepath, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"  CSV saved → {filepath}")


# ──────────────────────────────────────────────────────────────────────────────
#  Plotting
# ──────────────────────────────────────────────────────────────────────────────

def bar_plot_with_dots(
    ax,
    x_labels: list[str],
    data: dict,            # {sid: list[list[float|None]]}  shape [n_x][n_runs]
    selected: list[str],
    title: str,
    xlabel: str,
    ylabel: str = "Time (ms)",
):
    """
    Grouped bar chart.
    Each group = one x-label value.
    Each bar   = one algorithm (average across runs).
    Dots       = individual run times (translucent, random jitter).
    Dashed line= connects runs within a bar (translucent).
    """
    n_groups = len(x_labels)
    n_bars   = len(selected)
    width    = min(0.8 / max(n_bars, 1), 0.18)
    xs       = np.arange(n_groups)
    colors   = PALETTE[:n_bars]

    handles = []

    for bi, (sid, color) in enumerate(zip(selected, colors)):
        label = BINARIES[sid][1]
        offset = (bi - (n_bars - 1) / 2) * width

        means, lows, highs = [], [], []
        all_runs_per_x     = []

        for xi in range(n_groups):
            vals = [v for v in data[sid][xi] if v is not None]
            mean = np.mean(vals) if vals else 0.0
            lo   = mean - np.min(vals)  if vals else 0
            hi   = np.max(vals) - mean  if vals else 0
            means.append(mean)
            lows.append(lo)
            highs.append(hi)
            all_runs_per_x.append(vals)

        x_pos = xs + offset

        # Bars (average)
        bars = ax.bar(
            x_pos, means, width * 0.88,
            color=color, alpha=0.80, zorder=2, label=label,
            edgecolor="white", linewidth=0.5,
        )

        # Error bars
        ax.errorbar(
            x_pos, means, yerr=[lows, highs],
            fmt="none", ecolor=color, elinewidth=1.2,
            capsize=3, alpha=0.6, zorder=3,
        )

        # Per-run dots + connecting dashed line
        rng = np.random.default_rng(42 + bi)
        for xi, (xc, runs) in enumerate(zip(x_pos, all_runs_per_x)):
            if not runs:
                continue
            jitter = rng.uniform(-width * 0.3, width * 0.3, size=len(runs))
            jx = xc + jitter
            ax.scatter(
                jx, runs,
                color=color, alpha=0.35, s=22, zorder=4,
                marker="o", linewidths=0,
            )
            if len(runs) > 1:
                # dashed translucent line connecting dots (ordered by value)
                order = np.argsort(runs)
                ax.plot(
                    jx[order], np.array(runs)[order],
                    color=color, alpha=0.18, linewidth=0.9,
                    linestyle="--", zorder=3,
                )

        handles.append(mpatches.Patch(color=color, label=label, alpha=0.85))

    ax.set_xticks(xs)
    ax.set_xticklabels(x_labels, rotation=30, ha="right", fontsize=9)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.legend(handles=handles, fontsize=8, loc="upper left",
              framealpha=0.9, ncol=max(1, n_bars // 5))
    ax.grid(axis="y", alpha=0.3, linestyle="--", zorder=0)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)


def save_figure(fig, name: str):
    for ext in ("png", "pdf"):
        path = f"{OUTPUT_DIR}/{name}.{ext}"
        fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"  Plot saved → {OUTPUT_DIR}/{name}.png  /  .pdf")


# ──────────────────────────────────────────────────────────────────────────────
#  Experiments
# ──────────────────────────────────────────────────────────────────────────────

def run_experiment(
    label: str,
    selected: list[str],
    n_runs: int,
    configs: list[tuple],    # [(n, m, extra_label), ...]
    threads_per_config: list[int | None],
    x_labels: list[str],
    title: str,
    xlabel: str,
    csv_header: list[str],
    csv_extra_fn,            # fn(config_idx) → extra csv columns
    csv_name: str,
    plot_name: str,
):
    n_configs = len(configs)
    # data[sid][config_idx] = [t1, t2, …]
    data = {sid: [[] for _ in range(n_configs)] for sid in selected}
    csv_rows = []

    for run_i in range(n_runs):
        print(f"\n  ── Run {run_i + 1}/{n_runs} ──")
        for ci, (n, m, _) in enumerate(configs):
            threads = threads_per_config[ci]
            print(f"    n={n:,}  m={m:,}  threads={threads}  → generating … ", end="", flush=True)
            ok = generate_graph(n, m)
            print("ok" if ok else "FAILED")
            if not ok:
                continue

            for sid in selected:
                _, _, thread_arg_mode, _ = BINARIES[sid]
                t = run_binary(sid, GRAPH_OUT, threads if thread_arg_mode != "none" else None)
                data[sid][ci].append(t)
                print(f"      [{sid}] {t:.1f} ms" if t else f"      [{sid}] FAIL")
                csv_rows.append([sid, run_i + 1, *csv_extra_fn(ci, n, m, threads), t])

    save_csv(f"{OUTPUT_DIR}/{csv_name}.csv", csv_header, csv_rows)

    fig, ax = plt.subplots(figsize=(max(10, n_configs * 1.6 + 2), 6))
    bar_plot_with_dots(ax, x_labels, data, selected, title, xlabel)
    plt.tight_layout()
    save_figure(fig, plot_name)
    plt.close(fig)

    # Print summary table
    print(f"\n  {'Binary':<25} " + "  ".join(f"{xl:>9}" for xl in x_labels))
    print("  " + "-" * (25 + 11 * n_configs))
    for sid in selected:
        row_str = f"  {BINARIES[sid][1]:<25}"
        for ci in range(n_configs):
            vals = [v for v in data[sid][ci] if v is not None]
            row_str += f"  {np.mean(vals):>8.1f}" if vals else f"  {'—':>8}"
        print(row_str)

    return data


def experiment1(selected, n_runs):
    print("\n" + "━"*55)
    print("  EXPERIMENT 1 — Time vs Number of Nodes")
    print("  edges = n·log₂n,  threads = 16")
    print("━"*55)

    nodes   = list(range(10_000, 110_000, 10_000))
    THREADS = 16
    configs = [(n, nlogn(n), "") for n in nodes]
    t_list  = [THREADS] * len(configs)
    xlabs   = [f"{n//1000}K" for n in nodes]

    run_experiment(
        label          = "exp1",
        selected       = selected,
        n_runs         = n_runs,
        configs        = configs,
        threads_per_config = t_list,
        x_labels       = xlabs,
        title          = "Time vs Number of Nodes\n(edges = n·log₂n,  threads = 16)",
        xlabel         = "Nodes",
        csv_header     = ["binary", "run", "n", "m", "threads", "time_ms"],
        csv_extra_fn   = lambda ci, n, m, t: [n, m, t],
        csv_name       = "exp1_time_vs_nodes",
        plot_name      = "graph1_time_vs_nodes",
    )


def experiment2(selected, n_runs):
    print("\n" + "━"*55)
    print("  EXPERIMENT 2 — Time vs Number of Edges")
    print("  n = 500 000,  threads = 16")
    print("━"*55)

    N       = 50_000
    THREADS = 16
    base    = nlogn(N)
    mults   = list(range(1, 11))
    configs = [(N, base * k, f"{k}x") for k in mults]
    t_list  = [THREADS] * len(configs)
    xlabs   = [f"{k}·nlog₂n" for k in mults]

    run_experiment(
        label          = "exp2",
        selected       = selected,
        n_runs         = n_runs,
        configs        = configs,
        threads_per_config = t_list,
        x_labels       = xlabs,
        title          = "Time vs Number of Edges\n(n = 50K,  threads = 16)",
        xlabel         = "Edge Count",
        csv_header     = ["binary", "run", "n", "m", "multiplier", "threads", "time_ms"],
        csv_extra_fn   = lambda ci, n, m, t: [n, m, ci + 1, t],
        csv_name       = "exp2_time_vs_edges",
        plot_name      = "graph2_time_vs_edges",
    )


def experiment3(selected, n_runs):
    print("\n" + "━"*55)
    print("  EXPERIMENT 3 — Time vs Number of Threads")
    print("  n = 50 000,  edges = 5·n·log₂n")
    print("━"*55)

    # Only parallel binaries are meaningful here
    par_sel = [sid for sid in selected if BINARIES[sid][2]]
    if not par_sel:
        print("  No parallel binaries in selection — skipping Experiment 3.")
        return

    N       = 50_000
    m       = 5 * nlogn(N)
    threads = [1, 2, 4, 8, 16, 32]
    configs = [(N, m, f"t={t}") for t in threads]
    xlabs   = [str(t) for t in threads]

    run_experiment(
        label          = "exp3",
        selected       = par_sel,
        n_runs         = n_runs,
        configs        = configs,
        threads_per_config = threads,
        x_labels       = xlabs,
        title          = "Time vs Number of Threads\n(n = 50K,  edges = 5·n·log₂n)",
        xlabel         = "Threads",
        csv_header     = ["binary", "run", "n", "m", "threads", "time_ms"],
        csv_extra_fn   = lambda ci, n, m, t: [n, m, t],
        csv_name       = "exp3_time_vs_threads",
        plot_name      = "graph3_time_vs_threads",
    )


# ──────────────────────────────────────────────────────────────────────────────
#  Interactive menu
# ──────────────────────────────────────────────────────────────────────────────

def select_binaries() -> list[str]:
    print("\n┌─────────────────────────────────────────────────────┐")
    print("│               Available Binaries                    │")
    print("├────┬───────────────────────────────┬────────────────┤")
    print("│ ID │ Description                   │ Type           │")
    print("├────┼───────────────────────────────┼────────────────┤")
    for sid, (exe, desc, _, is_parallel) in BINARIES.items():
        kind = "Parallel (OMP)" if is_parallel else "Sequential    "
        exist = "✓" if os.path.exists(f"./{exe}") else "✗ (not built)"
        print(f"│ {sid} │ {desc:<29} │ {kind}  {exist} │")
    print("└────┴───────────────────────────────┴────────────────┘")

    raw = input("\nEnter binary IDs to test (comma-separated, e.g. 01,04,05,08): ").strip()
    selected = [s.strip().zfill(2) for s in raw.split(",")]
    valid    = [s for s in selected if s in BINARIES]

    if not valid:
        print("[ERROR] No valid binary IDs.")
        return []

    missing = [s for s in valid if not os.path.exists(f"./{BINARIES[s][0]}")]
    if missing:
        print(f"[WARN] These binaries are not built yet: {missing}")
        cont = input("       Continue anyway? [y/N]: ").strip().lower()
        if cont != "y":
            return []

    print(f"\nSelected: {', '.join(f'{s} ({BINARIES[s][1]})' for s in valid)}")
    return valid


def main():
    ensure_dirs()

    print("╔══════════════════════════════════════════════════════╗")
    print("║    Closeness Centrality — Experiment Runner          ║")
    print("╚══════════════════════════════════════════════════════╝")

    # Graph generator check
    if not os.path.exists(GRAPH_GEN):
        print(f"\n[ERROR] {GRAPH_GEN} not found.")
        print("  Build with:  make connected_graph_gen")
        return

    selected = select_binaries()
    if not selected:
        return

    try:
        n_runs = int(input("\nHow many times to run each configuration? "))
        assert n_runs >= 1
    except (ValueError, AssertionError):
        print("[ERROR] Enter a positive integer.")
        return

    print("\nWhich experiments to run?")
    print("  1 — Time vs Nodes   (n=100K…1M, edges=n·log₂n, threads=16)")
    print("  2 — Time vs Edges   (n=500K, edges=1…10×n·log₂n, threads=16)")
    print("  3 — Time vs Threads (n=500K, edges=5·n·log₂n, threads=1…32)")
    print("  all — Run all three")
    choice = input("Choice: ").strip().lower()

    if choice == "all":
        run_set = {"1", "2", "3"}
    else:
        run_set = set(x.strip() for x in choice.split(","))

    np.random.seed(0)   # reproducible jitter

    if "1" in run_set:
        experiment1(selected, n_runs)
    if "2" in run_set:
        experiment2(selected, n_runs)
    if "3" in run_set:
        experiment3(selected, n_runs)

    print(f"\n✓ Done!  All results saved to ./{OUTPUT_DIR}/")
    print("  Files: graph1_time_vs_nodes.png/pdf")
    print("         graph2_time_vs_edges.png/pdf")
    print("         graph3_time_vs_threads.png/pdf")
    print("         exp1_time_vs_nodes.csv")
    print("         exp2_time_vs_edges.csv")
    print("         exp3_time_vs_threads.csv")


if __name__ == "__main__":
    main()