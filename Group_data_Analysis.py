#!/usr/bin/env python3
import os, re, sys, argparse, csv
import numpy as np
from datetime import datetime

# ----- AXIS LIMITS -----
AX_XMIN, AX_XMAX = -5.0, 24.0
MOV_YMIN, MOV_YMAX = -6.0, 12.0
# -----------------------

INJECTION_MIN = 10.0

# Patterns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- GLOBAL FONT CONFIG ----------------
plt.rcParams["font.family"] = "Calibri"
plt.rcParams["font.size"] = 16
plt.rcParams["axes.titlesize"] = 16
plt.rcParams["axes.labelsize"] = 16
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.rcParams["legend.fontsize"] = 16
plt.rcParams["figure.titlesize"] = 16
# -----------------------------------------------------

ROI_PAT = re.compile(r'^roi[-_](?P<virus>.+)_(?P<analyte>[^_]+)_(?P<side>[^_.]+)\.txt$', re.IGNORECASE)
STATIC_PAT = re.compile(r'(?i)static[_-]map[_-](?P<start>\d+)[_-]to[_-](?P<end>\d+)', re.IGNORECASE)

def group_from_psc_filename(fn):
    base = re.sub(r'(?i)^psc_', '', fn)
    m = ROI_PAT.match(base)
    if not m:
        return None, None
    virus, analyte = m.group("virus"), m.group("analyte")
    group_key = f"{virus} | {analyte}"
    col_name = os.path.splitext(fn)[0]
    return group_key, col_name

def load_psc_series(path):
    try:
        arr = np.loadtxt(path, dtype=float)
    except Exception:
        arr = np.loadtxt(path, dtype=float, delimiter=",")
    if arr.ndim == 1:
        return arr.astype(float)
    return np.nanmean(arr.astype(float), axis=1)

def sanitize(name):
    return re.sub(r'[^A-Za-z0-9._-]+', "_", name.strip())[:200].strip("_")

def movmean_centered(x: np.ndarray, w: int) -> np.ndarray:
    if w <= 1:
        return x.copy()
    n = x.size
    y = np.empty_like(x, dtype=float)
    c = np.empty(n + 1, dtype=float)
    c[0] = 0.0
    c[1:] = np.cumsum(x)
    hl = (w - 1) // 2
    hr = w // 2
    for i in range(n):
        a = max(0, i - hl)
        b = min(n, i + hr + 1)
        y[i] = (c[b] - c[a]) / (b - a)
    return y

def minutes_axis(n_points: int, tr_sec: float) -> np.ndarray:
    return (np.arange(n_points) * tr_sec) / 60.0

def apply_x_limits(ax):
    ax.set_xlim(AX_XMIN, AX_XMAX)

def apply_mov_y_limits(ax):
    ax.set_ylim(MOV_YMIN, MOV_YMAX)

def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(2.0)
    ax.spines["bottom"].set_linewidth(2.0)
    ax.tick_params(axis="both")

def clip_xy_to_xlim(x, *ys, xmin=AX_XMIN, xmax=AX_XMAX):
    """
    Clip x and corresponding y-arrays to [xmin, xmax].
    Returns (x_clipped, y1_clipped, y2_clipped, ...).
    Any y that is None will be passed through as None.
    """
    x = np.asarray(x)
    mask = (x >= xmin) & (x <= xmax)
    x_c = x[mask]
    ys_c = []
    for y in ys:
        if y is None:
            ys_c.append(None)
        else:
            y_arr = np.asarray(y)
            ys_c.append(y_arr[mask])
    return (x_c, *ys_c)

# ---------------- PSC generation helpers ----------------

def find_baseline_window_in_dir(dirpath, verbose=False):
    candidates = [f for f in os.listdir(dirpath) if re.search(r"(?i)static[_-]map", f)]
    for fn in candidates:
        m = STATIC_PAT.search(fn)
        if m:
            a, b = int(m.group("start")), int(m.group("end"))
            return (min(a, b), max(a, b))
    return None

def compute_psc_from_series(series, base_start, base_end):
    n = len(series)
    a = max(0, min(base_start, n - 1))
    b = max(0, min(base_end, n - 1))
    if b < a:
        a, b = b, a
    baseline = float(np.nanmean(series[a:b+1]))
    if baseline == 0 or np.isnan(baseline):
        return None
    return ((series - baseline) / baseline) * 100.0

def generate_psc_files(root, verbose=False):
    made, skipped = 0, 0
    for dirpath, _, files in os.walk(root):
        roi_files = [f for f in files if re.match(r"(?i)^roi[-_].*\.txt$", f)]
        if not roi_files:
            continue
        win = find_baseline_window_in_dir(dirpath, verbose)
        if not win:
            skipped += len(roi_files)
            continue
        base_start, base_end = win
        for fn in roi_files:
            inp = os.path.join(dirpath, fn)
            outp = os.path.join(dirpath, f"psc_{fn}")
            if os.path.exists(outp):
                continue
            try:
                arr = np.loadtxt(inp, dtype=float)
            except Exception:
                arr = np.loadtxt(inp, dtype=float, delimiter=",")
            if arr.ndim == 1:
                arr = arr.reshape(-1, 1)
            psc_cols, invalid = [], False
            for j in range(arr.shape[1]):
                psc = compute_psc_from_series(arr[:, j], base_start, base_end)
                if psc is None:
                    invalid = True
                    break
                psc_cols.append(psc)
            if invalid:
                skipped += 1
                continue
            psc_mat = np.column_stack(psc_cols)
            np.savetxt(outp, psc_mat, fmt="%.6f", delimiter="\t")
            made += 1
    return made, skipped

def needs_autogen_psc(root, verbose=False):
    roi_nii = re.compile(r"(?i)^roi[-_].*\.nii(\.gz)?$")
    roi_txt = re.compile(r"(?i)^roi[-_].*\.txt$")
    psc_txt = re.compile(r"(?i)^psc_.*\.txt$")
    for dirpath, _, files in os.walk(root):
        if any(roi_nii.match(f) for f in files):
            if any(roi_txt.match(f) for f in files) and not any(psc_txt.match(f) for f in files):
                return True
    return False

# ---------------- Main analysis ----------------
def main():
    ap = argparse.ArgumentParser(
        description="Generate PSC, group averages, movmean, SEM, and unified PSC plots."
    )
    ap.add_argument("root")
    ap.add_argument("--tr", type=float, required=True)
    ap.add_argument("--outdir")
    ap.add_argument("--min-series", type=int, default=2)
    ap.add_argument("--strict-length", action="store_true")
    ap.add_argument("--movmean", type=int, default=0)
    ap.add_argument("--generate-psc", action="store_true")
    ap.add_argument("--verbose", "-v", action="store_true")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    if not os.path.isdir(root):
        sys.exit(f"Error: '{root}' is not a directory.")

    if args.tr <= 0:
        sys.exit("Error: --tr must be > 0.")

    if args.outdir:
        outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(root, args.outdir)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = os.path.join(root, f"{ts}_Group_Data")
    os.makedirs(outdir, exist_ok=True)

    if args.generate_psc:
        generate_psc_files(root, verbose=args.verbose)
    else:
        if needs_autogen_psc(root, verbose=args.verbose):
            generate_psc_files(root, verbose=args.verbose)

    # ----- Collect PSC filenames -----
    groups = {}
    psc_paths = []

    for dirpath, _, files in os.walk(root):
        for fn in files:
            if not re.match(r"(?i)^psc_.*\.txt$", fn):
                continue
            full = os.path.join(dirpath, fn)
            gkey, colname = group_from_psc_filename(fn)
            if not gkey:
                continue
            groups.setdefault(gkey, []).append((full, colname))
            psc_paths.append(full)

    if not psc_paths:
        print("[ERROR] No PSC files found.")
        return

    # Color cycle
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0', 'C1', 'C2', 'C3'])
    sorted_keys = sorted(groups.keys(), key=lambda s: s.lower())
    group_colors = {k: color_cycle[i % len(color_cycle)] for i, k in enumerate(sorted_keys)}

    # ----------------------------------------------------------------------
    # Individual PSC plots
    # ----------------------------------------------------------------------
    for full in psc_paths:
        fn = os.path.basename(full)
        gkey, _ = group_from_psc_filename(fn)
        if not gkey:
            continue

        color = group_colors.get(gkey, "black")

        try:
            s_file = load_psc_series(full)
            raw_t = minutes_axis(len(s_file), args.tr)
            tmin = raw_t - INJECTION_MIN

            if args.movmean and args.movmean > 1:
                s_mov = movmean_centered(s_file, args.movmean)
            else:
                s_mov = None

            # Clip to [AX_XMIN, AX_XMAX] before plotting
            t_plot, s_file_plot, s_mov_plot = clip_xy_to_xlim(tmin, s_file, s_mov)
            if t_plot.size == 0:
                continue

            svg_path = os.path.splitext(full)[0] + ".svg"

            fig, ax = plt.subplots(figsize=(7, 5), dpi=300)

            ax.plot(t_plot, s_file_plot, linewidth=2, color=color, alpha=0.7, label="Raw PSC")

            if s_mov_plot is not None:
                ax.plot(t_plot, s_mov_plot, linewidth=2.5, color=color,
                        label=f"Movmean (w={args.movmean})")

            xticks = np.arange(int(AX_XMIN), int(AX_XMAX) + 1, 5)
            ax.set_xticks(xticks)

            ax.axvline(0, color='black', linestyle='--', linewidth=1.5)
            ax.axvspan(0.0, 10.0, alpha=0.15, color='0.2')

            ax.set_xlabel("Time (min)")
            ax.set_ylabel("Percent Signal Change (%)")
            ax.set_title(f"Individual PSC — {fn}")

            apply_x_limits(ax)
            apply_mov_y_limits(ax)
            style_axes(ax)
            ax.legend()

            fig.tight_layout()
            fig.savefig(svg_path, format="svg")
            plt.close(fig)

        except Exception as e:
            print(f"[WARN] Could not render PSC for {full}: {e}")

    # ----------------------------------------------------------------------
    # GROUP-LEVEL PROCESSING
    # ----------------------------------------------------------------------
    combined_series = []
    scatter_data = []
    group_stats = []

    for gkey in sorted_keys:
        items = groups[gkey]

        series_list = []
        names = []
        lengths = []

        for path, colname in sorted(items, key=lambda x: x[1].lower()):
            try:
                s = load_psc_series(path)
            except Exception:
                continue
            series_list.append(s)
            names.append(colname)
            lengths.append(len(s))

        if not series_list:
            continue

        # Align length
        if len(set(lengths)) != 1:
            min_len = min(lengths)
            series_list = [s[:min_len] for s in series_list]
        else:
            min_len = lengths[0]

        mat = np.column_stack(series_list)
        M = mat.shape[1]
        safe_group = sanitize(gkey)

        # Save raw matrix CSV
        csv_path = os.path.join(outdir, f"{safe_group}.csv")
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["index"] + names)
            for i in range(min_len):
                w.writerow([i] + [f"{val:.6f}" for val in mat[i, :]])

        # Movmean columnwise
        if args.movmean and args.movmean > 1:
            sm_mat = np.column_stack([
                movmean_centered(mat[:, j], args.movmean) for j in range(M)
            ])
        else:
            sm_mat = mat.copy()

        movmean_csv = os.path.join(outdir, f"mov_mean_{safe_group}.csv")
        with open(movmean_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["index"] + names)
            for i in range(min_len):
                w.writerow([i] + [f"{val:.6f}" for val in sm_mat[i, :]])

        # Average + SEM
        avg = np.nanmean(mat, axis=1)
        sem = np.nanstd(mat, axis=1, ddof=1) / np.sqrt(M) if M >= 2 else None

        avg_path = os.path.join(outdir, f"average_{safe_group}.txt")
        np.savetxt(avg_path, avg, fmt="%.6f", delimiter="\t")

        t_raw = minutes_axis(min_len, args.tr)
        tmin = t_raw - INJECTION_MIN

        # ----------------------------------------------------------------------
        # GROUP AVERAGE PLOT (clipped)
        # ----------------------------------------------------------------------
        t_plot, avg_plot, sem_plot = clip_xy_to_xlim(tmin, avg, sem)
        if t_plot.size == 0:
            mov_txt = ""
            mov_svg = ""
            group_stats.append((gkey, M, min_len, csv_path, avg_path, "", "", ""))
            continue

        avg_svg = os.path.join(outdir, f"Average_time_Course_{safe_group}.svg")
        fig, ax = plt.subplots(figsize=(5, 8), dpi=1200)

        ax.plot(t_plot, avg_plot, linewidth=2.2, color=group_colors[gkey], label="Average")

        if sem_plot is not None:
            ax.fill_between(t_plot, avg_plot - sem_plot, avg_plot + sem_plot,
                            alpha=0.20, color=group_colors[gkey])

        xticks = np.arange(int(AX_XMIN), int(AX_XMAX) + 1, 5)
        ax.set_xticks(xticks)

        ax.axvline(0, color="black", linestyle="--", linewidth=1.5)
        ax.axvspan(0.0, 10.0, alpha=0.15, color="0.2")

        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Percent Signal Change (%)")
        ax.set_title(f"Average Time Course — {gkey}")
        ax.legend()

        apply_x_limits(ax)
        apply_mov_y_limits(ax)
        style_axes(ax)

        fig.tight_layout()
        fig.savefig(avg_svg, format="svg")
        plt.close(fig)

        # ----------------------------------------------------------------------
        # GROUP MOVMEAN PLOT (clipped)
        # ----------------------------------------------------------------------
        mov_avg = mov_sem = None
        mov_txt = ""
        mov_svg = ""

        if args.movmean and args.movmean > 1:
            mov_avg = movmean_centered(avg, args.movmean)
            mov_sem = movmean_centered(sem, args.movmean) if sem is not None else None

            t_mov_plot, mov_avg_plot, mov_sem_plot = clip_xy_to_xlim(tmin, mov_avg, mov_sem)
            if t_mov_plot.size > 0:
                mov_svg = os.path.join(outdir, f"Movmean_Average_{safe_group}.svg")
                mov_txt = os.path.join(outdir, f"Movmean_Average_{safe_group}.txt")

                # Save full-length mov_avg to txt (not clipped)
                np.savetxt(mov_txt, mov_avg, fmt="%.6f", delimiter="\t")

                fig, ax = plt.subplots(figsize=(5, 8), dpi=1200)

                ax.plot(t_mov_plot, mov_avg_plot, linewidth=2.2, color=group_colors[gkey],
                        label=f"Movmean (w={args.movmean})")

                if mov_sem_plot is not None:
                    ax.fill_between(t_mov_plot,
                                    mov_avg_plot - mov_sem_plot,
                                    mov_avg_plot + mov_sem_plot,
                                    alpha=0.20, color=group_colors[gkey])

                xticks = np.arange(int(AX_XMIN), int(AX_XMAX) + 1, 5)
                ax.set_xticks(xticks)

                ax.axvline(0, color="black", linestyle="--", linewidth=1.5)
                ax.axvspan(0.0, 10.0, alpha=0.15, color="0.2")

                ax.set_xlabel("Time (min)")
                ax.set_ylabel("% Signal Change")
                ax.set_title(f"Movmean — {gkey}")
                ax.legend()

                apply_x_limits(ax)
                apply_mov_y_limits(ax)
                style_axes(ax)

                fig.tight_layout()
                fig.savefig(mov_svg, format="svg")
                plt.close(fig)

                # store clipped series for combined plot
                combined_series.append((gkey, t_mov_plot, mov_avg_plot, mov_sem_plot))
            else:
                mov_avg = mov_sem = None  # nothing in visible range

        # ----- SCATTER STATISTICS -----
        mask_20_30 = (tmin >= 20.0) & (tmin <= 30.0)
        maxima_per_file = [float(np.nanmax(sm_mat[mask_20_30, j]))
                           for j in range(sm_mat.shape[1])]

        if mov_avg is not None:
            group_summary = float(np.nanmax(mov_avg[mask_20_30]))
        else:
            group_summary = float(np.nanmax(avg[mask_20_30]))

        scatter_data.append((gkey, group_colors[gkey], maxima_per_file, group_summary))

        group_stats.append(
            (
                gkey,
                M,
                min_len,
                csv_path,
                avg_path,
                avg_svg,
                mov_txt if args.movmean > 1 else "",
                mov_svg if args.movmean > 1 else "",
            )
        )

    # ----------------------------------------------------------------------
    # COMBINED MOVMEAN PLOT
    # ----------------------------------------------------------------------
    if combined_series:
        combined_svg = os.path.join(outdir, "Combined_Movmean_Averages.svg")
        fig, ax = plt.subplots(figsize=(5, 8), dpi=1200)

        for gkey, t_plot, mov_avg_plot, mov_sem_plot in sorted(combined_series):
            col = group_colors[gkey]
            ax.plot(t_plot, mov_avg_plot, linewidth=2, label=gkey, color=col)
            if mov_sem_plot is not None:
                ax.fill_between(t_plot,
                                mov_avg_plot - mov_sem_plot,
                                mov_avg_plot + mov_sem_plot,
                                alpha=0.12, color=col)

        xticks = np.arange(int(AX_XMIN), int(AX_XMAX) + 1, 5)
        ax.set_xticks(xticks)

        ax.axvline(0, color="black", linestyle="--", linewidth=1.5)
        ax.axvspan(0.0, 10.0, alpha=0.15, color="0.2")

        ax.set_xlabel("Time (min)")
        ax.set_ylabel("% Signal Change")
        ax.set_title("Combined Movmean Averages")
        ax.legend()

        apply_x_limits(ax)
        apply_mov_y_limits(ax)
        style_axes(ax)

        fig.tight_layout()
        fig.savefig(combined_svg, format="svg")
        plt.close(fig)

    # ----------------------------------------------------------------------
    # SCATTER PLOT (20–30 min)
    # ----------------------------------------------------------------------
    if scatter_data:
        scatter_svg = os.path.join(outdir, "Scatter_Max_20_30min_by_Group.svg")
        fig, ax = plt.subplots(figsize=(5, 8), dpi=1200)

        x_positions = np.arange(1, len(scatter_data) + 1)

        for i, (gkey, color, maxima, group_summary) in enumerate(scatter_data, start=1):
            for y in maxima:
                ax.scatter(i, y, s=60, color=color)
            ax.hlines(y=group_summary, xmin=i - 0.25, xmax=i + 0.25,
                      colors=color, linewidth=2)

        ax.set_xticks(x_positions)
        ax.set_xticklabels([item[0] for item in scatter_data],
                           rotation=25, ha="right")

        ax.set_ylabel("% Signal Change")
        ax.set_title("Max (20–30 min): Per-file vs Group Peak")

        apply_mov_y_limits(ax)
        style_axes(ax)
        ax.margins(x=0.02)

        fig.tight_layout()
        fig.savefig(scatter_svg, format="svg")
        plt.close(fig)

    print("\n=== GROUP BUILD COMPLETE ===")
    print("Output folder:", outdir)

# ----------------------------------------------------------------------
if __name__ == "__main__":
    main()
