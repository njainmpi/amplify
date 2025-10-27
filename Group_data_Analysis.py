#!/usr/bin/env python3
import os, re, sys, argparse, csv
import numpy as np
from datetime import datetime

# ----- AXIS LIMITS -----
AX_XMIN, AX_XMAX = 0.0, 30.0     # minutes for ALL line plots
MOV_YMIN, MOV_YMAX = -6.0, 15.0  # y-limits for Movmean_* and Combined plot (and scatter)
# -----------------------

# Patterns
ROI_PAT = re.compile(r'^roi[-_](?P<virus>.+)_(?P<analyte>[^_]+)_(?P<side>[^_.]+)\.txt$', re.IGNORECASE)
STATIC_PAT = re.compile(r'(?i)static[_-]map[_-](?P<start>\d+)[_-]to[_-](?P<end>\d+)', re.IGNORECASE)

def group_from_psc_filename(fn):
    """Returns ('<virus> | <analyte>', column_name) or (None, None) if not matching."""
    base = re.sub(r'(?i)^psc_', '', fn)  # strip leading psc_
    m = ROI_PAT.match(base)
    if not m:
        return None, None
    virus, analyte = m.group("virus"), m.group("analyte")
    group_key = f"{virus} | {analyte}"
    col_name = os.path.splitext(fn)[0]
    return group_key, col_name

def load_psc_series(path):
    """Load PSC as 1D array. If multiple columns, average across columns."""
    try:
        arr = np.loadtxt(path, dtype=float)
    except Exception:
        arr = np.loadtxt(path, dtype=float, delimiter=',')
    if arr.ndim == 1:
        return arr.astype(float)
    return np.nanmean(arr.astype(float), axis=1)

def sanitize(name):
    return re.sub(r'[^A-Za-z0-9._-]+', '_', name.strip())[:200].strip('_')

def movmean_centered(x: np.ndarray, w: int) -> np.ndarray:
    """Centered moving mean with partial windows at edges (no zero padding)."""
    if w <= 1:
        return x.copy()
    n = x.size
    y = np.empty_like(x, dtype=float)
    c = np.empty(n + 1, dtype=float); c[0] = 0.0; c[1:] = np.cumsum(x, dtype=float)
    hl = (w - 1) // 2; hr = w // 2
    for i in range(n):
        a = max(0, i - hl); b = min(n, i + hr + 1)  # exclusive
        y[i] = (c[b] - c[a]) / (b - a)
    return y

def minutes_axis(n_points: int, tr_sec: float) -> np.ndarray:
    """Return an array of time (in minutes) for 0..n_points-1 given TR in seconds."""
    return (np.arange(n_points, dtype=float) * tr_sec) / 60.0

def apply_x_limits(ax):
    ax.set_xlim(AX_XMIN, AX_XMAX)

def apply_mov_y_limits(ax):
    ax.set_ylim(MOV_YMIN, MOV_YMAX)

def style_axes(ax):
    """Apply cosmetics: spines, linewidths, tick/font sizes."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(3.0)
    ax.spines['bottom'].set_linewidth(3.0)
    ax.tick_params(axis='both', labelsize=40)

# ---------------- PSC generation helpers ----------------

def find_baseline_window_in_dir(dirpath, verbose=False):
    """
    Search dir for a file like:
      sm_coreg_func_Static_Map_200_to_400_and_1500_to_1700.nii.gz
    and return (base_start, base_end) from the FIRST '<start>_to_<end>' after 'Static_Map'.
    """
    candidates = [fn for fn in os.listdir(dirpath) if re.search(r'(?i)static[_-]map', fn)]
    for fn in candidates:
        m = STATIC_PAT.search(fn)
        if m:
            a, b = int(m.group("start")), int(m.group("end"))
            if verbose:
                print(f"[PSC] Static_Map match in {dirpath}: {fn} -> baseline {a}..{b}")
            return (a, b) if a <= b else (b, a)
    return None

def compute_psc_from_series(series, base_start, base_end):
    """
    PSC = ((x - baseline)/baseline)*100, baseline is mean over rows [base_start..base_end] inclusive.
    Clips indices to series length. Returns None if baseline is 0 or NaN.
    """
    n = len(series)
    a = max(0, min(base_start, n-1))
    b = max(0, min(base_end,   n-1))
    if b < a:
        a, b = b, a
    baseline = float(np.nanmean(series[a:b+1]))
    if baseline == 0 or np.isnan(baseline):
        return None
    return ((series - baseline) / baseline) * 100.0

def generate_psc_files(root, verbose=False):
    """
    For each roi_*.txt, compute PSC using Static_Map_* in same folder; write psc_roi_*.txt if missing.
    Supports single- or multi-column roi files. Writes tab-delimited.
    """
    made, skipped = 0, 0
    for dirpath, _, files in os.walk(root):
        roi_files = [f for f in files if re.match(r'(?i)^roi[-_].*\.txt$', f)]
        if not roi_files:
            continue
        win = find_baseline_window_in_dir(dirpath, verbose=verbose)
        if not win:
            if verbose:
                print(f"[PSC] No Static_Map_* found in: {dirpath} -> skip {len(roi_files)} roi files")
            skipped += len(roi_files)
            continue
        base_start, base_end = win
        for fn in roi_files:
            in_path = os.path.join(dirpath, fn)
            out_name = f"psc_{fn}"
            out_path = os.path.join(dirpath, out_name)
            if os.path.exists(out_path):
                if verbose:
                    print(f"[PSC] Exists -> {out_path}")
                continue
            try:
                arr = np.loadtxt(in_path, dtype=float)
            except Exception:
                arr = np.loadtxt(in_path, dtype=float, delimiter=',')
            if arr.ndim == 1:
                arr = arr.reshape(-1, 1)
            psc_cols = []
            invalid = False
            for j in range(arr.shape[1]):
                psc = compute_psc_from_series(arr[:, j].astype(float), base_start, base_end)
                if psc is None:
                    invalid = True
                    if verbose:
                        print(f"[PSC] Invalid baseline (zero/NaN) in {in_path} col {j}; skipping file.")
                    break
                psc_cols.append(psc)
            if invalid:
                skipped += 1
                continue
            psc_mat = np.column_stack(psc_cols)
            np.savetxt(out_path, psc_mat, fmt="%.6f", delimiter="\t")
            made += 1
            if verbose:
                print(f"[PSC] Wrote -> {out_path}  (baseline {base_start}..{base_end})")
    return made, skipped

def needs_autogen_psc(root, verbose=False):
    """
    Return True if we find any folder that has roi*.nii.gz AND roi*.txt
    but lacks psc_roi*.txt — i.e., we should auto-generate PSC.
    """
    roi_nii_pat = re.compile(r'(?i)^roi[-_].*\.nii(\.gz)?$')
    roi_txt_pat = re.compile(r'(?i)^roi[-_].*\.txt$')
    psc_txt_pat = re.compile(r'(?i)^psc_.*\.txt$')

    for dirpath, _, files in os.walk(root):
        has_roi_nii = any(roi_nii_pat.match(f) for f in files)
        if not has_roi_nii:
            continue
        has_roi_txt = any(roi_txt_pat.match(f) for f in files)
        has_psc_txt = any(psc_txt_pat.match(f) for f in files)
        if has_roi_txt and not has_psc_txt:
            if verbose:
                print(f"[PSC-AUTO] {dirpath} has roi*.nii.gz and roi*.txt but no psc_roi*.txt")
            return True
    return False

# ---------------- Main analysis ----------------

def main():
    ap = argparse.ArgumentParser(
        description="(Optionally) generate PSC from roi*.txt using Static_Map_<start>_to_<end>, then build group matrices, averages, SEM, movmean, and SVG plots from psc_roi*.txt. Also saves per-file psc_roi*.svg next to each text file."
    )
    ap.add_argument("root", help="Root folder to search")
    ap.add_argument("--tr", type=float, required=True,
                    help="Repetition time (TR) in seconds; x-axis is shown in minutes")
    ap.add_argument("--outdir", help="Custom output subfolder (default: <YYYYMMDD_HHMMSS>_Group_Data under ROOT)")
    ap.add_argument("--min-series", type=int, default=2, help="Min files required to build a group (default: 2)")
    ap.add_argument("--strict-length", action="store_true",
                    help="Error if series lengths differ (default: truncate to shortest)")
    ap.add_argument("--movmean", type=int, default=0,
                    help="Centered moving-mean window length (e.g., 5, 11). If 0, movmean is skipped.")
    ap.add_argument("--generate-psc", action="store_true",
                    help="Create psc_roi*.txt from roi*.txt using Static_Map_<start>_to_<end> in the same folder, then proceed.")
    ap.add_argument("--verbose", "-v", action="store_true", help="Print progress")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    if not os.path.isdir(root):
        sys.exit(f"Error: '{root}' is not a directory.")
    if args.tr <= 0:
        sys.exit("Error: --tr must be > 0 seconds.")

    # Matplotlib setup
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams['xtick.labelsize'] = 40
    plt.rcParams['ytick.labelsize'] = 40

    # Output folder (under ROOT)
    if args.outdir:
        outdir = args.outdir if os.path.isabs(args.outdir) else os.path.join(root, args.outdir)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = os.path.join(root, f"{ts}_Group_Data")
    os.makedirs(outdir, exist_ok=True)

    if args.verbose:
        print(f"[INFO] ROOT = {root}")
        print(f"[INFO] Output folder = {outdir}")
        print(f"[INFO] TR = {args.tr} sec  -> x-axis in minutes")
        print(f"[INFO] X-limits on ALL plots: [{AX_XMIN}, {AX_XMAX}] minutes")
        print(f"[INFO] Y-limits on Movmean & Combined (and scatter): [{MOV_YMIN}, {MOV_YMAX}] %PSC")
        if args.movmean and args.movmean > 1:
            print(f"[INFO] Movmean window = {args.movmean}")
        else:
            print(f"[INFO] Movmean = disabled (mov_mean CSV will match raw)")

    # --- Optional step 1: generate PSC files in-place ---
    if args.generate_psc:
        made, skipped = generate_psc_files(root, verbose=args.verbose)
        if args.verbose:
            print(f"[PSC] Created: {made}, skipped (no Static_Map or invalid baseline): {skipped}")
    else:
        # Auto-generate when a folder has roi*.nii.gz and roi*.txt but no psc_roi*.txt
        if needs_autogen_psc(root, verbose=args.verbose):
            if args.verbose:
                print("[PSC-AUTO] Auto-generating psc_roi*.txt where roi*.txt present and psc missing...")
            made, skipped = generate_psc_files(root, verbose=args.verbose)
            if args.verbose:
                print(f"[PSC-AUTO] Created: {made}, skipped: {skipped}")
    # ----------------------------------------------------

    # Collect psc files
    groups = {}  # key -> list[(full_path, col_name)]
    psc_count = 0
    for dirpath, _, files in os.walk(root):
        for fn in files:
            if not re.match(r'(?i)^psc_.*\.txt$', fn):
                continue
            gkey, colname = group_from_psc_filename(fn)
            if not gkey:
                if args.verbose:
                    print(f"[WARN] Skipping (name not matching roi pattern): {os.path.join(dirpath, fn)}")
                continue
            full = os.path.join(dirpath, fn)
            groups.setdefault(gkey, []).append((full, colname))
            psc_count += 1
            if args.verbose:
                print(f"[FOUND] {full}  -> group = {gkey}")

            # --------- per-file psc SVG (saved next to the .txt) ----------
            try:
                s_file = load_psc_series(full)  # average across columns if multi-col
                t_file = minutes_axis(len(s_file), args.tr)
                svg_path = os.path.splitext(full)[0] + ".svg"  # same folder, psc_roi*.svg
                fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150)
                ax.plot(t_file, s_file, linewidth=3)
                ax.axvspan(10.0, 20.0, alpha=0.15, color='0.2')
                ax.set_xlabel("Time (in minutes)", fontsize=40)
                ax.set_ylabel("Percent Signal Change (%)", fontsize=40)
                ax.set_title(os.path.basename(svg_path))
                apply_x_limits(ax)          # X-limits on all plots
                style_axes(ax)
                fig.tight_layout()
                fig.savefig(svg_path, format="svg")
                plt.close(fig)
                if args.verbose:
                    print(f"[OK] Per-file PSC SVG -> {svg_path}")
            except Exception as e:
                print(f"[WARN] Could not render PSC SVG for {full}: {e}")
            # -------------------------------------------------------------------

    if psc_count == 0:
        print(f"[ERROR] No 'psc_roi*.txt' files under:\n  {root}\n"
              f"Tip: pass --generate-psc to build them from roi_*.txt first.")
        return

    if args.verbose:
        print(f"[INFO] Total PSC files found: {psc_count}")
        print(f"[INFO] Total groups: {len(groups)}")

    # Color & marker maps per group (consistent colors across Combined & Scatter)
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color',
                    ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9'])
    sorted_group_keys = [k for k,_ in sorted(groups.items(), key=lambda kv: kv[0].lower())]
    group_colors  = {gk: color_cycle[i % len(color_cycle)] for i, gk in enumerate(sorted_group_keys)}
    marker_cycle  = ['o','X','^','D','v','P','*','h','>','<','p']  # one per group; no squares
    group_markers = {gk: marker_cycle[i % len(marker_cycle)] for i, gk in enumerate(sorted_group_keys)}

    # For scatter stats
    scatter_data = []  # (group_key, color, marker, [per-file max], group_summary_value)

    combined_series = []   # for Combined plot
    group_stats = []

    for gkey in sorted_group_keys:
        items = groups[gkey]
        if args.verbose:
            print(f"\n[GROUP] {gkey}  (#files={len(items)})")
        series_list, names, lengths = [], [], []

        # Load each PSC file → 1 series/file
        for path, colname in sorted(items, key=lambda x: x[1].lower()):
            try:
                s = load_psc_series(path)
            except Exception as e:
                print(f"[SKIP] Read error: {path} -> {e}")
                continue
            series_list.append(s)
            names.append(colname)
            lengths.append(len(s))
            if args.verbose:
                rel = os.path.relpath(path, root)
                print(f"        + {rel} (N={len(s)})")

        if not series_list:
            if args.verbose:
                print("[WARN] No readable series in this group; skipping.")
            continue

        # Align lengths within group
        if len(set(lengths)) != 1:
            if args.strict_length:
                print(f"[ERROR] Unequal lengths in group '{gkey}': {lengths}. Use --strict-length off to truncate.")
                continue
            min_len = min(lengths)
            series_list = [s[:min_len] for s in series_list]
            if args.verbose:
                print(f"[INFO] Truncated all series to min length = {min_len}")
        else:
            min_len = lengths[0]

        # Matrix & outputs
        mat = np.column_stack(series_list)   # shape: (min_len, M)
        M = mat.shape[1]
        safe_group = sanitize(gkey)

        # Per-group CSV (matrix)
        csv_path = os.path.join(outdir, f"{safe_group}.csv")
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["index"] + names)
            for i in range(min_len):
                w.writerow([i] + [f"{val:.6f}" for val in mat[i, :]])

        # Moving-mean matrix CSV (column-wise movmean)
        if args.movmean and args.movmean > 1:
            sm_mat = np.column_stack([movmean_centered(mat[:, j], args.movmean) for j in range(M)])
        else:
            sm_mat = mat.copy()
        movmean_csv_path = os.path.join(outdir, f"mov_mean_{safe_group}.csv")
        with open(movmean_csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["index"] + names)
            for i in range(min_len):
                w.writerow([i] + [f"{val:.6f}" for val in sm_mat[i, :]])

        # Average and SEM across columns
        avg = np.nanmean(mat, axis=1)        # (min_len,)
        sem = None
        if M >= 2:
            sem = np.nanstd(mat, axis=1, ddof=1) / np.sqrt(M)
            sem_txt = os.path.join(outdir, f"sem_{safe_group}.txt")
            np.savetxt(sem_txt, sem, fmt="%.6f", delimiter="\t")
        avg_txt = os.path.join(outdir, f"average_{safe_group}.txt")
        np.savetxt(avg_txt, avg, fmt="%.6f", delimiter="\t")

        # Time axis (minutes)
        tmin = minutes_axis(min_len, args.tr)

        # Plot average series (+ SEM if present)
        import matplotlib.pyplot as plt
        avg_svg = os.path.join(outdir, f"Average_time_Course_{safe_group}.svg")
        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=150)
        line_avg, = ax.plot(tmin, avg, linewidth=3, label="Average")
        if sem is not None:
            ax.fill_between(tmin, avg - sem, avg + sem, alpha=0.2,
                            color=line_avg.get_color(), label="_nolegend_")
        ax.axvspan(10.0, 20.0, alpha=0.15, color='0.2')
        ax.set_xlabel("Time (in minutes)", fontsize=40)
        ax.set_ylabel("Percent Signal Change (%)", fontsize=40)
        ax.set_title(f"Average Time Course — {gkey}")
        if sem is not None:
            ax.legend(loc="best", frameon=False)
        apply_x_limits(ax)
        style_axes(ax)
        fig.tight_layout()
        fig.savefig(avg_svg, format="svg")
        plt.close(fig)

        # Movmean (optional)
        mov_txt = mov_svg = mov_sem_txt = ""
        mov_avg = mov_sem = None
        if args.movmean and args.movmean > 1:
            mov_avg = movmean_centered(avg, args.movmean)
            mov_txt = os.path.join(outdir, f"Movmean_Average_{safe_group}.txt")
            np.savetxt(mov_txt, mov_avg, fmt="%.6f", delimiter="\t")
            if sem is not None:
                mov_sem = movmean_centered(sem, args.movmean)
                mov_sem_txt = os.path.join(outdir, f"Movmean_sem_{safe_group}.txt")
                np.savetxt(mov_sem_txt, mov_sem, fmt="%.6f", delimiter="\t")

            mov_svg = os.path.join(outdir, f"Movmean_Average_{safe_group}.svg")
            fig, ax = plt.subplots(figsize=(20, 12), dpi=1200)
            line_mov, = ax.plot(tmin, mov_avg, linewidth=3, label="Movmean Average", color=group_colors[gkey])
            if mov_sem is not None:
                ax.fill_between(tmin, mov_avg - mov_sem, mov_avg + mov_sem, alpha=0.2,
                                color=line_mov.get_color(), label="_nolegend_")
            ax.axvspan(10.0, 20.0, alpha=0.15, color='0.2')
            ax.set_xlabel("Time (in minutes)", fontsize=40)
            ax.set_ylabel("% Signal Change", fontsize=30)
            ax.set_title(f"Movmean (w={args.movmean}) — {gkey}")
            ax.legend(loc="best", frameon=False)
            apply_x_limits(ax)
            apply_mov_y_limits(ax)  # y-limits only on Movmean plots
            style_axes(ax)
            fig.tight_layout()
            fig.savefig(mov_svg, format="svg")
            plt.close(fig)

            combined_series.append((gkey, tmin, mov_avg, mov_sem))

        # ----- Scatter stats (20–30 min window) -----
        mask_20_30 = (tmin >= 20.0) & (tmin <= 30.0)
        per_file_max = [float(np.nanmax(sm_mat[mask_20_30, j])) for j in range(sm_mat.shape[1])]
        if len(series_list) >= 2:
            if args.movmean and args.movmean > 1:
                if mov_avg is None:
                    mov_avg = movmean_centered(avg, args.movmean)
                group_summary = float(np.nanmax(mov_avg[mask_20_30]))
            else:
                group_summary = float(np.nanmax(avg[mask_20_30]))
        else:
            group_summary = per_file_max[0]
        scatter_data.append((gkey, group_colors[gkey], group_markers[gkey], per_file_max, group_summary))
        # -------------------------------------------

        group_stats.append((gkey, len(series_list), min_len, csv_path, avg_txt, avg_svg, mov_txt, mov_svg, (sem is not None)))
        if args.verbose:
            print(f"[OK] CSV  -> {csv_path}")
            print(f"[OK] MOV_MEAN CSV -> {movmean_csv_path}")
            print(f"[OK] AVG  -> {avg_txt}")
            if sem is not None: print(f"[OK] SEM  -> {os.path.join(outdir, f'sem_{safe_group}.txt')}")
            print(f"[OK] AVG SVG -> {avg_svg}")
            if mov_txt: print(f"[OK] MOV AVG TXT -> {mov_txt}")
            if mov_sem is not None: print(f"[OK] MOV SEM TXT -> {os.path.join(outdir, f'Movmean_sem_{safe_group}.txt')}")
            if mov_svg: print(f"[OK] MOV AVG SVG -> {mov_svg}")

    # Combined Movmean plot
    if combined_series:
        combined_svg = os.path.join(outdir, "Combined_Movmean_Averages.svg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(20, 12), dpi=1200)
        for gkey, tmin, mov_avg, mov_sem in sorted(combined_series, key=lambda x: x[0].lower()):
            col = group_colors[gkey]
            ax.plot(tmin, mov_avg, linewidth=3, label=gkey, color=col)
            if mov_sem is not None:
                ax.fill_between(tmin, mov_avg - mov_sem, mov_avg + mov_sem, alpha=0.12, color=col, label='_nolegend_')
        ax.axvspan(10.0, 20.0, alpha=0.15, color='0.2')
        ax.set_xlabel("Time (in minutes)", fontsize=40)
        ax.set_ylabel("% Signal Change", fontsize=30)
        ax.set_title("Combined Movmean Averages (per group)")
        ax.legend(loc="best", frameon=False, ncol=1)
        apply_x_limits(ax)
        apply_mov_y_limits(ax)
        style_axes(ax)
        fig.tight_layout()
        fig.savefig(combined_svg, format="svg")
        plt.close(fig)
        print(f"[OK] Combined movmean plot -> {combined_svg}")
    else:
        print("[INFO] No movmean provided or produced; combined plot not created.")

    # Scatter plot of 20–30 min maxima per group (+ group max as flat line)
    if scatter_data:
        scatter_svg = os.path.join(outdir, "Scatter_Max_20_30min_by_Group.svg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(20, 12), dpi=1200)

        x_positions = np.arange(1, len(scatter_data) + 1)  # 1..G

        for i, (gkey, color, marker, maxima, group_summary) in enumerate(scatter_data, start=1):
            for y in maxima:
                ax.scatter(i, y, marker=marker, s=140, color=color, linewidths=2, edgecolors='none')
            ax.hlines(y=group_summary, xmin=i-0.22, xmax=i+0.22, colors=color, linewidth=2)

        ax.set_xticks(x_positions)
        ax.set_xticklabels([gkey for gkey, _, _, _, _ in scatter_data], rotation=25, ha='right', fontsize=24)
        ax.set_ylabel("% Signal Change", fontsize=30)
        ax.set_title("Max (20–30 min): Per-file vs Group-Average Peak", fontsize=32)
        apply_mov_y_limits(ax)
        style_axes(ax)
        ax.margins(x=0.02)

        from matplotlib.lines import Line2D
        legend_elems = [
            Line2D([0], [0], marker='o', color='black', linestyle='None', markersize=10,
                   label='Per-file Max (20–30 min)', markerfacecolor='black'),
            Line2D([0], [0], color='black', linestyle='-', linewidth=2,
                   label='Group-Average Max (20–30 min)'),
        ]
        ax.legend(handles=legend_elems, loc='best', frameon=False)

        fig.tight_layout()
        fig.savefig(scatter_svg, format="svg")
        plt.close(fig)
        print(f"[OK] Scatter (max 20–30 min) -> {scatter_svg}")

    # Final summary
    print("\n=== GROUP BUILD COMPLETE ===")
    print("Output folder:", outdir)
    if not group_stats:
        print("[WARN] No groups were produced (check input files).")
        return
    print("\nGroup\t#Files\tLength\tMatrix CSV\tMovMean CSV\tAverage TXT\tAverage SVG\tMovmean AVG TXT\tMovmean AVG SVG\tHas SEM?")
    for (gkey, M, L, csvp, avgt, avgs, mtxt, msvg, has_sem) in group_stats:
        movmean_csv = os.path.join(outdir, f"mov_mean_{sanitize(gkey)}.csv")
        print(f"{gkey}\t{M}\t{L}\t{csvp}\t{movmean_csv}\t{avgt}\t{avgs}\t{mtxt}\t{msvg}\t{has_sem}")

if __name__ == "__main__":
    main()
