#!/usr/bin/env python3
"""
Static-map movie from a 4D NIfTI using a sliding 200-volume window.

New in this version
-------------------
--rotate {0,90,180,270} : rotate each slice at DISPLAY time (CCW). Saved 4D NIfTI (if any) is NOT rotated.

Dependencies
------------
pip install nibabel numpy matplotlib imageio imageio-ffmpeg tqdm
"""
import argparse, os, math, tempfile
import numpy as np
import nibabel as nib
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
import imageio.v2 as iio  # v2 API (FFmpeg backend)



# -------------------------------------------------------------------
# Custom "blackbody" colormap registration (black → red → orange → yellow → white)
# -------------------------------------------------------------------
if "blackbody" not in plt.colormaps():
    blackbody = LinearSegmentedColormap.from_list(
        "blackbody",
        [
            (0.00, "#000000"),  # black
            (0.25, "#1c0000"),  # very dark red
            (0.45, "#6b0000"),  # deep red
            (0.65, "#ff3c00"),  # orange-red
            (0.82, "#ffae00"),  # yellow-orange
            (1.00, "#ffffff"),  # white
        ]
    )
    plt.register_cmap("blackbody", blackbody)
    
def parse_args():
    p = argparse.ArgumentParser(
        description="Make sliding-window static-map movie from 4D NIfTI (axial montage)."
    )
    p.add_argument("nii4d", help="Input 4D NIfTI (e.g., cleaned_mc_func.nii.gz)")

    # Data / math
    p.add_argument("--baseline-start", type=int, required=True,
                   help="Baseline start vol (0-based, inclusive)")
    p.add_argument("--baseline-end", type=int, required=True,
                   help="Baseline end vol (0-based, inclusive)")
    p.add_argument("--underlay-first-n", type=int, default=600,
                   help="Average first N vols for underlay (default: 600)")
    p.add_argument("--window", type=int, default=200,
                   help="Sliding window size (default: 200)")
    p.add_argument("--mask-zero", action="store_true",
                   help="Mask baseline==0 voxels to avoid divide-by-zero")
    p.add_argument("--save-nii", action="store_true",
                   help="Also save a 4D NIfTI of all sliding static maps (can be HUGE)")

    # Display / rendering
    p.add_argument("--mode", choices=["both", "pos", "neg"], default="both",
                   help="Overlay mode: both signs, positive-only, or negative-only")
    p.add_argument("--vmax", type=float, default=20.0,
                   help="Color limit for |pct-change| (default: 20)")
    p.add_argument("--cmap", default="seismic",
                   help="Matplotlib colormap for overlay (default: seismic)")
    p.add_argument("--cols", type=int, default=0,
                   help="Montage columns (default: 8 if slices>=8; else=#slices)")
    p.add_argument("--figw", type=float, default=14.0,
                   help="Figure width in inches (default: 14.0)")
    p.add_argument("--dpi", type=int, default=120,
                   help="Figure DPI (default: 120)")
    p.add_argument("--rotate", type=int, choices=[0,90,180,270], default=0,
                   help="Rotate each slice counter-clockwise by this many degrees (display only)")

    # Movie
    p.add_argument("--fps", type=int, default=10,
                   help="Frames per second for MP4 (default: 10)")
    p.add_argument("--tr", type=float, default=None,
                   help="TR in seconds for time label (optional)")
    p.add_argument("--out-prefix", default="Static_Map",
                   help="Output prefix (default: Static_Map)")
    p.add_argument("--no-video", action="store_true",
                   help="Skip MP4 creation")

    return p.parse_args()

def montage_grid(n_slices, user_cols):
    """
    If user_cols>0, honor it.
    Otherwise use 8 columns when n_slices>=8; else use n_slices columns.
    """
    if user_cols and user_cols > 0:
        cols = min(user_cols, n_slices)
    else:
        cols = 8 if n_slices >= 8 else n_slices
    cols = max(cols, 1)
    rows = int(math.ceil(n_slices / cols))
    return rows, cols

def compute_underlay_limits(underlay):
    # Robust grayscale limits (2nd–98th percentile)
    finite = underlay[np.isfinite(underlay)]
    lo = float(np.percentile(finite, 2)) if finite.size else 0.0
    hi = float(np.percentile(finite, 98)) if finite.size else 1.0
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = float(np.nanmin(underlay)), float(np.nanmax(underlay))
        if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
            lo, hi = 0.0, 1.0
    return lo, hi

def _rot_k_from_degrees(deg: int) -> int:
    # np.rot90 rotates CCW by 90° per k
    return (deg // 90) % 4

def _rotate2d(arr2d, k: int):
    return np.rot90(arr2d, k) if k else arr2d

def main():
    a = parse_args()
    rot_k = _rot_k_from_degrees(a.rotate)

    # Load data
    img = nib.load(a.nii4d)
    data = img.get_fdata(dtype=np.float32)
    if data.ndim != 4:
        raise ValueError("Input must be 4D.")
    X, Y, Z, T = data.shape

    # Underlay: avg first N vols
    if not (1 <= a.underlay_first_n <= T):
        raise ValueError(f"--underlay-first-n must be in [1, {T}]")
    underlay = np.nanmean(data[..., :a.underlay_first_n], axis=3).astype(np.float32)

    # Baseline: mean of [bs..be]
    bs, be = a.baseline_start, a.baseline_end
    if not (0 <= bs <= be < T):
        raise ValueError(f"Invalid baseline range: start={bs}, end={be}, T={T}")
    baseline = np.nanmean(data[..., bs:be+1], axis=3).astype(np.float32)
    base_mask = (baseline != 0) if a.mask_zero else np.ones_like(baseline, bool)

    # Sliding windows
    W = int(a.window)
    if W <= 0 or W > T:
        raise ValueError("--window must be in [1, T]")
    n_frames = T - W + 1  # sliding: inclusive
    print(f"Sliding window: {W} vols, frames: {n_frames} (0..{n_frames-1})")

    # NIfTI saving (optional, CAUTION: can be very large)
    if a.save_nii:
        maps_4d = np.zeros((X, Y, Z, n_frames), dtype=np.float32)

    # Prepare display scaling
    umin, umax = compute_underlay_limits(underlay)
    # Overlay norm
    if a.mode == "both":
        ov_norm = Normalize(vmin=-abs(a.vmax), vmax=abs(a.vmax))
    elif a.mode == "pos":
        ov_norm = Normalize(vmin=0.0, vmax=abs(a.vmax))
    else:  # "neg"
        ov_norm = Normalize(vmin=-abs(a.vmax), vmax=0.0)

    # If only NIfTI requested, compute and exit early
    if a.no_video and not a.save_nii:
        for k in tqdm(range(n_frames), ncols=80, desc="Computing maps"):
            s, e = k, k + W
            avgW = np.nanmean(data[..., s:e], axis=3)
            stat = np.full_like(avgW, np.nan, dtype=np.float32)
            np.divide((avgW - baseline) * 100.0, baseline, out=stat, where=base_mask)

            if a.mode == "pos":
                stat[stat < 0] = np.nan
            elif a.mode == "neg":
                stat[stat > 0] = np.nan
        print("Done (no outputs requested except computation).")
        return

    # Setup video writer (unless --no-video)
    writer = None
    if not a.no_video:
        mp4 = f"{a.out_prefix}_sliding_movie.mp4"
        try:
            writer = iio.get_writer(mp4, fps=a.fps, codec="libx264", quality=8)
        except Exception:
            writer = iio.get_writer(mp4, fps=a.fps)  # fallback
        print(f"Writing video to: {mp4} (fps={a.fps})")

    # Figure layout — enforce 8 columns (or fewer if Z<8)
    rows, cols = montage_grid(Z, a.cols)
    figw = a.figw
    tile_h = (figw / max(cols, 1))         # keep tiles ~square
    cbar_frac = 0.18                        # fraction of one tile for colorbar row
    figh = tile_h * (rows + cbar_frac) + 1.2  # extra headroom for suptitle

    # Colorbar mappable
    sm = ScalarMappable(norm=ov_norm, cmap=a.cmap)

    # Main loop
    for k in tqdm(range(n_frames), ncols=80, desc="Rendering"):
        s, e = k, k + W
        avgW = np.nanmean(data[..., s:e], axis=3)

        # Percent change
        stat = np.full_like(avgW, np.nan, dtype=np.float32)
        np.divide((avgW - baseline) * 100.0, baseline, out=stat, where=base_mask)

        # Apply mode
        if a.mode == "pos":
            stat[stat < 0] = np.nan
        elif a.mode == "neg":
            stat[stat > 0] = np.nan

        # Save into 4D stack if requested (NOTE: not rotated)
        if a.save_nii:
            maps_4d[..., k] = np.nan_to_num(stat, nan=0.0, posinf=0.0, neginf=0.0)

        # Skip frame drawing if no video
        if a.no_video:
            continue

        # -------- draw frame (all z slices montage) --------
        fig = plt.figure(figsize=(figw, figh), constrained_layout=False)

        # GridSpec: rows of slices + 1 short row for the horizontal colorbar
        height_ratios = [1.0] * rows + [cbar_frac]
        gs = fig.add_gridspec(rows + 1, cols, wspace=0.02, hspace=0.02,
                              height_ratios=height_ratios)

        # Time label
        if a.tr is not None:
            t0 = s * a.tr
            t1 = (e - 1) * a.tr
            subtitle = f"Window: vols {s+1}–{e}   Time: {t0:.2f}–{t1:.2f} s"
        else:
            subtitle = f"Window: vols {s+1}–{e}   Frame: {k+1}/{n_frames}"

        fig.suptitle(
            "Percent Change ((Mean(window) − Baseline) / Baseline) × 100\n" + subtitle,
            fontsize=11, fontweight="bold", y=0.98
        )

        # Draw slices (rows x cols)
        zi = 0
        for r in range(rows):
            for c in range(cols):
                ax = fig.add_subplot(gs[r, c])
                ax.axis("off")
                if zi < Z:
                    # Underlay (grayscale)
                    u = underlay[:, :, zi]
                    u = _rotate2d(u, rot_k)
                    ax.imshow(u, cmap="gray", vmin=umin, vmax=umax,
                              origin="lower", interpolation="nearest")

                    # Overlay (percent change)
                    o = stat[:, :, zi]
                    o = _rotate2d(o, rot_k)
                    o_draw = np.ma.masked_invalid(o)
                    ax.imshow(o_draw, cmap=a.cmap, norm=ov_norm, alpha=0.7,
                              origin="lower", interpolation="nearest")

                    ax.set_title(f"z={zi}", fontsize=8, pad=1.0)
                    zi += 1

        # Horizontal colorbar on its own full-width row (below the montage)
        cbar_ax = fig.add_subplot(gs[rows, :])
        cb = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
        if a.mode == "both":
            cb.set_label("Pct change (±)", fontsize=9)
        elif a.mode == "pos":
            cb.set_label("Pct change (+)", fontsize=9)
        else:
            cb.set_label("Pct change (−)", fontsize=9)

        # Render to PNG in-memory then append
        tmp_png = None
        try:
            tmp_png = os.path.join(tempfile.gettempdir(), f"_frame_{k:06d}.png")
            fig.savefig(tmp_png, dpi=a.dpi, bbox_inches="tight")
            plt.close(fig)
            writer.append_data(iio.imread(tmp_png))
        finally:
            if tmp_png and os.path.exists(tmp_png):
                try:
                    os.remove(tmp_png)
                except Exception:
                    pass

    # Finish up
    if writer is not None:
        writer.close()

    if a.save_nii:
        out_nii = f"{a.out_prefix}_sliding_pct_static_maps.nii.gz"
        nib.save(nib.Nifti1Image(maps_4d, img.affine, img.header), out_nii)
        print(f"Wrote: {out_nii}")

    print("Done.")

if __name__ == "__main__":
    main()
