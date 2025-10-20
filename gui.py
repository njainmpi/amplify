# app.py
import io, os, csv, glob, shutil, pathlib, urllib.parse
from typing import List, Tuple, Dict, Any, Optional

import numpy as np
import pandas as pd
import nibabel as nib
import streamlit as st
from string import Template
from streamlit.components.v1 import html as st_html
import matplotlib.cm as cm

st.set_page_config(page_title="fMRI Pipeline GUI", layout="wide")

# ==== PORTAL CHROME (header/footer text) ====
PORTAL_TITLE = "Database Management Portal"
PORTAL_VERSION = "v0.1"
FOOTER_LEFT = "© 2025 — Your Lab"
FOOTER_RIGHT = "Root: {root}"

# ==== PERMANENT ROOTS ====
ROOT = "/run/user/6631/gvfs/smb-share:server=sk3,share=pr_ohlendorf/fMRI"
DEFAULT_CSV_PATH = f"{ROOT}/RawData/Animal_Experiments_Sequences_v1.csv"
DEFAULT_ANALYSED_ROOT = f"{ROOT}/AnalysedData"

# ==== UI CONSTANTS ====
HEADER_LINES_DEFAULT = 2
FIRST8_HEADERS = ["Col1","Col2","Col3","Col4","Col5","Col6","Col7","Col8"]
TILE_COLS = 8
TILE_ROWS = 4
PROJECT_TILE_COLS = 4
MATCH_TILE_COLS = 4

EXPECTED_FILES = [
    "mc_func.nii.gz",
    "mask_mean_mc_func.nii.gz",
    "cleaned_mc_func.nii.gz",
    "cleaned_anatomy.nii.gz",   # (special Col7 location rule)
    "sm_fMRI_for_scm.nii.gz",
    "mask_mean_fMRI_coregistered_to_struct.nii.gz",
    "sm_coreg_func_Static_Map_*.nii.gz",
]

# ---------- helpers ----------
def expand_path(p: str) -> str:
    return os.path.abspath(os.path.expanduser(os.path.expandvars(p or ".")))

def try_decode(raw: bytes) -> str:
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")

def read_csv_any(path_or_buffer, header_lines: int) -> Tuple[pd.DataFrame, List[List[str]]]:
    if hasattr(path_or_buffer, "read"):
        raw = path_or_buffer.read()
        text = try_decode(raw)
    else:
        with open(path_or_buffer, "rb") as f:
            text = try_decode(f.read())
    rows = list(csv.reader(io.StringIO(text)))
    df = pd.DataFrame(rows)
    ncols = max(8, df.shape[1])
    if df.shape[1] < ncols:
        for _ in range(ncols - df.shape[1]):
            df[_] = ""
    return df, rows

def build_preview(rows: List[List[str]], header_lines: int) -> pd.DataFrame:
    preview_records = []
    for i in range(header_lines, len(rows)):
        r = rows[i]
        if not r or all((c is None) or (str(c).strip() == "") for c in r):
            continue
        first8 = [(r[k] if k < len(r) else "").replace("|", "¦") for k in range(8)]
        preview_records.append({"CSV Line": i+1, **{FIRST8_HEADERS[k]: first8[k] for k in range(8)}})
    return pd.DataFrame(preview_records)

def extract_meta(row: List[str]) -> Dict[str, Any]:
    get = lambda n: (row[n-1].strip() if n-1 < len(row) and row[n-1] is not None else "")
    return {
        "Project_Name": get(3),
        "Sub_project_Name": get(4),
        "Dataset_Name": get(2),
        "Structural_Run": get(5),
        "Functional_Run": get(6),
        "Structural_for_CoReg": get(7),
        "Baseline_Duration_min": get(8),
        "Injection_Duration_min": get(9),
        "Injected_Left": get(22),
        "Injected_Right": get(23),
    }

def meta_field(label: str, value: str):
    c1, c2 = st.columns([0.7, 1.3])
    c1.caption(label)
    c2.write(value if value else "—")

def unique_projects(rows: List[List[str]], header_lines: int) -> List[str]:
    projects = []
    for i in range(header_lines, len(rows)):
        r = rows[i]
        if not r:
            continue
        val = (r[2].strip() if len(r) > 2 and r[2] is not None else "")
        if val:
            projects.append(val)
    return sorted(set(projects), key=str.lower)

def lines_for_project(rows: List[List[str]], header_lines: int, project: Optional[str]) -> List[int]:
    line_nums = []
    for i in range(header_lines, len(rows)):
        r = rows[i]
        if not r:
            continue
        if project in (None, ""):
            line_nums.append(i+1)
        else:
            val = (r[2].strip() if len(r) > 2 and r[2] is not None else "")
            if val == project:
                line_nums.append(i+1)
    return line_nums

def analysed_base_path(root: str, meta: Dict[str, Any]) -> str:
    return os.path.join(
        root, "AnalysedData",
        meta.get("Project_Name",""),
        meta.get("Sub_project_Name",""),
        meta.get("Dataset_Name",""),
    )

def find_matches_containing_token(base: str, token_value: str) -> list[str]:
    matches: list[str] = []
    token = (token_value or "").strip().lower()
    try:
        if os.path.isdir(base):
            for entry in os.scandir(base):
                if entry.is_dir():
                    name_lower = entry.name.lower()
                    if (not token) or (token in name_lower):
                        matches.append(entry.path)
    except PermissionError:
        pass
    return sorted(matches)

def list_subdirs(path: str) -> list[str]:
    dirs: list[str] = []
    try:
        if os.path.isdir(path):
            for entry in os.scandir(path):
                if entry.is_dir():
                    dirs.append(entry.path)
    except PermissionError:
        pass
    return sorted(dirs)

def check_expected_file_at(base_dir: str, pattern: str) -> Dict[str, Any]:
    abs_pat = os.path.join(base_dir, pattern)
    if any(ch in pattern for ch in "*?[]"):
        matches = sorted(glob.glob(abs_pat))
        return {"pattern": pattern, "base_dir": base_dir, "found": len(matches) > 0, "matches": matches}
    else:
        exists = os.path.isfile(abs_pat)
        return {"pattern": pattern, "base_dir": base_dir, "found": exists, "matches": [abs_pat] if exists else []}

# ---------- Optional static serving for WebGL ----------
STATIC_DIR = pathlib.Path(__file__).parent / "static" / "niis"
STATIC_DIR.mkdir(parents=True, exist_ok=True)

def serve_local_file_as_static_url(src_path: str) -> str:
    src = pathlib.Path(src_path)
    dst = STATIC_DIR / src.name
    try:
        if dst.exists():
            if dst.stat().st_size == src.stat().st_size and int(dst.stat().st_mtime) == int(src.stat().st_mtime):
                pass
            else:
                dst.unlink()
        if not dst.exists():
            try:
                os.link(src, dst)
            except Exception:
                try:
                    dst.symlink_to(src)
                except Exception:
                    shutil.copy2(src, dst)
    except Exception:
        shutil.copy2(src, dst)
    return "/static/niis/" + urllib.parse.quote(src.name)

# ---------- NIfTI I/O ----------
@st.cache_resource(show_spinner=False)
def _open_nifti_lazy(path: str):
    return nib.load(path, mmap=True)

def _get_file_size_mb(path: str) -> float:
    try:
        return os.path.getsize(path) / (1024*1024)
    except Exception:
        return 0.0

def _slice_from_arrayproxy(proxy, axis: str, index: int, vol_index: Optional[int], full_shape) -> np.ndarray:
    """Return 2D slice from 3D/4D NIfTI, raw (no interpolation)."""
    if len(full_shape) >= 4:
        if vol_index is None:
            vol_index = 0
        vol_view = proxy[..., vol_index]  # XYZ
    else:
        vol_view = proxy

    if len(full_shape) == 2:
        arr = np.asanyarray(vol_view)
        return arr.astype(np.float32, copy=False)

    if axis == "axial":      # Z
        arr = np.asanyarray(vol_view[:, :, index])
    elif axis == "coronal":  # Y
        arr = np.asanyarray(vol_view[:, index, :])
        arr = np.rot90(arr, k=1)[:, ::-1]
    elif axis == "sagittal": # X
        arr = np.asanyarray(vol_view[index, :, :])
        arr = np.rot90(arr, k=1)
    else:
        arr = np.asanyarray(vol_view[:, :, index])

    return arr.astype(np.float32, copy=False)

def _auto_window_from_slices(slices: List[np.ndarray], p_low: float, p_high: float) -> Tuple[float, float]:
    all_vals = []
    max_collect = 3_000_000
    for sl in slices:
        vals = sl[np.isfinite(sl)]
        if vals.size == 0:
            continue
        cap = max_collect // max(1, len(slices))
        if vals.size > cap:
            choice = np.random.choice(vals.size, size=cap, replace=False)
            vals = vals[choice]
        all_vals.append(vals)
    if all_vals:
        all_vals = np.concatenate(all_vals)
        vmin = float(np.percentile(all_vals, p_low))
        vmax = float(np.percentile(all_vals, p_high))
        if vmin == vmax:
            vmax = vmin + 1.0
        return vmin, vmax
    return 0.0, 1.0

# ---------- CPU multi-slice viewer ----------
def render_cpu_slice_viewer(path: str):
    st.markdown("<style>.stImage img{image-rendering: pixelated !important;}</style>", unsafe_allow_html=True)

    size_mb = _get_file_size_mb(path)
    img = _open_nifti_lazy(path)
    shape = img.shape
    hdr = img.header
    dtype = str(hdr.get_data_dtype())
    prox = img.dataobj

    st.write(
        f"**File:** `{os.path.basename(path)}`  •  **Size:** {size_mb:.0f} MB  •  "
        f"**Shape:** {shape}  •  **dtype:** {dtype}"
    )

    vol_index = None
    if len(shape) >= 4:
        tmax = shape[3] - 1
        vol_index = st.slider("Time / Volume (T)", 0, int(tmax), 0, 1, key="cpu_vol")

    axis = st.radio("View", ["axial", "coronal", "sagittal"], horizontal=True, index=0, key="cpu_axis")

    if len(shape) == 2:
        spatial = (shape[0], shape[1], 1)
    else:
        spatial = shape[:3]

    if axis == "axial":
        max_idx = max(0, spatial[2] - 1)
    elif axis == "coronal":
        max_idx = max(0, spatial[1] - 1)
    else:
        max_idx = max(0, spatial[0] - 1)

    idx = st.slider("Center slice index", 0, int(max_idx), int(max_idx // 2) if max_idx > 0 else 0, key="cpu_center_idx")

    c1, c2, c3 = st.columns([1,1,1])
    with c1:
        rot_deg = st.selectbox("Rotate", [0, 90, 180, 270], index=0, key="cpu_rotate")
    with c2:
        cols_in_grid = st.number_input("Grid columns", 1, 12, 8, 1, key="cpu_cols")
    with c3:
        rows_in_grid = st.number_input("Grid rows", 1, 8, 2, 1, key="cpu_rows")

    step = st.number_input("Slice step", 1, max(1, int(max_idx) or 1), 1, 1, key="cpu_step")

    n_slices = int(cols_in_grid * rows_in_grid)
    step_int = int(step)
    if n_slices == 1:
        start = min(max(idx, 0), max_idx)
    else:
        half = (n_slices - 1) // 2
        start = idx - half * step_int
        start = max(0, min(start, max_idx - (n_slices - 1) * step_int))
    indices = start + np.arange(n_slices) * step_int
    indices = np.clip(indices, 0, max_idx).astype(int)

    slices = []
    k_rot = (rot_deg // 90) % 4
    for sidx in indices:
        arr2d = _slice_from_arrayproxy(prox, axis, int(sidx), vol_index, shape)
        if k_rot:
            arr2d = np.rot90(arr2d, k=k_rot)
        slices.append(arr2d)

    vmin_init, vmax_init = _auto_window_from_slices(slices, 2.0, 98.0)

    i1, i2 = st.columns(2)
    with i1:
        vmin = st.number_input("Lower intensity", value=float(vmin_init), format="%.6g", key="cpu_vmin")
    with i2:
        vmax = st.number_input("Upper intensity", value=float(vmax_init), format="%.6g", key="cpu_vmax")

    if vmax <= vmin:
        vmax = vmin + 1e-6

    norm_slices = [np.clip((sl - vmin) / max(1e-12, (vmax - vmin)), 0.0, 1.0) for sl in slices]

    st.caption(
        f"{n_slices} slices • axis **{axis}** "
        + (f"• T={vol_index} " if vol_index is not None else "")
        + f"• window [{vmin:.3g}, {vmax:.3g}] • rotation {rot_deg}°"
    )

    ptr = 0
    for _ in range(int(rows_in_grid)):
        cols = st.columns(int(cols_in_grid), gap="small")
        for c in cols:
            if ptr >= len(norm_slices):
                c.empty()
            else:
                sidx = int(indices[ptr])
                c.image(
                    norm_slices[ptr],
                    caption=f"{axis} {sidx}",
                    use_container_width=True,
                    clamp=True,
                    channels="L",
                    output_format="PNG",
                )
            ptr += 1

# ---------- Overlay viewer (multi-slice, same controls) ----------
def render_overlay_viewer_grid(base_path: str, overlay_path: str):
    """Overlay one NIfTI on another in a multi-slice grid (no interpolation)."""
    st.markdown("<style>.stImage img{image-rendering: pixelated !important;}</style>", unsafe_allow_html=True)

    base_img = _open_nifti_lazy(base_path)
    over_img = _open_nifti_lazy(overlay_path)

    base_shape = base_img.shape
    over_shape = over_img.shape

    base_xyz = base_shape[:3] if len(base_shape) >= 3 else (base_shape[0], base_shape[1], 1)
    over_xyz  = over_shape[:3] if len(over_shape)  >= 3 else (over_shape[0],  over_shape[1],  1)

    if base_xyz != over_xyz:
        st.error(f"Spatial shapes differ: base {base_xyz} vs overlay {over_xyz}. Use coregistered images.")
        return

    st.write(
        f"**Base:** `{os.path.basename(base_path)}`  •  Shape: {base_shape}  "
        f"•  **Overlay:** `{os.path.basename(overlay_path)}`  •  Shape: {over_shape}"
    )

    base_t = st.slider("Base volume (if 4D)", 0, max(0, (base_shape[3]-1) if len(base_shape)>=4 else 0), 0, 1, key="ovlm_base_t") if len(base_shape)>=4 else None
    over_t = st.slider("Overlay volume (if 4D)", 0, max(0, (over_shape[3]-1) if len(over_shape)>=4 else 0), 0, 1, key="ovlm_over_t") if len(over_shape)>=4 else None

    axis = st.radio("Axis", ["axial", "coronal", "sagittal"], horizontal=True, index=0, key="ovlm_axis")

    if axis == "axial":
        max_idx = base_xyz[2] - 1
    elif axis == "coronal":
        max_idx = base_xyz[1] - 1
    else:
        max_idx = base_xyz[0] - 1

    center_idx = st.slider("Center slice index", 0, int(max_idx), int(max_idx//2) if max_idx>0 else 0, key="ovlm_center_idx")

    c1, c2, c3 = st.columns([1,1,1])
    with c1:
        rot_deg = st.selectbox("Rotate", [0, 90, 180, 270], index=0, key="ovlm_rotate")
    with c2:
        cols_in_grid = st.number_input("Grid columns", 1, 12, 8, 1, key="ovlm_cols")
    with c3:
        rows_in_grid = st.number_input("Grid rows", 1, 8, 2, 1, key="ovlm_rows")

    step = st.number_input("Slice step", 1, max(1, int(max_idx) or 1), 1, 1, key="ovlm_step")

    # windows, colormap, opacity, mask
    st.markdown("**Intensity windows** (base & overlay)")
    b1, b2 = st.columns(2)
    o1, o2 = st.columns(2)

    # sample one pair to auto-init window
    k_rot = (rot_deg // 90) % 4
    sample_b = _slice_from_arrayproxy(base_img.dataobj, axis, int(center_idx), base_t, base_shape)
    sample_o = _slice_from_arrayproxy(over_img.dataobj, axis, int(center_idx), over_t, over_shape)
    if k_rot:
        sample_b = np.rot90(sample_b, k=k_rot)
        sample_o = np.rot90(sample_o, k=k_rot)
    bmin_i, bmax_i = _auto_window_from_slices([sample_b], 2.0, 98.0)
    omin_i, omax_i = _auto_window_from_slices([sample_o], 2.0, 98.0)

    with b1:
        bmin = st.number_input("Base lower", value=float(bmin_i), format="%.6g", key="ovlm_bmin")
    with b2:
        bmax = st.number_input("Base upper", value=float(bmax_i), format="%.6g", key="ovlm_bmax")
    if bmax <= bmin: bmax = bmin + 1e-6

    with o1:
        omin = st.number_input("Overlay lower", value=float(omin_i), format="%.6g", key="ovlm_omin")
    with o2:
        omax = st.number_input("Overlay upper", value=float(omax_i), format="%.6g", key="ovlm_omax")
    if omax <= omin: omax = omin + 1e-6

    cmap_name = st.selectbox(
        "Overlay colormap",
        ["viridis","plasma","hot","turbo","magma","inferno","coolwarm","bwr","jet","gray"],
        index=0,
        key="ovlm_cmap"
    )
    opacity = st.slider("Overlay opacity", 0.0, 1.0, 0.6, 0.05, key="ovlm_opacity")
    mask_zeros = st.checkbox("Mask overlay zeros", value=True, key="ovlm_mask_zeros")

    # compute slice indices like CPU viewer
    n_slices = int(cols_in_grid * rows_in_grid)
    step_int = int(step)
    if n_slices == 1:
        start = min(max(center_idx, 0), max_idx)
    else:
        half = (n_slices - 1) // 2
        start = center_idx - half * step_int
        start = max(0, min(start, max_idx - (n_slices - 1) * step_int))
    indices = start + np.arange(n_slices) * step_int
    indices = np.clip(indices, 0, max_idx).astype(int)

    # render grid
    st.caption(
        f"{n_slices} slices • axis **{axis}** • window base [{bmin:.3g}, {bmax:.3g}] "
        f"overlay [{omin:.3g}, {omax:.3g}] • rotation {rot_deg}°"
    )

    ptr = 0
    for _ in range(int(rows_in_grid)):
        cols = st.columns(int(cols_in_grid), gap="small")
        for c in cols:
            if ptr >= len(indices):
                c.empty()
            else:
                sidx = int(indices[ptr])
                b2d = _slice_from_arrayproxy(base_img.dataobj, axis, sidx, base_t, base_shape)
                o2d = _slice_from_arrayproxy(over_img.dataobj, axis, sidx, over_t, over_shape)
                if k_rot:
                    b2d = np.rot90(b2d, k=k_rot)
                    o2d = np.rot90(o2d, k=k_rot)

                # base → grayscale
                b = np.clip((b2d - bmin) / max(1e-12, (bmax - bmin)), 0, 1)
                base_rgb = np.repeat(b[..., None], 3, axis=2)

                # overlay → cmap
                on = np.clip((o2d - omin) / max(1e-12, (omax - omin)), 0, 1)
                cmap = cm.get_cmap(cmap_name)
                over_rgba = cmap(on)

                alpha = over_rgba[..., 3] * opacity
                if mask_zeros:
                    alpha = np.where(o2d == 0, 0.0, alpha)

                out = over_rgba[..., :3] * alpha[..., None] + base_rgb * (1.0 - alpha[..., None])

                c.image(out, caption=f"{axis} {sidx}", use_container_width=True, clamp=True)
            ptr += 1

# ---------- WebGL viewer ----------
def render_nifti_viewer_url(url: str, name: str = "volume.nii.gz", height: int = 640):
    html = f"""
<div id="viewer-holder" style="width: 100%; height: {height}px; background: #000; border-radius: 8px; overflow: hidden; position: relative;">
  <div id="nv-wrap" style="position:absolute; inset:0; display:flex;">
    <canvas id="gl" style="width:100%; height:100%; display:block;"></canvas>
  </div>
</div>
<script src="https://unpkg.com/niivue@latest/dist/niivue.umd.js"></script>
<script>
(function(){{
  const url = "{url}";
  const name = "{name}";
  const holder = document.getElementById("viewer-holder");
  const canvas = document.getElementById("gl");
  function fit() {{
    const dpr = window.devicePixelRatio||1;
    canvas.width = Math.max(1, holder.clientWidth) * dpr;
    canvas.height = Math.max(1, holder.clientHeight) * dpr;
  }}
  try {{
    fit();
    window.addEventListener("resize", fit);
    const nv = new niivue.Niivue({{ isColorbar: true }});
    nv.attachToCanvas(canvas);
    nv.loadVolumes([{{url, name, colormap:"gray"}}]).then(() => {{
      nv.setSliceType(nv.sliceTypeMultiplanar);
      nv.drawScene();
    }});
  }} catch(e) {{
    console.warn("NiiVue init failed", e);
  }}
}})();
</script>
"""
    st_html(html, height=height+8, scrolling=False)

# ---------- header/footer ----------
def render_portal_chrome():
    if st.session_state.get("rows") and st.session_state.get("line_selected"):
        ln = st.session_state["line_selected"]
        row = st.session_state["rows"][ln - 1]
        meta = extract_meta(row)
        parts = []
        if meta.get("Project_Name"): parts.append(meta["Project_Name"])
        if meta.get("Sub_project_Name"): parts.append(meta["Sub_project_Name"])
        if meta.get("Dataset_Name"): parts.append(meta["Dataset_Name"])
        parts = " / ".join(parts)
        sel_status = f"Selected: Line {ln}" + (f" — {parts}" if parts else "")
    else:
        sel_status = "No selection"
    footer_text = f"{FOOTER_LEFT} • {FOOTER_RIGHT.format(root=ROOT)}"
    tpl = Template(r"""
<style>
html, body, #root, .stApp { margin:0!important; padding:0!important; }
.stApp > header, header[role="banner"] { display:none!important; }
:root { --header-h:92px; --footer-h:56px; --band-bg:#0f172a; --band-border:#1f2937; --text:#e5e7eb; --muted:#9ca3af; --dot-bg:#111827; --dot-hover:#1f2937; }
.block-container { padding-top:calc(var(--header-h) + 24px)!important; padding-bottom:calc(var(--footer-h) + 24px)!important; }
.app-header { position:fixed; top:0; left:0; right:0; z-index:1000; height:var(--header-h); background:var(--band-bg); border-bottom:1px solid var(--band-border); display:flex; align-items:center; }
.app-header .inner { max-width:1200px; width:100%; margin:0 auto; padding:0 20px; color:var(--text); display:flex; justify-content:space-between; align-items:center; font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, Cantarell, Noto Sans, Helvetica Neue, Arial; }
.app-header .title { font-weight:800; font-size:24px; }
.app-header .right { display:flex; align-items:center; gap:14px; opacity:.95; }
.app-header .status { font-size:13px; color:var(--muted); white-space:nowrap; overflow:hidden; text-overflow:ellipsis; max-width:52vw; }
.app-header .dots { display:inline-flex; align-items:center; justify-content:center; width:36px; height:36px; border-radius:8px; background:var(--dot-bg); color:var(--text); border:1px solid var(--band-border); font-size:20px; text-decoration:none; }
.app-header .dots:hover { background:var(--dot-hover); }
.app-footer { position:fixed; bottom:0; left:0; right:0; z-index:1000; height:var(--footer-h); background:#0f172a; border-top:1px solid var(--band-border); display:flex; align-items:center; }
.app-footer .inner { max-width:1200px; width:100%; margin:0 auto; padding:0 20px; color:var(--text); display:flex; justify-content:space-between; align-items:center; font-size:13px; opacity:.95; }
</style>
<div class="app-header"><div class="inner"><div class="title">${title}</div><div class="right"><div class="status">${status}</div><a href="#settings-panel" class="dots" title="Settings">⋮</a></div></div></div>
<div class="app-footer"><div class="inner"><div>${footer}</div><div>Streamlit UI</div></div></div>
""")
    st.markdown(tpl.substitute(title=PORTAL_TITLE, status=sel_status, footer=footer_text), unsafe_allow_html=True)

# ---------- state ----------
if "rows" not in st.session_state: st.session_state.rows = None
if "df_all" not in st.session_state: st.session_state.df_all = None
if "line_selected" not in st.session_state: st.session_state.line_selected = None
if "header_lines" not in st.session_state: st.session_state.header_lines = HEADER_LINES_DEFAULT
if "loaded_name" not in st.session_state: st.session_state.loaded_name = None
if "tile_page" not in st.session_state: st.session_state.tile_page = 1
if "project_filter" not in st.session_state: st.session_state.project_filter = None
if "show_logs" not in st.session_state: st.session_state.show_logs = False
if "logs_match_sel" not in st.session_state: st.session_state.logs_match_sel = None
if "logs_subdir_sel" not in st.session_state: st.session_state.logs_subdir_sel = None
if "tested_found_paths" not in st.session_state: st.session_state.tested_found_paths = []
if "tested_mapping" not in st.session_state: st.session_state.tested_mapping = {}

# ---------- inject header/footer ----------
render_portal_chrome()

# ---------- SETTINGS PANEL ----------
st.markdown('<div id="settings-panel"></div>', unsafe_allow_html=True)
with st.expander("Settings", expanded=False):
    st.caption("General")
    st.text_input("Root (read-only)", value=ROOT, disabled=True)
    st.caption("About")
    st.write(f"Portal version: **{PORTAL_VERSION}**")

# ---------- body ----------
st.title("Database Management Portal")
st.markdown("Prototype • Step 1–4 layout (Python/Streamlit)")

# ====== ROW 1: three columns (Step 1 | Loaded CSV | Meta+Progress) ======
col1, col2, col3 = st.columns(3, gap="large")

# --- Step 1 (load CSV + project tiles) ---
with col1:
    st.subheader("Step 1 — Load CSV")
    with st.container(border=True):
        st.caption("Choose analysis mode & CSV source")
        mode = st.radio("Mode", ["New Data Analysis", "Load Existing"], index=0, horizontal=True)
        st.session_state.header_lines = st.number_input("Header lines to skip (CSV)", 0, 20, st.session_state.header_lines, 1)

        if mode == "New Data Analysis":
            uploaded = st.file_uploader("Upload CSV", type=["csv"])
            st.caption("…or load from a local path")
            csv_path = st.text_input("CSV path", value=DEFAULT_CSV_PATH)
            if st.button("Load from path"):
                p = expand_path(csv_path)
                if os.path.isfile(p):
                    df_all, rows = read_csv_any(p, st.session_state.header_lines)
                    st.session_state.df_all, st.session_state.rows = df_all, rows
                    st.session_state.loaded_name = p
                    st.session_state.tile_page = 1
                    st.session_state.project_filter = None
                    st.session_state.line_selected = None
                    st.session_state.show_logs = False
                    st.session_state.logs_match_sel = None
                    st.session_state.logs_subdir_sel = None
                    st.session_state.tested_found_paths = []
                    st.session_state.tested_mapping = {}
                else:
                    st.error(f"File not found: {p}")
            if uploaded is not None:
                df_all, rows = read_csv_any(uploaded, st.session_state.header_lines)
                st.session_state.df_all, st.session_state.rows = df_all, rows
                st.session_state.loaded_name = getattr(uploaded, "name", "uploaded.csv")
                st.session_state.tile_page = 1
                st.session_state.project_filter = None
                st.session_state.line_selected = None
                st.session_state.show_logs = False
                st.session_state.logs_match_sel = None
                st.session_state.logs_subdir_sel = None
                st.session_state.tested_found_paths = []
                st.session_state.tested_mapping = {}
        else:
            st.info("Load an existing analysis directory (placeholder).")
            st.text_input("Analysis root", value=DEFAULT_ANALYSED_ROOT)

        if st.session_state.rows:
            st.markdown("---")
            st.caption("Filter by **Project (Column 3)**")
            projects = unique_projects(st.session_state.rows, st.session_state.header_lines)

            is_all = st.session_state.project_filter in (None, "")
            if st.button("All projects", key="proj_all", type=("primary" if is_all else "secondary"), use_container_width=True):
                st.session_state.project_filter = None
                st.session_state.tile_page = 1
                st.session_state.line_selected = None
                st.session_state.show_logs = False
                st.session_state.logs_match_sel = None
                st.session_state.logs_subdir_sel = None
                st.session_state.tested_found_paths = []
                st.session_state.tested_mapping = {}
                st.rerun()

            for i in range(0, len(projects), PROJECT_TILE_COLS):
                row_items = projects[i:i+PROJECT_TILE_COLS]
                cols = st.columns(len(row_items))
                for c, proj in zip(cols, row_items):
                    with c:
                        is_sel = (st.session_state.project_filter == proj)
                        if st.button(proj, key=f"proj_{i}_{proj}", type=("primary" if is_sel else "secondary"), use_container_width=True):
                            st.session_state.project_filter = proj
                            st.session_state.tile_page = 1
                            st.session_state.line_selected = None
                            st.session_state.show_logs = False
                            st.session_state.logs_match_sel = None
                            st.session_state.logs_subdir_sel = None
                            st.session_state.tested_found_paths = []
                            st.session_state.tested_mapping = {}
                            st.rerun()

# --- Middle column: CSV preview + tiles ---
with col2:
    st.subheader("Loaded CSV")
    with st.container(border=True):
        if st.session_state.rows:
            st.write(f"**File:** {st.session_state.loaded_name or '(in-memory)'}")
            preview_df_all = build_preview(st.session_state.rows, st.session_state.header_lines)
            filtered_lines = lines_for_project(st.session_state.rows, st.session_state.header_lines, st.session_state.project_filter)
            if filtered_lines:
                preview_df = preview_df_all[preview_df_all["CSV Line"].isin(filtered_lines)].copy()
                st.caption(f"Showing {len(preview_df)} rows" + (f" for project **{st.session_state.project_filter}**" if st.session_state.project_filter else " (all projects)"))
                st.dataframe(preview_df, use_container_width=True, hide_index=True, height=300)
            else:
                st.warning("No rows match the selected project.")
                preview_df = preview_df_all.iloc[0:0].copy()

            all_lines = preview_df["CSV Line"].tolist()
            total = len(all_lines)
            per_page = TILE_COLS * TILE_ROWS
            total_pages = max(1, (total + per_page - 1) // per_page)
            if st.session_state.tile_page > total_pages:
                st.session_state.tile_page = 1
            page = st.number_input("Page", 1, total_pages, st.session_state.tile_page, 1, help=f"{total} lines • {per_page} per page")
            st.session_state.tile_page = page

            start = (page - 1) * per_page
            end = min(start + per_page, total)
            subset = all_lines[start:end]

            st.caption("Select CSV line to analyze (click a tile)")
            for r in range(TILE_ROWS):
                row_items = subset[r * TILE_COLS : (r + 1) * TILE_COLS]
                if not row_items: break
                cols = st.columns(len(row_items))
                for c, ln in zip(cols, row_items):
                    with c:
                        selected = (st.session_state.get("line_selected") == ln)
                        if st.button(str(ln), key=f"line_tile_{ln}", type=("primary" if selected else "secondary"), use_container_width=True):
                            st.session_state.line_selected = ln
                            st.session_state.show_logs = False
                            st.session_state.logs_match_sel = None
                            st.session_state.logs_subdir_sel = None
                            st.session_state.tested_found_paths = []
                            st.session_state.tested_mapping = {}
                            st.rerun()
        else:
            st.caption("No CSV loaded yet. Load it in Step 1.")

# --- Right column: Meta + Progress / Logs selection ---
with col3:
    st.subheader("Step 2 — Meta Data")
    with st.container(border=True):
        if st.session_state.rows and st.session_state.line_selected:
            ln = st.session_state.line_selected
            row = st.session_state.rows[ln - 1]
            meta = extract_meta(row)
            meta_field("CSV Line", str(ln))
            meta_field("Project Name", meta["Project_Name"])
            meta_field("Subproject Name", meta["Sub_project_Name"])
            meta_field("Dataset Name", meta["Dataset_Name"])
            meta_field("Functional Run Number", meta["Functional_Run"])
            meta_field("Structural Run Number", meta["Structural_Run"])
            meta_field("Structural for Coregistration", meta["Structural_for_CoReg"])
            meta_field("Baseline Duration (min)", meta["Baseline_Duration_min"])
            meta_field("Injection Duration (min)", meta["Injection_Duration_min"])
            meta_field("Injected Left", meta["Injected_Left"])
            meta_field("Injected Right", meta["Injected_Right"])
        else:
            st.caption("Select a CSV line to populate metadata.")

    st.subheader("Step 3 — Progress")
    with st.container(border=True):
        if st.session_state.line_selected:
            ln = st.session_state.line_selected
            row = st.session_state.rows[ln - 1]
            meta = extract_meta(row)

            st.markdown("**Current status:** Ready to run motion correction (Step 1)")
            st.progress(25)
            b1, b2 = st.columns(2)
            b1.button("Start Analysis (placeholder)")

            if b2.button("View Logs"):
                st.session_state.show_logs = True

            if st.session_state.show_logs:
                base = analysed_base_path(ROOT, meta)
                col6 = meta.get("Functional_Run", "")
                matches = find_matches_containing_token(base, col6)

                st.markdown("**Logs — folder selection (Col6)**")
                st.caption(f"Base: `{base}`  •  Contains: `{col6}`")
                if not os.path.isdir(base):
                    st.warning(f"Base path does not exist yet:\n{base}")
                elif not matches:
                    st.info("No matching folders found.")
                else:
                    if st.session_state.logs_match_sel not in matches:
                        st.session_state.logs_match_sel = matches[0] if len(matches) == 1 else None
                    if len(matches) > 1:
                        for i in range(0, len(matches), MATCH_TILE_COLS):
                            row_items = matches[i:i+MATCH_TILE_COLS]
                            cols = st.columns(len(row_items))
                            for c, p in zip(cols, row_items):
                                with c:
                                    is_sel = (st.session_state.logs_match_sel == p)
                                    label = os.path.basename(p)
                                    if st.button(label, key=f"match_tile_{i}_{label}", type=("primary" if is_sel else "secondary"), use_container_width=True):
                                        st.session_state.logs_match_sel = p
                                        st.session_state.logs_subdir_sel = None
                                        st.rerun()
                    sel_match = st.session_state.logs_match_sel
                    if sel_match:
                        st.markdown(f"**Inside:** `{sel_match}`")
                        subdirs = list_subdirs(sel_match)
                        if not subdirs:
                            st.info("No subfolders in this directory.")
                        else:
                            for i in range(0, len(subdirs), MATCH_TILE_COLS):
                                row_items = subdirs[i:i+MATCH_TILE_COLS]
                                cols = st.columns(len(row_items))
                                for c, p in zip(cols, row_items):
                                    with c:
                                        is_sel = (st.session_state.logs_subdir_sel == p)
                                        if st.button(os.path.basename(p), key=f"subdir_{i}_{os.path.basename(p)}", type=("primary" if is_sel else "secondary"), use_container_width=True):
                                            st.session_state.logs_subdir_sel = p
                                            st.rerun()
        else:
            st.write("Idle")

# ====== Step 4 ======
st.subheader("Step 4 — Viewer & Outputs")
with st.container(border=True):
    tabs = st.tabs([
        "Checklist", "QC", "Static Maps", "Logs",
        "Viewer (CPU slices)", "Viewer (WebGL)", "Overlay (CPU)"
    ])

    # Checklist
    with tabs[0]:
        st.write("Expected outputs (placeholders):")
        for idx, pat in enumerate(EXPECTED_FILES, start=1):
            st.write(f"Row {idx}: {pat}")

    with tabs[1]:
        st.caption("(Future) mean/std/tSNR/motion plots")

    with tabs[2]:
        st.caption("(Future) baseline/signal windows & SCM overlays")

    # Logs tab: compute results + store TESTED FOUND PATHS for viewer
    with tabs[3]:
        sel = st.session_state.logs_subdir_sel
        if not sel:
            st.info("Select a folder in **Step 3 → View Logs** to check for files.")
        else:
            ln = st.session_state.line_selected
            row = st.session_state.rows[ln - 1]
            meta = extract_meta(row)

            base_common = analysed_base_path(ROOT, meta)
            col7 = meta.get("Structural_for_CoReg", "")
            col7_matches = find_matches_containing_token(base_common, col7)
            selected_leaf = os.path.basename(sel)

            def col7_base_dir_for_cleaned_anat() -> Optional[str]:
                if not col7_matches: return None
                chosen = col7_matches[0]
                leaf_path = os.path.join(chosen, selected_leaf)
                if os.path.isdir(leaf_path): return leaf_path
                subs = list_subdirs(chosen)
                if subs: return subs[0]
                return chosen

            st.markdown(f"**Checking files relative to selected subfolder:** `{sel}`")
            st.caption(f"Col7 token: `{col7}`  •  Base: `{base_common}`")

            results: List[Dict[str, Any]] = []
            for pat in EXPECTED_FILES:
                if pat == "cleaned_anatomy.nii.gz":
                    special_base = col7_base_dir_for_cleaned_anat()
                    if special_base is None:
                        intended_base = os.path.join(base_common, "(no-Col7-match)")
                        results.append({"pattern": pat, "base_dir": intended_base, "found": False, "matches": []})
                    else:
                        results.append(check_expected_file_at(special_base, pat))
                else:
                    results.append(check_expected_file_at(sel, pat))

            st.markdown("### Exact paths that will be tested")
            for res in results:
                pattern = res["pattern"]; found = res["found"]; base_dir = res["base_dir"]; matches = res.get("matches", [])
                if any(ch in pattern for ch in "*?[]"):
                    if matches:
                        for m in matches: st.markdown(f":green[`{m}`]")
                    else:
                        st.markdown(f":red[`{os.path.join(base_dir, pattern)}`]")
                else:
                    abs_path = os.path.join(base_dir, pattern)
                    if found:
                        st.markdown(f":green[`{abs_path}`]")
                    else:
                        st.markdown(f":red[`{abs_path}`]")

            st.markdown("### Status")
            found_paths: List[str] = []
            mapping: Dict[str, List[str]] = {}
            for res in results:
                pattern = res["pattern"]; found = res["found"]; matches = res.get("matches", [])
                cols = st.columns([2.5, 1.5, 3])
                with cols[0]: st.write(pattern)
                with cols[1]:
                    if found: st.success("File found")
                    else:     st.error("File missing")
                with cols[2]:
                    if any(ch in pattern for ch in "*?[]"):
                        st.caption("Matches:")
                        st.write(", ".join(os.path.basename(p) for p in matches) if matches else "—")
                        if matches:
                            mapping[pattern] = matches
                            found_paths.extend(matches)
                    else:
                        abs_path = os.path.join(res["base_dir"], pattern)
                        mapping[pattern] = [abs_path] if found else []
                        if found: found_paths.append(abs_path)

            st.session_state.tested_found_paths = sorted(set(found_paths))
            st.session_state.tested_mapping = mapping

    # Viewer (CPU slices) — ONLY tested found paths
    with tabs[4]:
        paths = st.session_state.get("tested_found_paths") or []
        if not paths:
            st.info("No found NIfTI paths yet. Go to **Logs** tab, run the checks, and ensure files are found.")
        else:
            chosen = st.selectbox("Choose a tested NIfTI file (found)", paths, index=0, key="cpu_file_sel", format_func=os.path.basename)
            size_mb = _get_file_size_mb(chosen)
            if chosen.lower().endswith(".nii.gz") and size_mb > 200:
                st.warning("Large .nii.gz — CPU slice streaming works, but uncompressed .nii will be faster.")
            render_cpu_slice_viewer(chosen)

    # Viewer (WebGL) — ONLY tested found paths
    with tabs[5]:
        paths = st.session_state.get("tested_found_paths") or []
        if not paths:
            st.info("No found NIfTI paths yet. Go to **Logs** tab and run the checks.")
        else:
            chosen = st.selectbox("Choose a tested NIfTI file (WebGL)", paths, index=0, key="webgl_file_sel", format_func=os.path.basename)
            size_mb = _get_file_size_mb(chosen)
            if size_mb > 200:
                st.info("This file is large; the CPU slice viewer may be smoother.")
            url = serve_local_file_as_static_url(chosen)
            render_nifti_viewer_url(url, name=os.path.basename(chosen), height=640)

    # Overlay (CPU) — multi-slice grid; ONLY tested found paths
    with tabs[6]:
        paths = st.session_state.get("tested_found_paths") or []
        if not paths:
            st.info("No found NIfTI paths yet. Go to **Logs** tab and run the checks.")
        else:
            c1, c2 = st.columns(2)
            with c1:
                base_path = st.selectbox("Base image", paths, index=0, key="overlay_base_sel", format_func=os.path.basename)
            with c2:
                overlay_path = st.selectbox("Overlay image", paths, index=min(1, len(paths)-1), key="overlay_over_sel", format_func=os.path.basename)
            if base_path == overlay_path:
                st.warning("Base and overlay are the same file. You can still proceed if you want.")
            render_overlay_viewer_grid(base_path, overlay_path)
