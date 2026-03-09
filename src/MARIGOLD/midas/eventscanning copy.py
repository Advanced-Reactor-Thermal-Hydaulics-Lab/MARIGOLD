from ..config import *



from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess
import random
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt



# Root directory (same meaning as in remaster_v3)
BASE_DIR = Path(r"D:\NEUP")

# Angles: folders whose names end with "deg"
ANGLE_SUFFIX = "deg"

# Subdirectory that contains conductivity data under each ANGLE folder
CONDUCTIVITY_SUBDIR = "Conductivity"

# Name of the input file in each RUN directory
INP_FILENAME = "Input.inp"

# Sampling frequency (Hz)
SAMPLE_FREQ = 50000
import matplotlib.pyplot as plt

# Azimuth folder names (update if your dataset differs)
AZIMUTHS = None  # set to e.g. ["90","67.5","45","22.5","00"] to restrict; None auto-detects

# --- Optional traversal filters (for runtime control) ---
# Use None to disable a filter. Use a string or list of strings.
# Matching rule: case-insensitive substring match.
ANGLE_INCLUDE = None   # e.g. ["90deg", "60deg"] or ["90"]
ANGLE_EXCLUDE = None
PORT_INCLUDE  = None   # e.g. ["P1_"] or ["P2_25D"]
PORT_EXCLUDE  = None
AZIMUTH_INCLUDE = None # e.g. ["90", "45"]
AZIMUTH_EXCLUDE = None
POSITION_INCLUDE = None # e.g. ["r9o", "r85o"]
POSITION_EXCLUDE = None

# --- Selection mode (toggleable) ---
# MODE options:
#   'all'      : process every discovered .DAT (ignores include/exclude filters below)
#   'filtered' : process only .DATs that pass include/exclude filters below
#   'random'   : process a random sample of .DATs (optionally after filtering)
PROCESS_MODE = 'random'

# Random sampling controls (only used when PROCESS_MODE == 'random')
RANDOM_SAMPLE_N = 25        # number of .DAT files to randomly select
APPLY_FILTERS_IN_RANDOM = True  # if True, sample from the FILTERED pool; else sample from ALL discovered .DATs
RANDOM_SAMPLE_SEED = 12345



# Executable to copy into each RUN directory (set to your MIDAS/APAC build)
# Option A (APAC.exe build tree)
APAC_EXE = (
    BASE_DIR
    / "MIDASv1.14_4_21_14"
    / "MIDASv1.14_4_21_14"
    / "Microsoft Visual Studio"
    / "APAC"
    / "x64"
    / "Release"
    / "APAC.exe"
)
# Option B (single executable)
# APAC_EXE = BASE_DIR / "MIDASv1.14_4_21_14" / "MIDASv1.14d.exe"

# Output folder (created next to this notebook by default)
OUTPUT_ROOT = Path.cwd() / f"SignalPlots_2026-01-23_12-00"

# Plot controls (tune as needed)
PLOT_MAX_POINTS = 100_000   # downsample for plotting if file is huge

# --- Plot styling ---
FIGSIZE = (4, 4)  # inches
FIGURE_FONT_FAMILY = "Times New Roman"

# Show point markers on plotted lines (useful for downsampled windows).
PLOT_SHOW_MARKERS = True
PLOT_MARKER = "."
PLOT_MARKERSIZE = 2.0
# Rough upper bound on how many markers appear per trace (matplotlib markevery).
PLOT_MAX_MARKERS = 250

# Apply global font family for all plots in this notebook/script.
plt.rcParams.update({"font.family": FIGURE_FONT_FAMILY})
SAVE_DPI = 200

OUTPUT_ROOT

# Plot styling / channel selection
# Four columns correspond to sensors 1-4. Colors are applied in order.
SENSOR_COLORS = ["red", "blue", "green", "black"]

# Toggle which sensors to plot (sensor1..sensor4). Default: ON, ON, OFF, OFF.
DEFAULT_SENSOR_MASK = (True, True, False, False)

# Y-axis bounds for normalized voltage plots
Y_LIM = (-0.2, 1.2)


# --- Plotting mode options ---
# If True, make "event window" plots around randomly sampled coincident activations (recommended for long records).
SAVE_EVENT_PLOTS = True
# If True, also save the full-length Norm/Square plots (can be slow / not very informative for 30s@50kHz).
SAVE_FULL_LENGTH_PLOTS = False

# Event scan parameters (for SquareSig)
EVENT_N = 8                 # number of random events to sample per DAT
EVENT_WINDOW_DT = 0.002     # half-window [s] around the event center to plot (e.g., 0.05 -> +/-50 ms)
EVENT_PROXIMITY_DT = 0.002  # [s] sensors considered "near each other" if their 1's occur within +/- this window
EVENT_MIN_SEGMENT_DT = 0.0005  # [s] minimum combined-true segment duration to consider an event

# For event plots, default is to show all four sensors regardless of DEFAULT_SENSOR_MASK.
EVENT_PLOT_SENSOR_MASK = (True, True, True, True)

# If you prefer a custom label in filenames, set this; otherwise the executable file stem is used.
EXECUTABLE_LABEL_OVERRIDE = None  # e.g., "MIDASv1.14d" or "APAC"

RANDOM_SEED = 12345  # for reproducible random event sampling

# Optional: restrict port folder name prefixes (None = accept all)
PORT_PREFIXES = None  # e.g. ('P1_','P2_','P3_')



# --- Helpers for parsing / editing / running ---

_SIGNAL_OUT_RE = re.compile(r"^\s*SignalOutput\s*=\s*\d+\s*$", re.IGNORECASE)

def _as_list(x):
    if x is None:
        return None
    if isinstance(x, (list, tuple, set)):
        return [str(v) for v in x]
    return [str(x)]

def _matches(name: str, include=None, exclude=None) -> bool:
    n = str(name).lower()
    inc = _as_list(include)
    exc = _as_list(exclude)
    if inc is not None and len(inc) > 0:
        inc_l = [s.lower() for s in inc]
        if not any(s in n for s in inc_l):
            return False
    if exc is not None and len(exc) > 0:
        exc_l = [s.lower() for s in exc]
        if any(s in n for s in exc_l):
            return False
    return True


def parse_radius_name_from_dat(dat_path: Path) -> str:
    """Example: ...\r9o.DAT -> 'r9o'"""
    return dat_path.stem

def radius_decimal_from_name(radius_name: str) -> float:
    """Convert rXo / rXi / rX into decimal r/R (e.g., r85o -> 0.85)."""
    name = radius_name.lower()
    if not name.startswith("r"):
        raise ValueError(f"Unexpected radius name: {radius_name}")
    core = name[1:]
    if core.endswith(("o","i")):
        digits = core[:-1]
    else:
        digits = core
    if not digits.isdigit():
        raise ValueError(f"Unexpected radius digits in name: {radius_name}")
    if len(digits) == 1:
        return float(f"0.{digits}")
    return float(f"{digits[0]}.{digits[1:]}")

def compute_sample_time(dat_path: Path, fs: int = SAMPLE_FREQ) -> int:
    data = np.loadtxt(dat_path)
    n = data.shape[0] if hasattr(data, "shape") else len(data)
    return int(round(n / fs))

def edit_inp_file_for_single_dat(
    inp_path: Path,
    r_over_R: float,
    sample_time: int,
    filename_base: str,
    signal_output: int = 1,
) -> None:
    """
    Update an existing Input.inp so that MIDAS/APAC processes ONLY one DAT and emits signal outputs.

    Edits (line-start match, ignoring leading whitespace):
      - r/R=
      - measuretime=
      - Filename=
      - SignalOutput=
    """
    lines_out: list[str] = []
    found_signal_output = False

    with inp_path.open("r") as f:
        for line in f:
            stripped = line.lstrip()

            if stripped.startswith("r/R="):
                lines_out.append(f"r/R={r_over_R:.2f}\n")

            elif stripped.lower().startswith("measuretime="):
                lines_out.append(f"measuretime={sample_time}\n")

            elif stripped.lower().startswith("filename="):
                lines_out.append(f"Filename={filename_base}\n")

            elif _SIGNAL_OUT_RE.match(stripped):
                lines_out.append(f"SignalOutput={int(signal_output)}\n")
                found_signal_output = True
            else:
                lines_out.append(line)

    if not found_signal_output:
        # Add it at the end (safe for MIDAS/APAC templates that don't include this key)
        if lines_out and not lines_out[-1].endswith("\n"):
            lines_out[-1] = lines_out[-1] + "\n"
        lines_out.append(f"SignalOutput={int(signal_output)}\n")

    with inp_path.open("w") as f:
        f.writelines(lines_out)

def ensure_run_dir_has_tools(run_dir: Path, azimuth_dir: Path) -> Path:
    """Ensure APAC/MIDAS exe and Input.inp exist in the run_dir. Return the exe path within run_dir."""
    if not APAC_EXE.is_file():
        raise FileNotFoundError(f"APAC executable not found at: {APAC_EXE}")

    target_exe = run_dir / APAC_EXE.name
    if not target_exe.is_file():
        shutil.copy2(APAC_EXE, target_exe)

    inp_path = run_dir / INP_FILENAME
    if not inp_path.is_file():
        template_inp = azimuth_dir / INP_FILENAME
        if template_inp.is_file():
            shutil.copy2(template_inp, inp_path)
        else:
            raise FileNotFoundError(f"No {INP_FILENAME} found in {run_dir} and no template in {azimuth_dir}")

    return target_exe

def run_midas_apac(cwd: Path, exe_path: Path) -> None:
    subprocess.run([str(exe_path)], cwd=cwd, check=True)

def is_original_dat(dat_path: Path) -> bool:
    """Avoid re-processing generated files like *_NormSig.dat etc."""
    n = dat_path.stem.lower()
    return not (n.endswith("_mediansig") or n.endswith("_normsig") or n.endswith("_squaresig"))

def load_signal_dat(path: Path) -> np.ndarray:
    return np.loadtxt(path)

def downsample_for_plot(y: np.ndarray, max_points: int) -> tuple[np.ndarray, slice]:
    n = y.shape[0]
    if n <= max_points:
        return y, slice(None)
    step = int(np.ceil(n / max_points))
    return y[::step, :], slice(0, n, step)

def plot_4ch_signal(
    ts: np.ndarray,
    y: np.ndarray,
    title: str,
    out_png: Path,
    sensor_mask: tuple[bool, bool, bool, bool] = DEFAULT_SENSOR_MASK,
    colors: list[str] = SENSOR_COLORS,
    y_lim: tuple[float, float] = Y_LIM,
    dpi: int = SAVE_DPI,
    figsize: tuple[float, float] = FIGSIZE,
    show_markers: bool = PLOT_SHOW_MARKERS,
    marker: str = PLOT_MARKER,
    markersize: float = PLOT_MARKERSIZE,
    max_markers: int = PLOT_MAX_MARKERS,
) -> None:
    """Plot a 4-column signal array with fixed per-sensor colors and an on/off mask.

    Parameters
    ----------
    ts : (N,) ndarray
        Time vector [s].
    y : (N, 4) ndarray
        Four sensor channels.
    sensor_mask : (4,) tuple[bool]
        Whether to plot each channel (sensors 1..4).
    colors : list[str]
        Matplotlib color names for sensors 1..4.
    y_lim : (min, max)
        Y-axis bounds.
    """
    out_png.parent.mkdir(parents=True, exist_ok=True)

    y = np.asarray(y)
    if y.ndim == 1:
        y = y.reshape(-1, 1)

    if y.shape[1] < 4:
        # Be permissive, but keep legend consistent with available cols
        n = y.shape[1]
    else:
        n = 4

    plt.figure(figsize=figsize, constrained_layout=True)

    labels = ["sensor 1", "sensor 2", "sensor 3", "sensor 4"]
    for j in range(n):
        if j < len(sensor_mask) and sensor_mask[j]:
            col = colors[j] if j < len(colors) else None
            if show_markers:
                # Avoid over-plotting: show at most ~max_markers per trace.
                me = max(1, int(np.ceil(len(ts) / max_markers))) if len(ts) else 1
                plt.plot(ts, y[:, j], linewidth=0.8, color=col, label=labels[j],
                         marker=marker, markersize=markersize, markevery=me)
            else:
                plt.plot(ts, y[:, j], linewidth=0.8, color=col, label=labels[j])

    plt.axhline(0.0, linewidth=0.8)
    plt.title(title)
    plt.xlabel("Time [s]")
    plt.ylabel("Normalized voltage")
    plt.ylim(y_lim[0], y_lim[1])
    # Tight x-limits to the plotted window (no whitespace on either side).
    if len(ts) > 0:
        plt.xlim(float(ts[0]), float(ts[-1]))
        plt.margins(x=0)
    plt.legend(loc="upper right")
    plt.savefig(out_png, dpi=dpi)
    plt.close()

def delete_if_exists(path: Path) -> None:
    try:
        if path.is_file():
            path.unlink()
    except Exception as e:
        print(f"          !! Could not delete {path}: {e}")


def _format_dt_tag(dt_s: float) -> str:
    """Format dt for filenames."""
    if dt_s < 1.0:
        ms = int(round(dt_s * 1000))
        return f"{ms}ms"
    return f"{dt_s:.2f}s".replace(".", "p")

def _binary_dilate_1d(x: np.ndarray, half_window: int) -> np.ndarray:
    """Binary dilation via convolution with a flat window."""
    if half_window <= 0:
        return x.astype(bool)
    k = 2 * half_window + 1
    w = np.ones(k, dtype=int)
    # Using int convolution is fast enough for 1D signals
    return (np.convolve(x.astype(int), w, mode="same") > 0)

def find_coincident_event_segments(square_y: np.ndarray, fs: int, proximity_dt: float, min_segment_dt: float) -> list[tuple[int,int]]:
    """Return [(start_idx, end_idx), ...] where all 4 sensors are 'active near each other'.

    'Near each other' is implemented by dilating each sensor's boolean activation by +/- proximity_dt,
    then AND'ing across sensors.

    Parameters
    ----------
    square_y : (N,4)
        Square signal (values expected ~0/1).
    fs : int
        Sampling frequency [Hz].
    proximity_dt : float
        +/- time window [s] used to declare temporal proximity.
    min_segment_dt : float
        Minimum segment duration [s] for the combined coincidence mask.
    """
    y = np.asarray(square_y)
    if y.ndim != 2 or y.shape[1] < 4:
        raise ValueError("SquareSig must be a 2D array with at least 4 columns.")

    # treat >0.5 as "1"
    b = (y[:, :4] > 0.5)

    half = int(round(proximity_dt * fs))
    b_dil = np.column_stack([_binary_dilate_1d(b[:, j], half) for j in range(4)])
    both = np.all(b_dil, axis=1)

    if not np.any(both):
        return []

    # Find contiguous True runs in `both`
    idx = np.flatnonzero(both)
    # breaks where diff > 1
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = [idx[0]] + [idx[i+1] for i in breaks]
    ends   = [idx[i] for i in breaks] + [idx[-1]]

    min_len = max(1, int(round(min_segment_dt * fs)))
    segs = [(s, e) for s, e in zip(starts, ends) if (e - s + 1) >= min_len]
    return segs

def sample_event_centers(segments: list[tuple[int,int]], n: int, rng: np.random.Generator) -> list[int]:
    """Randomly sample up to n segments; return their center indices."""
    if len(segments) == 0:
        return []
    n = min(n, len(segments))
    picks = rng.choice(len(segments), size=n, replace=False)
    centers = []
    for k in picks:
        s, e = segments[int(k)]
        centers.append((s + e) // 2)
    centers.sort()
    return centers

def plot_event_window_pair(
    ctx: "DatContext",
    position: str,
    exe_label: str,
    center_idx: int,
    dt_s: float,
    norm_y: np.ndarray,
    square_y: np.ndarray,
    fs: int,
    sensor_mask: tuple[bool, bool, bool, bool],
) -> None:
    """Plot and save Norm + Square windows around a center index."""
    w = int(round(dt_s * fs))
    n = norm_y.shape[0]
    i0 = max(0, center_idx - w)
    i1 = min(n, center_idx + w)
    sl = slice(i0, i1)

    ts = (np.arange(i0, i1) / fs).astype(float)
    dt_tag = _format_dt_tag(dt_s)

    # Save into stage folders, with the requested naming scheme
    out_norm = output_png_path(ctx, position=position, stage="Norm", filename=f"{center_idx}_{dt_tag}_{exe_label}.png")
    out_sq   = output_png_path(ctx, position=position, stage="Square", filename=f"{center_idx}_{dt_tag}_{exe_label}.png")

    plot_4ch_signal(ts, norm_y[sl, :4], title=f"{ctx.port} | {ctx.azimuth} | {position} | Norm | idx={center_idx}", out_png=out_norm, sensor_mask=sensor_mask)
    plot_4ch_signal(ts, square_y[sl, :4], title=f"{ctx.port} | {ctx.azimuth} | {position} | Square | idx={center_idx}", out_png=out_sq, sensor_mask=sensor_mask)




# --- Directory traversal (mirrors remaster_v3 structure) ---

@dataclass(frozen=True)
class DatContext:
    angle: str
    port: str
    azimuth: str
    run_dir: Path
    dat_path: Path


def _get_angle_dirs(base_dir: Path) -> list[Path]:
    """Robustly interpret BASE_DIR.
    Accepts:
      - a NEUP root that contains <ANGLE>deg folders
      - a specific <ANGLE>deg folder
      - a Conductivity folder (contains PORT folders)
    """
    if base_dir.is_dir():
        # If user points directly at an ANGLE folder
        if base_dir.name.endswith(ANGLE_SUFFIX) and (base_dir / CONDUCTIVITY_SUBDIR).is_dir():
            return [base_dir]
        # If user points directly at the Conductivity folder
        if base_dir.name.lower() == CONDUCTIVITY_SUBDIR.lower():
            return [base_dir.parent]
        # Otherwise scan children for ANGLE folders
        angle_dirs = [p for p in base_dir.iterdir() if p.is_dir() and p.name.endswith(ANGLE_SUFFIX)]
        return sorted(angle_dirs, key=lambda p: p.name)
    return []
def iter_dat_contexts(
    base_dir: Path = BASE_DIR,
    angle_include=ANGLE_INCLUDE,
    angle_exclude=ANGLE_EXCLUDE,
    port_include=PORT_INCLUDE,
    port_exclude=PORT_EXCLUDE,
    azimuth_include=AZIMUTH_INCLUDE,
    azimuth_exclude=AZIMUTH_EXCLUDE,
    position_include=POSITION_INCLUDE,
    position_exclude=POSITION_EXCLUDE,
):
    """Yield DatContext items for all original DAT files under the NEUP Conductivity tree."""
    for angle_dir in _get_angle_dirs(base_dir):
        if not angle_dir.is_dir():
            continue
        if not angle_dir.name.endswith(ANGLE_SUFFIX):
            continue
        if not _matches(angle_dir.name, angle_include, angle_exclude):
            continue

        conductivity_dir = angle_dir / CONDUCTIVITY_SUBDIR
        if not conductivity_dir.is_dir():
            continue

        for port_dir in conductivity_dir.iterdir():
            if not port_dir.is_dir():
                continue
            # Port directory filter (optional)
            if PORT_PREFIXES is not None:
                if not port_dir.name.startswith(tuple(PORT_PREFIXES)):
                    continue
            if not _matches(port_dir.name, port_include, port_exclude):
                continue

            # Determine azimuth directories
            if AZIMUTHS is None:
                az_dirs = [p for p in port_dir.iterdir() if p.is_dir()]
                az_dirs = sorted(az_dirs, key=lambda p: p.name)
            else:
                az_dirs = [port_dir / az for az in AZIMUTHS]

            for azimuth_dir in az_dirs:
                azimuth = azimuth_dir.name
                if not azimuth_dir.is_dir():
                    continue
                if not _matches(azimuth, azimuth_include, azimuth_exclude):
                    continue

                run_dirs = [p for p in azimuth_dir.iterdir() if p.is_dir() and any(list(p.glob("*.DAT")) + list(p.glob("*.dat")))]
                run_dirs = sorted(run_dirs, key=lambda p: p.name)

                # Root-level DATs (some datasets place r*o.dat directly under AZIMUTH)
                for dat_path in sorted(list(azimuth_dir.glob("*.DAT")) + list(azimuth_dir.glob("*.dat")), key=lambda p: p.name):
                    if is_original_dat(dat_path) and _matches(dat_path.stem, position_include, position_exclude):
                        # Treat the AZIMUTH dir itself as a 'run_dir' for execution purposes
                        yield DatContext(angle_dir.name, port_dir.name, azimuth, azimuth_dir, dat_path)

                # Run-dir DATs
                for run_dir in run_dirs:
                    for dat_path in sorted(run_dir.glob("*.DAT"), key=lambda p: p.name):
                        if is_original_dat(dat_path) and _matches(dat_path.stem, position_include, position_exclude):
                            yield DatContext(angle_dir.name, port_dir.name, azimuth, run_dir, dat_path)

                # Extra DATs deeper in tree whose stems are not present at root (same rule as remaster_v3)
                root_stems = {p.stem for p in list(azimuth_dir.glob("*.DAT")) + list(azimuth_dir.glob("*.dat"))}
                for dat_path in azimuth_dir.rglob("*.DAT"):
                    if dat_path.parent == azimuth_dir:
                        continue
                    if dat_path.stem in root_stems:
                        continue
                    if is_original_dat(dat_path) and _matches(dat_path.stem, position_include, position_exclude):
                        yield DatContext(angle_dir.name, port_dir.name, azimuth, dat_path.parent, dat_path)


def build_context_from_dat(dat_path: Path, base_dir: Path = BASE_DIR) -> DatContext:
    """Build a DatContext from a single DAT path (best-effort parsing from path structure).

    Expected pattern (typical):
        .../<ANGLE*deg>/Conductivity/<PORT>/<AZIMUTH>/.../<RUN_DIR>/<FILE>.DAT
    """
    dat_path = Path(dat_path)
    if not dat_path.is_file():
        raise FileNotFoundError(dat_path)

    # Find ANGLE folder (closest ancestor ending with ANGLE_SUFFIX)
    angle_dir = None
    for p in [dat_path.parent, *dat_path.parents]:
        if p.name.endswith(ANGLE_SUFFIX):
            angle_dir = p
            break
    if angle_dir is None:
        raise ValueError(f"Could not infer ANGLE folder for: {dat_path}")

    # Find Conductivity, then infer PORT and AZIMUTH
    conductivity_dir = None
    port = None
    azimuth = None
    parts = list(dat_path.parts)
    try:
        idx = parts.index(CONDUCTIVITY_SUBDIR)
        port = parts[idx + 1]
        azimuth = parts[idx + 2]
        conductivity_dir = Path(*parts[: idx + 1])
    except Exception:
        # Fallback: walk down from angle_dir
        c = angle_dir / CONDUCTIVITY_SUBDIR
        if c.exists():
            conductivity_dir = c

    if port is None or azimuth is None:
        # Fallback: locate conductivity_dir in parents, then take next two segments
        for i, p in enumerate(dat_path.parents):
            if p.name == CONDUCTIVITY_SUBDIR:
                # p is .../Conductivity, so dat_path relative to p should begin with PORT/AZIMUTH/...
                rel = dat_path.relative_to(p)
                rel_parts = rel.parts
                if len(rel_parts) >= 3:
                    port = rel_parts[0]
                    azimuth = rel_parts[1]
                break

    if port is None or azimuth is None:
        raise ValueError(f"Could not infer PORT/AZIMUTH for: {dat_path}")

    # Choose run_dir: directory that contains the DAT (execution happens here)
    run_dir = dat_path.parent

    return DatContext(angle=angle_dir.name, port=port, azimuth=str(azimuth), run_dir=run_dir, dat_path=dat_path)


def output_png_path(ctx: DatContext, position: str, stage: str, filename: str | None = None) -> Path:
    """Return output path for a PNG.

    Tree: OUTPUT_ROOT/<ANGLE>/<PORT>/<AZIMUTH>/<POSITION>/<STAGE>/<filename>
    If filename is None, defaults to '<POSITION>_<STAGE>.png'.
    """
    if filename is None:
        filename = f"{position}_{stage}.png"
    return OUTPUT_ROOT / ctx.angle / ctx.port / ctx.azimuth / position / stage / filename

def process_one_dat(ctx: DatContext,  sensor_mask: tuple[bool, bool, bool, bool] = DEFAULT_SENSOR_MASK) -> None:
    run_dir = ctx.run_dir
    # Best-effort: identify the AZIMUTH directory (used only for locating a template Input.inp)
    # Expected structure: .../Conductivity/<PORT>/<AZIMUTH>/.../<RUN_DIR>/<file>.DAT
    azimuth_dir = ctx.run_dir
    p = ctx.run_dir
    while p != p.parent:
        # If this directory is directly under the PORT directory, treat it as AZIMUTH
        if p.parent.name == ctx.port:
            azimuth_dir = p
            break
        p = p.parent

    # Ensure tools exist in the execution directory
    exe_in_run = ensure_run_dir_has_tools(run_dir, azimuth_dir)
    inp_path = run_dir / INP_FILENAME

    dat_path = ctx.dat_path
    position = dat_path.stem

    # Edit Input.inp to process this single DAT and output signals
    r_over_R = radius_decimal_from_name(parse_radius_name_from_dat(dat_path))
    sample_time = compute_sample_time(dat_path)

    edit_inp_file_for_single_dat(
        inp_path=inp_path,
        r_over_R=r_over_R,
        sample_time=sample_time,
        filename_base=position,
        
        signal_output=1,
    )

    # Run MIDAS/APAC
    try:
        run_midas_apac(cwd=run_dir, exe_path=exe_in_run)

    except subprocess.CalledProcessError as e:
        print(f"[APAC ERROR] Exit code: {e.returncode}")
        print(f"Command: {e.cmd}")
        pass

    # Generated files (by convention from signalquality_master.ipynb)
    norm_path = run_dir / f"{position}_NormSig.dat"
    sqr_path  = run_dir / f"{position}_SquareSig.dat"
    mdn_path  = run_dir / f"{position}_MedianSig.dat"
    sqpr_path = run_dir / f"{position}_SquarePairSig.dat"

    # Executable label for filenames
    exe_label = EXECUTABLE_LABEL_OVERRIDE or Path(exe_in_run).stem

    # Load full signals once (needed for event scanning)
    norm_y = load_signal_dat(norm_path) if norm_path.is_file() else None
    sqr_y  = load_signal_dat(sqr_path)  if sqr_path.is_file()  else None

    # --- Full-length plots (optional) ---
    if SAVE_FULL_LENGTH_PLOTS and norm_y is not None:
        y_ds, s = downsample_for_plot(norm_y, PLOT_MAX_POINTS)
        ts = (np.arange(norm_y.shape[0])[s] / SAMPLE_FREQ).astype(float)
        out_png = output_png_path(ctx, position=position, stage="Norm")
        plot_4ch_signal(ts, y_ds, title=f"{ctx.angle} | {ctx.port} | az={ctx.azimuth} | {position} | Norm", out_png=out_png, sensor_mask=sensor_mask)

    if SAVE_FULL_LENGTH_PLOTS and sqr_y is not None:
        y_ds, s = downsample_for_plot(sqr_y, PLOT_MAX_POINTS)
        ts = (np.arange(sqr_y.shape[0])[s] / SAMPLE_FREQ).astype(float)
        out_png = output_png_path(ctx, position=position, stage="Square")
        plot_4ch_signal(ts, y_ds, title=f"{ctx.angle} | {ctx.port} | az={ctx.azimuth} | {position} | Square", out_png=out_png, sensor_mask=sensor_mask)

    # --- Event-window plots (recommended for long records) ---
    if SAVE_EVENT_PLOTS and (norm_y is not None) and (sqr_y is not None):
        try:
            segs = find_coincident_event_segments(
                square_y=sqr_y,
                fs=SAMPLE_FREQ,
                proximity_dt=EVENT_PROXIMITY_DT,
                min_segment_dt=EVENT_MIN_SEGMENT_DT,
            )
            rng = np.random.default_rng(RANDOM_SEED)
            centers = sample_event_centers(segs, n=EVENT_N, rng=rng)

            # If no events found, fall back to a single centered window at mid-record (optional)
            if len(centers) == 0 and norm_y.shape[0] > 0:
                centers = [norm_y.shape[0] // 2]

            for cidx in centers:
                plot_event_window_pair(
                    ctx=ctx,
                    position=position,
                    exe_label=exe_label,
                    center_idx=int(cidx),
                    dt_s=float(EVENT_WINDOW_DT),
                    norm_y=norm_y,
                    square_y=sqr_y,
                    fs=SAMPLE_FREQ,
                    sensor_mask=EVENT_PLOT_SENSOR_MASK,
                )
        except Exception as e:
            print(f"          !! Event plotting failed for {position}: {e}")

    # Delete generated signal outputs immediately (keep original DAT and any TAB)
    delete_if_exists(mdn_path)
    delete_if_exists(norm_path)
    delete_if_exists(sqr_path)
    delete_if_exists(sqpr_path)


def process_single_dat(dat_path: str | Path,  sensor_mask: tuple[bool, bool, bool, bool] = DEFAULT_SENSOR_MASK) -> None:
    """Process a single original DAT file: run MIDAS/APAC with SignalOutput=1, plot Norm/Square, cleanup outputs."""
    dat_path = Path(dat_path)
    if not is_original_dat(dat_path):
        raise ValueError(f"Refusing to process non-original DAT: {dat_path.name}")
    ctx = build_context_from_dat(dat_path, base_dir=BASE_DIR)
    process_one_dat(ctx,  sensor_mask=sensor_mask)


def diagnose_tree(max_show: int = 10) -> None:
    """Print what the traversal sees under BASE_DIR."""
    angle_dirs = _get_angle_dirs(BASE_DIR)
    print(f"BASE_DIR={BASE_DIR}")
    print(f"Found {len(angle_dirs)} angle dir(s):")
    for p in angle_dirs[:max_show]:
        print("  -", p)
    dats = []
    for i, ctx in enumerate(iter_dat_contexts(BASE_DIR)):
        dats.append(ctx.dat_path)
        if len(dats) >= max_show:
            break
    print(f"First {len(dats)} DAT(s) discovered:")
    for d in dats:
        print("  -", d)

def run_all(
    sensor_mask: tuple[bool, bool, bool, bool] = DEFAULT_SENSOR_MASK,
    limit: int | None = None,
    mode: str = PROCESS_MODE,
    random_n: int | None = RANDOM_SAMPLE_N,
    random_seed: int = RANDOM_SAMPLE_SEED,
    apply_filters_in_random: bool = APPLY_FILTERS_IN_RANDOM,
    angle_include=ANGLE_INCLUDE,
    angle_exclude=ANGLE_EXCLUDE,
    port_include=PORT_INCLUDE,
    port_exclude=PORT_EXCLUDE,
    azimuth_include=AZIMUTH_INCLUDE,
    azimuth_exclude=AZIMUTH_EXCLUDE,
    position_include=POSITION_INCLUDE,
    position_exclude=POSITION_EXCLUDE,
) -> None:
    """Process DAT files according to a toggleable selection mode.

    mode:
      - 'all'      : ignore include/exclude filters and process every discovered DAT
      - 'filtered' : process only DATs that pass include/exclude filters
      - 'random'   : randomly sample DATs either from the filtered pool (apply_filters_in_random=True)
                    or from all discovered DATs (apply_filters_in_random=False)
    """
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    mode_norm = str(mode).strip().lower()
    if mode_norm not in {"all", "filtered", "random"}:
        raise ValueError(f"Unknown mode={mode!r}. Use 'all', 'filtered', or 'random'.")

    # --- Discover contexts ---
    if mode_norm == "all":
        contexts = list(iter_dat_contexts(BASE_DIR))

    elif mode_norm == "filtered":
        contexts = list(
            iter_dat_contexts(
                BASE_DIR,
                angle_include=angle_include,
                angle_exclude=angle_exclude,
                port_include=port_include,
                port_exclude=port_exclude,
                azimuth_include=azimuth_include,
                azimuth_exclude=azimuth_exclude,
                position_include=position_include,
                position_exclude=position_exclude,
            )
        )

    else:  # random
        # Always discover broadly first
        all_contexts = list(iter_dat_contexts(BASE_DIR))
        if apply_filters_in_random:
            pool = list(
                iter_dat_contexts(
                    BASE_DIR,
                    angle_include=angle_include,
                    angle_exclude=angle_exclude,
                    port_include=port_include,
                    port_exclude=port_exclude,
                    azimuth_include=azimuth_include,
                    azimuth_exclude=azimuth_exclude,
                    position_include=position_include,
                    position_exclude=position_exclude,
                )
            )
        else:
            pool = all_contexts

        if not pool:
            contexts = []
        else:
            rng = random.Random(int(random_seed))
            n = len(pool) if random_n is None else int(random_n)
            if n <= 0:
                contexts = []
            else:
                n = min(n, len(pool))
                contexts = rng.sample(pool, n)

    # --- Apply limit (post-selection) ---
    if limit is not None and limit > 0:
        contexts = contexts[: int(limit)]

    if not contexts:
        print(
            "WARNING: No .DAT files were selected.\n"
            "  - Verify BASE_DIR points to the tree that contains your <ANGLE> folders\n"
            "  - If using filters, relax ANGLE/PORT/AZIMUTH/POSITION include/exclude\n"
            "  - If using random mode, ensure RANDOM_SAMPLE_N > 0 and that the pool is non-empty"
        )
        print(f"Done. PNGs saved under: {OUTPUT_ROOT}")
        return

    # --- Execute ---
    for k, ctx in enumerate(contexts, 1):
        try:
            print(f"[{k}/{len(contexts)}] {ctx.angle} / {ctx.port} / az={ctx.azimuth} :: {ctx.dat_path}")
            process_one_dat(ctx, sensor_mask=sensor_mask)
        except Exception as e:
            print(f"      !! ERROR: {e}")

    print(f"Done. PNGs saved under: {OUTPUT_ROOT}")
    return