from ..config import *


# FUNCTION LIST
#
#
def _as_list(x):
    return

def _matches(name: str, include=None, exclude=None) -> bool:
    return

def ensure_run_dir_has_tools(run_dir: Path, azimuth_dir: Path) -> Path:
    return

def is_original_dat(dat_path: Path) -> bool:
    return

def load_signal_dat(path: Path) -> np.ndarray:
    return

def downsample_for_plot(y: np.ndarray, max_points: int) -> tuple[np.ndarray, slice]:
    return

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
    return

def delete_if_exists(path: Path) -> None:
    return

def _format_dt_tag(dt_s: float) -> str:
    return

def _binary_dilate_1d(x: np.ndarray, half_window: int) -> np.ndarray:
    return

def find_coincident_event_segments(square_y: np.ndarray, fs: int, proximity_dt: float, min_segment_dt: float) -> list[tuple[int,int]]:
    return

def sample_event_centers(segments: list[tuple[int,int]], n: int, rng: np.random.Generator) -> list[int]:
    return

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
    return






def parse_radius_name_from_dat(dat_path: Path) -> str:
    # KEEP EVENTSCANNING VERSION
    return

def radius_decimal_from_name(radius_name: str) -> float:
    # RM
    return

def find_run_dirs_with_dats(azimuth_dir: Path) -> list[Path]:
    return

def find_extra_dats_by_dir_datname_uniqueness(azimuth_dir: Path) -> dict[Path, list[Path]]:
    return

def process_run_dir_with_midas(
    run_dir: Path,
    azimuth_dir: Path,
    distance_bias: float,
    dat_paths: list[Path] | None = None,
) -> None:
    return

def compute_sample_time(dat_path: Path, fs: int = SAMPLE_FREQ) -> int:
    # ES
    return

def edit_inp_file(
    inp_path: Path,
    r_over_R: float,
    sample_time: int,
    filename_base: str,
    distance_bias: float = 0.0,
):
    # ES
    return

def run_midas_for_dat(dat_path: Path, apac_exe: Path):
    # def run_midas_apac(cwd: Path, exe_path: Path) -> None:
    # ES
    return

def filename_order_from_inp(inp_path: Path) -> list[str]:
    return

def compile_azimuth_tab(
    azimuth_dir: Path,
    azimuth_name: str,
    executable_name: str,
    extra_run_dirs: list[Path] | None = None,
    radius_order: list[str] | None = None,
):
    return

def excel_range_from_top_left(top_left: str, nrows: int, ncols: int):
    return

def coerce_cell(value: str):
    return

def read_tab_as_2d_list_numeric(tab_path: Path):
    return

def paste_tab_into_sheet(ws, azimuth: str, tab_path: Path, clear_pad=(0, 0)):
    return

def phase1_run_midas_and_build_tabs(distance_bias: float, max_workers: int = 4) -> None:
    return

def phase15_fix_r85_in_azimuth_tabs(make_backups: bool = True, tol: float = 1e-9):
    return

def copy_xlsm_to_modified(suffix):
    return

def find_condition_and_port_for_xlsm(xlsm_path: Path, suffix):
    return

def infer_tab_shape_for_port(port_dir: Path, azimuths, default_rows=21, default_cols=50):
    return

def clear_azimuth_block(ws, azimuth: str, nrows: int, ncols: int):
    return

def find_existing_azimuth_tabs(port_dir: Path, azimuths, id):
    return

def paste_data_into_sheet(ws, azimuth: str, data_2d):
    return

def phase2_edit_xlsm_from_tabs(suffix, id):
    return

def phase25_recalc_saveas_fallback(
    folder: Path,
    pattern="*.xlsm",
    visible=True,
    local_outdir: Path = Path(r"C:\Temp\neup_cachefix_out")
):
    return

