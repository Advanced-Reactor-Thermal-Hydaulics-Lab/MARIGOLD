from ..config import *

_DEFAULT_COLORS = [
    "#E84393",
    "#6BBBE8",
    "#7ED957",
    "#8E44AD",
    "#F39C12",
]
_DEFAULT_MARKERS = ["o"]
_DEFAULT_LINESTYLES = ["--", "-"]

def _safe_get(midas_dict, key, default=np.nan):
    try:
        return midas_dict[key]
    except Exception:
        return default

def plot_xy_profiles(
    cond,
    x_param: str,
    y_param: str,
    phis_to_plot=(90, 67.5, 45, 22.5, 0, 270, 247.5, 225, 202.5, 180),
    *,
    ax=None,
    show_markers: bool = True,
    show_lines: bool = True,
    x_label: str | None = None,
    y_label: str | None = None,
    set_xlim=None,                  # tuple (xmin, xmax) or None
    set_ylim=None,                  # tuple (ymin, ymax) or None
    percent_error: float = 0.0,     # e.g. 0.05 => 5% error bars; 0 => none
    legend_phi: bool = True,
    title: bool = True,
    fig_size=(4, 4),
    fs: int = 10,
    save_dir: str | None = None,
    show: bool = True,
    close: bool = True,
    db_label: str | None = None,    # used for overlay legend in multi-db plotting
    style=None,                     # dict to override color/marker/linestyle per call if desired
):
    if ax is None:
        fig, ax = plt.subplots(figsize=fig_size, dpi=300, layout="compressed")
    else:
        fig = ax.figure

    plt.rcParams.update({"font.size": fs})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"
    ax.tick_params(direction="in", which="both")
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")

    color_cyc = cycle(_DEFAULT_COLORS)
    marker_cyc = cycle(_DEFAULT_MARKERS)
    ls_cyc = cycle(_DEFAULT_LINESTYLES)

    # Map base_phi -> (color, marker, linestyle) so complements reuse identical styling
    style_by_base_phi: dict[float, tuple[str, str, str]] = {}

    def _base_phi(phi_val: float) -> float:
        # Treat complements as the same "family"
        return phi_val if phi_val <= 180 else phi_val - 180

    def _resolve_style_for_phi(phi_val: float):
        base = _base_phi(phi_val)

        if base not in style_by_base_phi:
            color = next(color_cyc)
            marker = next(marker_cyc)
            ls = next(ls_cyc)

            # DB-level overrides: keep constant across all phis in this call if provided
            if style:
                color = style.get("color", color)
                marker = style.get("marker", marker)
                ls = style.get("linestyle", ls)

            # Respect toggles
            if not show_lines:
                ls = "None"
            if not show_markers:
                marker = "None"

            style_by_base_phi[base] = (color, marker, ls)

        return base, style_by_base_phi[base]

    # Collect and plot each azimuth
    plotted_any = False
    for phi in phis_to_plot:
        if phi not in cond.data:
            continue

        r_dict = cond.data[phi]
        rs_sorted = sorted(r_dict.keys())

        xs, ys = [], []
        for r in rs_sorted:
            if abs(r - 1.0) < 1e-6:     # r* = 1.0 dummy
                continue

            midas_out = r_dict[r]
            x = _safe_get(midas_out, x_param, np.nan)
            y = _safe_get(midas_out, y_param, np.nan)

            # Omit y == 0 points and break line
            if y == 0.0:
                xs.append(np.nan)
                ys.append(np.nan)
            else:
                xs.append(x)
                ys.append(y)

        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)

        base, (color, marker, ls) = _resolve_style_for_phi(phi)

        # Label only the base angle; omit complements from legend
        phi_label = f"{base}°" if (legend_phi and phi == base) else None

        if percent_error and percent_error > 0:
            xerr = np.abs(xs) * percent_error
            yerr = np.abs(ys) * percent_error
            ax.errorbar(
                xs, ys,
                yerr=yerr, #xerr=xerr,
                capsize=3, ecolor="black",
                color=color, marker=marker, linestyle=ls,
                markeredgecolor="black", markeredgewidth=0.8,
                label=phi_label,
            )
        else:
            ax.plot(
                xs, ys,
                color=color, marker=marker, linestyle=ls,
                markeredgecolor="black", markeredgewidth=0.8,
                label=phi_label,
            )

        plotted_any = True

    if x_label is None:
        x_label = x_param
    if y_label is None:
        y_label = y_param
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    if set_xlim is not None:
        ax.set_xlim(*set_xlim)
    if set_ylim is not None:
        ax.set_ylim(*set_ylim)

    if title:
        ax.set_title(
            rf"$\theta = {cond.theta}^\circ,\; j_f = {cond.jf}\,[m/s],\; j_g = {cond.jgref}\,[m/s]$"
            + f", {cond.port}"
        )

    if legend_phi:
        ax.legend(loc="best", edgecolor="white")
    
    # Save/show
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        fname = f"{x_param}_vs_{y_param}_theta{cond.theta}_jf{cond.jf}_jg{cond.jgref}_{cond.port}.png"
        fig.savefig(os.path.join(save_dir, fname))
        plt.close(fig)
    else:
        if show:
            plt.show()
        elif close:
            plt.close(fig)

    return fig, ax, plotted_any

def _cond_match_key(cond):
    return (cond.theta, cond.jf, cond.jgref, cond.port)

def plot_two_conditions_xy_profiles(
    cond_a,
    cond_b,
    label_a: str,
    label_b: str,
    x_param: str,
    y_param: str,
    phis_to_plot=(90, 270),
    *,
    show_markers: bool = True,
    show_lines: bool = True,
    x_label: str | None = None,
    y_label: str | None = None,
    set_xlim=None,
    set_ylim=None,
    percent_error: float = 0.0,
    fig_size=(4, 4),
    fs: int = 10,
    show: bool = True,
    save_path: str | None = None,
    style_a: dict | None = None,
    style_b: dict | None = None,
):
    """
    """

    plt.rcParams.update({"font.size": fs})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    fig, ax = plt.subplots(figsize=fig_size, dpi=300, layout="compressed")
    ax.tick_params(direction="in", which="both")
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    
    # Auto styles if none provided
    if style_a is None and style_b is None:
        color_cyc = cycle(_DEFAULT_COLORS)
        marker_cyc = cycle(_DEFAULT_MARKERS)
        ls_cyc = cycle(_DEFAULT_LINESTYLES)

        style_a = {
            "color": next(color_cyc),
            "marker": next(marker_cyc),
            "linestyle": next(ls_cyc),
        }
        style_b = {
            "color": next(color_cyc),
            "marker": next(marker_cyc),
            "linestyle": next(ls_cyc),
        }

    elif style_a is None:
        color_cyc = cycle(_DEFAULT_COLORS)
        marker_cyc = cycle(_DEFAULT_MARKERS)
        ls_cyc = cycle(_DEFAULT_LINESTYLES)

        style_a = {
            "color": next(color_cyc),
            "marker": next(marker_cyc),
            "linestyle": next(ls_cyc),
        }

    elif style_b is None:
        color_cyc = cycle(_DEFAULT_COLORS)
        marker_cyc = cycle(_DEFAULT_MARKERS)
        ls_cyc = cycle(_DEFAULT_LINESTYLES)

        style_b = {
            "color": next(color_cyc),
            "marker": next(marker_cyc),
            "linestyle": next(ls_cyc),
        }

    fig, ax, plotted_a = plot_xy_profiles(
        cond_a,
        x_param, y_param,
        phis_to_plot=phis_to_plot,
        ax=ax,
        show_markers=show_markers,
        show_lines=show_lines,
        x_label=x_label,
        y_label=y_label,
        set_xlim=set_xlim,
        set_ylim=set_ylim,
        percent_error=percent_error,
        legend_phi=True,
        title=True,
        fig_size=fig_size,
        fs=fs,
        db_label=label_a,
        style=style_a,
        show=False,
        close=False,
        save_dir=None,
    )

    fig, ax, plotted_b = plot_xy_profiles(
        cond_b,
        x_param, y_param,
        phis_to_plot=phis_to_plot,
        ax=ax,
        show_markers=show_markers,
        show_lines=show_lines,
        x_label=x_label,
        y_label=y_label,
        set_xlim=set_xlim,
        set_ylim=set_ylim,
        percent_error=percent_error,
        legend_phi=False,
        title=True,
        fig_size=fig_size,
        fs=fs,
        db_label=label_b,
        style=style_b,
        show=False,
        close=False,
        save_dir=None,
    )

    bold_title = FontProperties(weight="bold")

    db_proxy_handles = [
    Line2D(
        [0], [0],
        color=(style_a or {}).get("color", "black"),
        linestyle=(style_a or {}).get("linestyle", "-"),
        marker=(style_a or {}).get("marker", "None"),
        linewidth=1.5,
        markeredgecolor="black",
        markeredgewidth=0.8,
        ),
        Line2D(
            [0], [0],
            color=(style_b or {}).get("color", "black"),
            linestyle=(style_b or {}).get("linestyle", "-"),
            marker=(style_b or {}).get("marker", "None"),
            linewidth=1.5,
            markeredgecolor="black",
            markeredgewidth=0.8,
        ),
    ]
    db_proxy_labels = [label_a, label_b]

    ax.legend(
        db_proxy_handles,
        db_proxy_labels,
        loc="lower left",
        title="Database",
        title_fontproperties=bold_title,
        edgecolor="white",
    )

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path)
        plt.close(fig)
        return None, None

    if show:
        plt.show()

    return fig, ax

def plot_matched_xy_profiles_across_databases(
    databases,
    db_labels,
    x_param: str,
    y_param: str,
    phis_to_plot=(90, 270),
    *,
    show_markers: bool = True,
    show_lines: bool = True,
    x_label: str | None = None,
    y_label: str | None = None,
    set_xlim=None,
    set_ylim=None,
    percent_error: float = 0.0,
    fig_size=(4, 4),
    fs: int = 10,
    save_dir: str | None = None,
    show: bool = True,
    db_styles: list[dict] | None = None,
):
    """
    """

    if len(databases) != len(db_labels):
        raise ValueError("databases and db_labels must have the same length.")
    if len(databases) != 2:
        raise ValueError("This function expects exactly two databases.")

    db1, db2 = databases
    label1, label2 = db_labels

    # Build indexes
    idx1 = {_cond_match_key(cond): cond for cond in db1}
    idx2 = {_cond_match_key(cond): cond for cond in db2}

    matched_keys = sorted(set(idx1.keys()) & set(idx2.keys()))
    if not matched_keys:
        return []
    
    # Styles
    if db_styles is None:
        color_cyc = cycle(_DEFAULT_COLORS)
        marker_cyc = cycle(_DEFAULT_MARKERS)
        ls_cyc = cycle(_DEFAULT_LINESTYLES)
        db_styles = [
            {"color": next(color_cyc), "marker": next(marker_cyc), "linestyle": next(ls_cyc)},
            {"color": next(color_cyc), "marker": next(marker_cyc), "linestyle": next(ls_cyc)},
        ]
    else:
        if len(db_styles) < 2:
            raise ValueError("db_styles must have at least two entries for two databases.")

    outputs = []

    for key in matched_keys:
        cond_a = idx1[key]
        cond_b = idx2[key]

        save_path = None
        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)
            theta, jf, jgref, port = key
            fname = f"{x_param}_vs_{y_param}_theta{theta}_jf{jf}_jg{jgref}_{port}.png"
            save_path = os.path.join(save_dir, fname)

        fig, ax = plot_two_conditions_xy_profiles(
            cond_a,
            cond_b,
            label1,
            label2,
            x_param,
            y_param,
            phis_to_plot=phis_to_plot,
            show_markers=show_markers,
            show_lines=show_lines,
            x_label=x_label,
            y_label=y_label,
            set_xlim=set_xlim,
            set_ylim=set_ylim,
            percent_error=percent_error,
            fig_size=fig_size,
            fs=fs,
            show=show if save_path is None else False,
            save_path=save_path,
            style_a=db_styles[0],
            style_b=db_styles[1],
        )

        # If strict plotting failed, fig/ax will be None and nothing saved/shown
        if save_path is not None:
            outputs.append(save_path)
        else:
            outputs.append((fig, ax))

    return outputs


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import os

def plot_contour2(
    cond,
    param: str,
    save_dir='.',
    show=True,
    set_max=None,
    set_min=None,
    fig_size=4,
    colorbar_label=None,
    suppress_colorbar=False,
    rot_angle=0,
    ngridr=50,
    ngridphi=50,
    colormap='plasma',
    num_levels=0,
    level_step=0.01,
    title=False,
    title_str='',
    extra_save_text='',
    annotate_h=False,
    cartesian=False,
    h_star_kwargs={'method': 'max_dsm', 'min_void': '0.05'},
    plot_measured_points=False,
    font_size=12,
    # NEW:
    ax=None,
    fig=None,
    add_colorbar=True,   # grid use-case: usually False, do one shared bar outside
):
    """
    Same as your function, but now can draw on a provided `ax`.
    Backward compatible: if ax is None, it creates its own fig/ax as before.
    Returns: ax, mpbl (QuadContourSet), and (optionally) cbar
    """
    created_fig = False
    cbar = None

    if ax is None:
        created_fig = True
        if cartesian:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300)
        else:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300,
                                   subplot_kw=dict(projection='polar'))
    else:
        if fig is None:
            fig = ax.figure

    plt.rcParams.update({'font.size': font_size})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    rs, phis, vals = [], [], []

    for phi_angle, r_dict in cond.data.items():
        for r, midas_output in r_dict.items():
            if r >= 0:
                rs.append(r)
                phis.append(phi_angle)
                try:
                    vals.append(midas_output[param])
                except Exception:
                    vals.append(0.0)
                    # keep your print if you want:
                    # print(f"Could not find {param} for φ = {phi_angle}, r = {r}. Substituting 0")

    rs = np.asarray(rs)
    phis = (np.asarray(phis) + rot_angle) * np.pi / 180
    vals = np.asarray(vals)

    # Periodic seam fix (highly recommended): duplicate points at ±2π
    phis0, rs0, vals0 = phis, rs, vals
    phis = np.concatenate([phis0, phis0 + 2*np.pi, phis0 - 2*np.pi])
    rs   = np.concatenate([rs0,   rs0,            rs0])
    vals = np.concatenate([vals0, vals0,          vals0])

    ri = np.linspace(0, 1, ngridr)
    rot = rot_angle * np.pi / 180
    phii = np.linspace(rot, 2*np.pi + rot, ngridphi, endpoint=False)

    triang = tri.Triangulation(phis, rs)
    interpolator = tri.LinearTriInterpolator(triang, vals)

    PHII, RI = np.meshgrid(phii, ri)
    XI = RI * np.cos(PHII)
    YI = RI * np.sin(PHII)

    parami = interpolator(PHII, RI)
    # handle masked arrays cleanly
    if np.ma.isMaskedArray(parami):
        parami_f = parami.filled(np.nan)
    else:
        parami_f = np.asarray(parami)

    if set_min is None:
        set_min = np.nanmin(parami_f)
    if set_max is None:
        pmax = np.nanmax(parami_f)
        set_max = pmax + 0.1 * pmax

    extend_min = set_min > np.nanmin(parami_f)
    extend_max = set_max < np.nanmax(parami_f)
    if extend_min and extend_max:
        extend_opt = 'both'
    elif extend_min:
        extend_opt = 'min'
    elif extend_max:
        extend_opt = 'max'
    else:
        extend_opt = 'neither'

    if abs(set_max - set_min) < level_step:
        level_step = 0.01 * (set_max - set_min)

    if num_levels:
        lvs = np.linspace(set_min, set_max, num_levels)
    else:
        if abs(set_max) < 1e-8:
            lvs = np.arange(set_max, abs(set_min) + 1e-8, level_step)
            lvs = np.flip(lvs) * -1
        else:
            lvs = np.arange(set_min, set_max + 1e-8, level_step)

    if cartesian:
        mpbl = ax.contourf(XI, YI, parami_f, levels=lvs,
                           vmin=set_min, vmax=set_max, cmap=colormap)

        x = np.linspace(-1, 1, 200)
        ax.plot(x, np.sqrt(1 - x**2), color='black', linewidth=1)
        ax.plot(x, -np.sqrt(1 - x**2), color='black', linewidth=1)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel(r'$x/R$ [-]')
        ax.set_ylabel(r'$y/R$ [-]')
        ax.grid(False)
    else:
        mpbl = ax.contourf(PHII, RI, parami_f, levels=lvs,
                           vmin=set_min, vmax=set_max, cmap=colormap, extend=extend_opt)
        ax.grid(False)
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    if annotate_h:
        hstar = cond.find_hstar_pos(**h_star_kwargs)
        # your original annotate coordinates assumed cartesian-ish scaling
        ax.annotate('h', (-1, 1-hstar), (1, 1-hstar),
                    arrowprops=dict(color='r', width=3, headwidth=3),
                    color='r', verticalalignment='center')

    if colorbar_label is None:
        colorbar_label = param
        if param == 'alpha':
            colorbar_label = r"$\alpha \ [-]$"
        elif param == 'ai':
            colorbar_label = r"$a_{i} \ [m^{-1}]$"
        elif param == 'Dsm1':
            colorbar_label = r"$D_{sm,1} \ [mm]$"
        elif param == 'ug1':
            colorbar_label = r"$v_{g} \ [m/s]$"

    # Only add colorbar when desired (for grids, usually make one shared bar)
    if (not suppress_colorbar) and add_colorbar:
        # robust ticks
        span = set_max - set_min
        if span > 0:
            tx_step = span / 5
            tx = np.linspace(set_min, set_max, 6)
        else:
            tx = [set_min]
        cbar = fig.colorbar(mpbl, ax=ax, label=colorbar_label, ticks=tx)

    if title_str != '':
        title = True
    if title:
        ax.set_title(title_str if title_str else cond.name)

    if created_fig:
        fig.tight_layout()
        if show:
            plt.show()
        else:
            os.makedirs(save_dir, exist_ok=True)
            fig.savefig(os.path.join(save_dir, f"{param}_contours_{cond.name + extra_save_text}.png"))
            plt.close(fig)

    return ax, mpbl, cbar


def plot_contour2_grid(
    database,
    param: str,
    row_org: str,
    col_org: str,
    *,
    cartesian=False,
    fig_size=3.0,
    ngridr=50,
    ngridphi=50,
    rot_angle=0,
    colormap='plasma',
    num_levels=0,
    level_step=0.01,
    set_min=None,
    set_max=None,
    font_size=10,
    show=True,
    save_path=None,
    shared_colorbar=True,
    colorbar_label=None,
    # optional formatting
    sort_key=None,  # function to sort unique org values; default: numeric then str
):
    """
    Builds a grid where rows are unique getattr(cond, row_org) and columns are unique getattr(cond, col_org).
    Calls plot_contour2 for each cell that has a matching condition.
    """

    def _default_sort(vals):
        # try numeric sort; fall back to string
        try:
            return sorted(vals, key=lambda x: float(x))
        except Exception:
            return sorted(vals, key=lambda x: str(x))

    if sort_key is None:
        sort_key = _default_sort

    # Collect unique row/col values
    row_vals = sort_key({getattr(c, row_org) for c in database})
    col_vals = sort_key({getattr(c, col_org) for c in database})

    nrows = len(row_vals)
    ncols = len(col_vals)

    subplot_kw = {} if cartesian else dict(projection='polar')
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(fig_size * ncols, fig_size * nrows),
        dpi=300,
        subplot_kw=subplot_kw,
        squeeze=False
    )

    # Map conditions to grid cells (if multiple match, last one wins)
    lookup = {}
    for cond in database:
        r = getattr(cond, row_org)
        c = getattr(cond, col_org)
        lookup[(r, c)] = cond

    # If you want consistent scaling across all plots, compute global min/max from data.
    # This is optional; if you pass set_min/set_max, those are used directly.
    # Here’s a lightweight strategy: evaluate each cond quickly by using its raw values
    # rather than re-interpolating; good enough for shared colorbar scaling.
    if set_min is None or set_max is None:
        all_vals = []
        for cond in database:
            # attempt to collect raw values for param
            for _, r_dict in cond.data.items():
                for r, midas_output in r_dict.items():
                    if r >= 0:
                        v = midas_output.get(param, np.nan) if hasattr(midas_output, "get") else midas_output.get(param, np.nan) if isinstance(midas_output, dict) else np.nan
                        all_vals.append(v)
        all_vals = np.asarray(all_vals, dtype=float)
        if set_min is None:
            set_min = np.nanmin(all_vals)
        if set_max is None:
            vmax = np.nanmax(all_vals)
            set_max = vmax + 0.1 * vmax

    mappable = None

    for i, rv in enumerate(row_vals):
        for j, cv in enumerate(col_vals):
            ax = axes[i, j]
            cond = lookup.get((rv, cv), None)

            # label outer edges
            if i == 0:
                ax.set_title(f"{col_org} = {cv}")
            if j == 0:
                # y-label space for leftmost column
                ax.text(-0.15, 0.5, f"{row_org} = {rv}", transform=ax.transAxes,
                        rotation=90, va='center', ha='right')

            if cond is None:
                # blank cell
                ax.axis('off')
                continue

            _, mpbl, _ = plot_contour2(
                cond, param,
                show=False,
                cartesian=cartesian,
                fig_size=fig_size,
                suppress_colorbar=True,    # no per-axes bars
                add_colorbar=False,
                rot_angle=rot_angle,
                ngridr=ngridr,
                ngridphi=ngridphi,
                colormap=colormap,
                num_levels=num_levels,
                level_step=level_step,
                set_min=set_min,
                set_max=set_max,
                font_size=font_size,
                ax=ax,
                fig=fig,
                title=False
            )
            mappable = mpbl  # keep last; OK for shared colorbar

    fig.tight_layout()

    if shared_colorbar and (mappable is not None):
        if colorbar_label is None:
            colorbar_label = param
        fig.colorbar(mappable, ax=axes, fraction=0.02, pad=0.02, label=colorbar_label)

    if save_path is not None:
        os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        fig.savefig(save_path, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

    return fig, axes