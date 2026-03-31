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

def plot_params(
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
    set_xlim=None,
    set_ylim=None,
    percent_error: float = 0.0,
    legend_phi: bool = True,
    title: bool = True,
    fig_size=(4, 4),
    fs: int = 10,
    show: bool = True,
    close: bool = True,
    db_label: str | None = None,    # currently unused
    style=None,
):
    created_fig = ax is None

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
        return phi_val if phi_val <= 180 else phi_val - 180

    def _resolve_style_for_phi(phi_val: float):
        base = _base_phi(phi_val)

        if base not in style_by_base_phi:
            color = next(color_cyc)
            marker = next(marker_cyc)
            ls = next(ls_cyc)

            if style:
                color = style.get("color", color)
                marker = style.get("marker", marker)
                ls = style.get("linestyle", ls)

            if not show_lines:
                ls = "None"
            if not show_markers:
                marker = "None"

            style_by_base_phi[base] = (color, marker, ls)

        return base, style_by_base_phi[base]

    plotted_any = False

    for phi in phis_to_plot:
        if phi not in cond.data:
            continue

        r_dict = cond.data[phi]
        rs_sorted = sorted(r_dict.keys())

        xs, ys = [], []

        for r in rs_sorted:
            if abs(r - 1.0) < 1e-6:   # r* = 1.0 dummy
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
            yerr = np.abs(ys) * percent_error
            ax.errorbar(
                xs,
                ys,
                yerr=yerr,
                capsize=3,
                ecolor="black",
                color=color,
                marker=marker,
                linestyle=ls,
                markeredgecolor="black",
                markeredgewidth=0.8,
                label=phi_label,
            )
        else:
            ax.plot(
                xs,
                ys,
                color=color,
                marker=marker,
                linestyle=ls,
                markeredgecolor="black",
                markeredgewidth=0.8,
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

    if legend_phi and plotted_any:
        ax.legend(loc="best", edgecolor="white")

    if created_fig:
        if show:
            plt.show()
        elif close:
            plt.close(fig)

    return fig, ax, plotted_any

def plot_stack(
    conds,
    x_param: str,
    y_param: str,
    labels=None,
    phis_to_plot=(90, 270),
    *,
    styles=None,
    ax=None,
    show_markers: bool = True,
    show_lines: bool = True,
    x_label: str | None = None,
    y_label: str | None = None,
    set_xlim=None,
    set_ylim=None,
    percent_error: float = 0.0,
    legend_phi: bool = False,
    legend_stack: bool = True,
    title: bool = True,
    fig_size=(4, 4),
    fs: int = 10,
    show: bool = True,
    close: bool = True,
):
    if labels is None:
        labels = [str(cond) for cond in conds]

    if len(labels) != len(conds):
        raise ValueError(
            f"labels must have the same length as conds; got {len(labels)} and {len(conds)}"
        )

    created_fig = ax is None
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

    # Build default per-condition styles if none are provided
    if styles is None:
        color_cyc = cycle(_DEFAULT_COLORS)
        marker_cyc = cycle(_DEFAULT_MARKERS)
        ls_cyc = cycle(_DEFAULT_LINESTYLES)

        styles = []
        for _ in conds:
            st = {
                "color": next(color_cyc),
                "marker": next(marker_cyc),
                "linestyle": next(ls_cyc),
            }

            if not show_lines:
                st["linestyle"] = "None"
            if not show_markers:
                st["marker"] = "None"

            styles.append(st)
    else:
        if len(styles) != len(conds):
            raise ValueError(
                f"styles must have the same length as conds; got {len(styles)} and {len(conds)}"
            )

        styles = [dict(s) if s is not None else {} for s in styles]
        for st in styles:
            if not show_lines:
                st["linestyle"] = "None"
            if not show_markers:
                st["marker"] = "None"

    plotted_any = False

    # Plot each condition by delegating to plot_params()
    for cond, label, style in zip(conds, labels, styles):
        _, _, cond_plotted = plot_params(
            cond,
            x_param=x_param,
            y_param=y_param,
            phis_to_plot=phis_to_plot,
            ax=ax,
            show_markers=show_markers,
            show_lines=show_lines,
            x_label=x_label,
            y_label=y_label,
            set_xlim=set_xlim,
            set_ylim=set_ylim,
            percent_error=percent_error,
            legend_phi=legend_phi,
            title=False,   # stack-level title handled below
            fig_size=fig_size,
            fs=fs,
            show=False,
            close=False,
            style=style,
        )

        # Add proxy handle for condition-level legend
        if cond_plotted and legend_stack:
            ax.plot([], [], label=label, **style)

        plotted_any = plotted_any or cond_plotted

    # Labels
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

    # Title
    if title and len(conds) > 0:
        ref = conds[0]
        ax.set_title(
            rf"$x = {x_param},\; y = {y_param}$"
            + "\n"
            + rf"$\theta = {ref.theta}^\circ,\; j_f = {ref.jf}\,[m/s],\; j_{{gref}} = {ref.jgref}\,[m/s]$"
        )

    # Legends
    if legend_phi and legend_stack and plotted_any:
        handles_all, labels_all = ax.get_legend_handles_labels()

        stack_label_set = set(labels)
        phi_handles, phi_labels = [], []
        stack_handles, stack_labels = [], []

        for h, lab in zip(handles_all, labels_all):
            if lab in stack_label_set:
                stack_handles.append(h)
                stack_labels.append(lab)
            else:
                phi_handles.append(h)
                phi_labels.append(lab)

        if phi_handles:
            leg1 = ax.legend(
                phi_handles,
                phi_labels,
                loc="best",
                edgecolor="white",
                title="Phi",
            )
            ax.add_artist(leg1)

        if stack_handles:
            ax.legend(
                stack_handles,
                stack_labels,
                loc="upper right",
                edgecolor="white",
                title="Condition",
            )

    elif legend_phi and plotted_any:
        ax.legend(loc="best", edgecolor="white")

    elif legend_stack and plotted_any:
        ax.legend(loc="best", edgecolor="white")

    if created_fig:
        if show:
            plt.show()
        elif close:
            plt.close(fig)

    return fig, ax, plotted_any

def plot_vary(
    database,
    variable: str,
    x_param: str,
    y_param: str,
    phis_to_plot=(90, 270),
    *,
    styles=None,
    ax=None,
    show_markers: bool = True,
    show_lines: bool = True,
    x_label: str | None = None,
    y_label: str | None = None,
    set_xlim=None,
    set_ylim=None,
    percent_error: float = 0.0,
    legend_phi: bool = False,
    legend_stack: bool = True,
    title: bool = True,
    fig_size=(4, 4),
    fs: int = 10,
    show: bool = True,
    close: bool = True,
    sort_vary_values: bool = True,
):
    """
    Group conditions by all tracked attributes except `variable`, and for each group
    call plot_stack() so that only `variable` changes within a figure.

    Returns
    -------
    results : list[tuple]
        List of (fig, ax, plotted_any, conds_group) for each generated figure.
    """

    vary_attrs = ["theta", "jf", "jgref", "port", "database"]

    def _format_attr(cond, attr):
        if attr == "theta":
            return rf"$\theta = {cond.theta}^\circ$"
        if attr == "jf":
            return rf"$j_f = {cond.jf:.2f}\,[m/s]$"
        if attr == "jgref":
            return rf"$j_{{gref}} = {cond.jgref:.2f}\,[m/s]$"
        if attr == "port":
            return f"{cond.port}"
        if attr == "database":
            return f"Author: {cond.database}"
        return f"{attr}={getattr(cond, attr)}"

    if variable not in vary_attrs:
        raise ValueError(
            f"variable must be one of {vary_attrs}; got {variable!r}"
        )

    if ax is not None:
        raise ValueError("plot_vary() creates multiple figures, so ax must be None.")

    fixed_attrs = [attr for attr in vary_attrs if attr != variable]

    def _key(cond):
        return tuple(getattr(cond, attr) for attr in fixed_attrs)

    def _vary_value(cond):
        return getattr(cond, variable)

    def _sort_key(cond):
        val = _vary_value(cond)
        try:
            return (0, float(val))
        except Exception:
            return (1, str(val))

    groups = {}
    for cond in database:
        key = _key(cond)
        if key not in groups:
            groups[key] = []
        groups[key].append(cond)

    results = []

    for _, conds_group in groups.items():
        if sort_vary_values:
            conds_group = sorted(conds_group, key=_sort_key)

        labels = [_format_attr(cond, variable) for cond in conds_group]

        fig, ax_out, plotted_any = plot_stack(
            conds=conds_group,
            x_param=x_param,
            y_param=y_param,
            labels=labels,
            phis_to_plot=phis_to_plot,
            styles=styles,
            ax=None,
            show_markers=show_markers,
            show_lines=show_lines,
            x_label=x_label,
            y_label=y_label,
            set_xlim=set_xlim,
            set_ylim=set_ylim,
            percent_error=percent_error,
            legend_phi=legend_phi,
            legend_stack=legend_stack,
            title=False,
            fig_size=fig_size,
            fs=fs,
            show=False,
            close=False,
        )

        if title and len(conds_group) > 0:
            cond0 = conds_group[0]
            parts = [_format_attr(cond0, attr) for attr in fixed_attrs]
            ax_out.set_title(", ".join(parts))

        if show:
            plt.show()
        elif close:
            plt.close(fig)

        results.append((fig, ax_out, plotted_any, conds_group))

    return results

def plot_contour2(
    cond,
    param: str,
    show=True,
    close=True,
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
    ax=None,
    fig=None,
    add_colorbar=True,
):
    """
    Draw a contour plot for `param` from a condition object.

    If `ax` is None, creates its own figure/axes.
    If `ax` is provided, draws on that axes and does not show/close the figure.

    Returns
    -------
    ax : matplotlib axes
    mpbl : QuadContourSet
    cbar : Colorbar or None
    """
    created_fig = False
    cbar = None

    if ax is None:
        created_fig = True
        if cartesian:
            fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300)
        else:
            fig, ax = plt.subplots(
                figsize=(fig_size, fig_size),
                dpi=300,
                subplot_kw=dict(projection='polar')
            )
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

    rs = np.asarray(rs, dtype=float)
    phis = (np.asarray(phis, dtype=float) + rot_angle) * np.pi / 180.0
    vals = np.asarray(vals, dtype=float)

    ri = np.linspace(0, 1, ngridr)
    phii = np.linspace(
        rot_angle * np.pi / 180.0,
        2 * np.pi + rot_angle * np.pi / 180.0,
        ngridphi
    )

    triang = tri.Triangulation(phis, rs)
    interpolator = tri.LinearTriInterpolator(triang, vals)

    PHII, RI = np.meshgrid(phii, ri)
    XI = RI * np.cos(PHII)
    YI = RI * np.sin(PHII)
    parami = interpolator(PHII, RI)

    # Fill masked -> NaN for nanmin/nanmax logic
    parami_f = parami.filled(np.nan) if np.ma.isMaskedArray(parami) else np.asarray(parami)

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
        mpbl = ax.contourf(
            XI, YI, parami_f,
            levels=lvs,
            vmin=set_min,
            vmax=set_max,
            cmap=colormap
        )

        x = np.linspace(-1, 1, 200)
        ax.plot(x, np.sqrt(1 - x**2), color='black', linewidth=1)
        ax.plot(x, -np.sqrt(1 - x**2), color='black', linewidth=1)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel(r'$x/R$ [-]')
        ax.set_ylabel(r'$y/R$ [-]')
        ax.grid(False)
    else:
        mpbl = ax.contourf(
            PHII, RI, parami_f,
            levels=lvs,
            vmin=set_min,
            vmax=set_max,
            cmap=colormap,
            extend=extend_opt
        )
        ax.grid(False)
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    if annotate_h:
        hstar = cond.find_hstar_pos(**h_star_kwargs)
        ax.annotate(
            'h',
            (-1, 1 - hstar),
            (1, 1 - hstar),
            arrowprops=dict(color='r', width=3, headwidth=3),
            color='r',
            verticalalignment='center'
        )

    if plot_measured_points:
        if cartesian:
            pts_x = rs * np.cos(phis)
            pts_y = rs * np.sin(phis)
            ax.plot(pts_x, pts_y, 'k.', markersize=2)
        else:
            ax.plot(phis, rs, 'k.', markersize=2)

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

    if (not suppress_colorbar) and add_colorbar:
        span = set_max - set_min
        if span > 0:
            tx = np.linspace(set_min, set_max, 6)
        else:
            tx = [set_min]

        cbar = fig.colorbar(mpbl, ax=ax, label=colorbar_label, ticks=tx)

    if title_str != '':
        title = True
    if title:
        ax.set_title(title_str if title_str else str(cond))

    if created_fig:
        fig.tight_layout()
        if show:
            plt.show()
        elif close:
            plt.close(fig)

    return ax, mpbl, cbar

def _make_key_fn(org_spec):
    """
    Returns a function cond -> key used to bin into rows/cols.

    org_spec can be:
      - "port" (attribute name)
      - "pair" (special: (jf, jgref))
      - ("jf","jgref") (tuple/list of attributes)
      - callable(cond) -> key
    """
    if callable(org_spec):
        return org_spec

    if isinstance(org_spec, str):
        s = org_spec.lower().strip()

        if s == "pair":
            attrs = ("jf", "jgref")
            return lambda cond: tuple(getattr(cond, a) for a in attrs)

        # extend with other shorthands if you want:
        # if s == "pair_loc": attrs = ("jf","jgref","port")

        return lambda cond: getattr(cond, org_spec)

    if isinstance(org_spec, (tuple, list)):
        attrs = tuple(org_spec)
        return lambda cond: tuple(getattr(cond, a) for a in attrs)

    raise TypeError(f"Unsupported org spec: {org_spec!r}")

def plot_contour_grid(
    database,
    param: str,
    row_org='pair',
    col_org='port',
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
    close=True,
    shared_colorbar=True,
    colorbar_label=None,
    sort_key=None,
    on_duplicate="last",
):
    import matplotlib.pyplot as plt
    import numpy as np
    from collections import defaultdict

    def _format_attr(cond, attr):
        if attr == "theta":
            return rf"$\theta = {cond.theta}^\circ$"
        if attr == "pair":
            return rf"$j_f = {cond.jf:.2f}\,[m/s], j_{{gref}} = {cond.jgref:.2f}\,[m/s]$"
        if attr == "port":
            return f"{cond.port}"
        if attr == "database":
            return f"Author: {cond.database}"
        return f"{attr}={getattr(cond, attr)}"

    def _default_sort(vals):
        try:
            return sorted(vals)
        except Exception:
            return sorted(vals, key=lambda x: str(x))

    def _coerce_numeric_array(x):
        """
        Return x as a 1D float array with only finite numeric values kept.
        Invalid entries like '#N/A', None, '', etc. are discarded.
        """
        if x is None:
            return np.array([], dtype=float)

        arr = np.asarray(x, dtype=object).ravel()
        cleaned = []

        for v in arr:
            try:
                fv = float(v)
                if np.isfinite(fv):
                    cleaned.append(fv)
            except (TypeError, ValueError):
                pass

        return np.asarray(cleaned, dtype=float)

    if sort_key is None:
        sort_key = _default_sort

    row_key = _make_key_fn(row_org)
    col_key = _make_key_fn(col_org)

    vary_attrs = ["theta", "pair", "port", "database"]

    if row_org not in vary_attrs:
        raise ValueError(f"row_org={row_org!r} is not supported")
    if col_org not in vary_attrs:
        raise ValueError(f"col_org={col_org!r} is not supported")
    if row_org == col_org:
        raise ValueError("row_org and col_org must be different")

    fixed_attrs = [attr for attr in vary_attrs if attr not in {row_org, col_org}]

    def _fixed_key(cond):
        vals = []
        for attr in fixed_attrs:
            if attr == "pair":
                vals.append((cond.jf, cond.jgref))
            else:
                vals.append(getattr(cond, attr))
        return tuple(vals)

    grouped = defaultdict(list)
    for cond in database:
        grouped[_fixed_key(cond)].append(cond)

    plt.rcParams.update({'font.size': font_size})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    # Compute global limits only from valid numeric entries
    if set_min is None or set_max is None:
        vals_all = []

        for cond in database:
            for _, r_dict in cond.data.items():
                for r, midas_output in r_dict.items():
                    if r < 0:
                        continue

                    try:
                        vals = _coerce_numeric_array(midas_output[param])
                        if vals.size > 0:
                            vals_all.extend(vals.tolist())
                    except Exception:
                        pass

        vals_all = np.asarray(vals_all, dtype=float)

        if vals_all.size > 0:
            if set_min is None:
                set_min = np.nanmin(vals_all)
            if set_max is None:
                vmax = np.nanmax(vals_all)
                set_max = vmax + 0.1 * vmax

    figures = {}

    for gkey, conds in grouped.items():
        row_vals = sort_key({row_key(c) for c in conds})
        col_vals = sort_key({col_key(c) for c in conds})

        nrows, ncols = len(row_vals), len(col_vals)
        subplot_kw = {} if cartesian else dict(projection='polar')

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(fig_size * ncols, fig_size * nrows),
            dpi=300,
            subplot_kw=subplot_kw,
            squeeze=False,
        )

        lookup = {}
        for cond in conds:
            rk = row_key(cond)
            ck = col_key(cond)
            key = (rk, ck)

            if key in lookup:
                if on_duplicate == "first":
                    continue
                if on_duplicate == "error":
                    raise ValueError(
                        f"Duplicate condition for cell {key}: {lookup[key]} and {cond}"
                    )

            lookup[key] = cond

        mappable = None

        for i, rv in enumerate(row_vals):
            for j, cv in enumerate(col_vals):
                ax = axes[i, j]
                cond = lookup.get((rv, cv))

                if i == 0:
                    ax.set_title(f"{cv}")

                if j == 0:
                    ax.text(
                        -0.15, 0.5, f"{rv}",
                        transform=ax.transAxes,
                        rotation=90,
                        va='center',
                        ha='right'
                    )

                if cond is None:
                    ax.axis('off')
                    continue

                try:
                    _, mpbl, _ = plot_contour2(
                        cond,
                        param,
                        show=False,
                        close=False,
                        cartesian=cartesian,
                        fig_size=fig_size,
                        suppress_colorbar=True,
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
                        title=False,
                    )
                    mappable = mpbl

                except Exception:
                    # This condition cannot be plotted for this param.
                    # Leave the cell empty.
                    ax.axis('off')
                    continue

        if conds:
            cond0 = conds[0]
            parts = [_format_attr(cond0, attr) for attr in fixed_attrs]
            fig.suptitle(", ".join(parts), y=0.995)

        fig.tight_layout(rect=[0, 0, 1, 0.97])

        if shared_colorbar and (mappable is not None):
            cbar_label = colorbar_label
            if cbar_label is None:
                cbar_label = param
                if param == 'alpha':
                    cbar_label = r"$\alpha \ [-]$"
                elif param == 'ai':
                    cbar_label = r"$a_{i} \ [m^{-1}]$"
                elif param == 'Dsm1':
                    cbar_label = r"$D_{sm,1} \ [mm]$"
                elif param == 'ug1':
                    cbar_label = r"$v_{g} \ [m/s]$"

            fig.colorbar(
                mappable,
                ax=axes,
                fraction=0.02,
                pad=0.02,
                label=cbar_label
            )

        figures[gkey] = (fig, axes)

        if show:
            plt.show()
        elif close:
            plt.close(fig)

    return figures