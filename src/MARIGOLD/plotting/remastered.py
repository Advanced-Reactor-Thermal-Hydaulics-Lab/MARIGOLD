from ..config import *

_DEFAULT_COLORS = [
    "#E84393",
    "#6BBBE8",
    "#7ED957",
    "#8E44AD",
    "#F39C12",
]
_DEFAULT_MARKERS = ["o"]
_DEFAULT_LINESTYLES = ["--"]

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

    return fig, ax, plotted_any

def _cond_match_key(cond):
    return (cond.theta, cond.jf, cond.jgref, cond.port)

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
    if len(databases) != len(db_labels):
        raise ValueError("databases and db_labels must have the same length.")

    if len(databases) < 2:
        raise ValueError("Provide at least two databases to compare.")

    indexed = []
    for db in databases:
        d = {}
        for cond in db:
            d[_cond_match_key(cond)] = cond
        indexed.append(d)

    common_keys = set(indexed[0].keys())
    for d in indexed[1:]:
        common_keys &= set(d.keys())
    matched_keys = sorted(common_keys)

    # If nothing matches, return nothing
    if not matched_keys:
        return []

    n_db = len(databases)

    if db_styles is None:
        color_cyc = cycle(_DEFAULT_COLORS)
        marker_cyc = cycle(_DEFAULT_MARKERS)
        ls_cyc = cycle(_DEFAULT_LINESTYLES)

        db_styles = []
        for _ in range(n_db):
            db_styles.append({
                "color": next(color_cyc),
                "marker": next(marker_cyc),
                "linestyle": next(ls_cyc),
            })
    else:
        if len(db_styles) < n_db:
            raise ValueError(
                f"db_styles has {len(db_styles)} entries but {n_db} databases were provided."
            )

    # Global rcParams for consistent typography
    plt.rcParams.update({"font.size": fs})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    bold_title = FontProperties(weight="bold")

    for key in matched_keys:
        conds = [d[key] for d in indexed]

        fig, ax = plt.subplots(figsize=fig_size, dpi=300, layout="compressed")
        ax.tick_params(direction="in", which="both")
        ax.yaxis.set_ticks_position("left")
        ax.xaxis.set_ticks_position("bottom")

        db_proxy_handles = []
        db_proxy_labels = []

        phi_legend_added = False
        plotted_all = True

        # Plot each database on the same axes
        for db_i, (cond, label) in enumerate(zip(conds, db_labels)):
            style = db_styles[db_i]

            before_lines = len(ax.lines)

            fig, ax, plotted_any = plot_xy_profiles(
                cond,
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
                legend_phi=not phi_legend_added,
                title=True,
                fig_size=fig_size,
                fs=fs,
                db_label=label,
                style=style,
            )

            after_lines = len(ax.lines)
            added_any_artists = after_lines > before_lines

            # STRICT: if this DB didn't actually add plotted data, kill the whole figure
            if not (plotted_any and added_any_artists):
                plotted_all = False
                break

            # Only add φ legend once (the first DB that successfully plotted something)
            if not phi_legend_added:
                phi_legend_added = True

            # Proxy handle for DB legend (neutral black line with DB linestyle)
            db_proxy_handles.append(
                Line2D([0], [0],
                       color="black",
                       linestyle=style.get("linestyle", "-"),
                       marker="None",
                       linewidth=1.5)
            )
            db_proxy_labels.append(label)

        # If strict match requirement fails, do not output/save a plot
        if not plotted_all:
            plt.close(fig)
            continue

        # db_leg = ax.legend(
        #     db_proxy_handles,
        #     db_proxy_labels,
        #     loc="lower left",
        #     title="Database",
        #     title_fontproperties=bold_title,
        #     edgecolor="white",
        #     # bbox_to_anchor=(0.0, 1.0),
        # )

        # phi_leg = ax.legend(
        #     loc="upper right",
        #     title=r"$\varphi$",
        #     title_fontproperties=bold_title,
        #     edgecolor="white",
        #     # bbox_to_anchor=(0.0, 0.72),
        # )

        phi_leg = ax.legend(
            labels=db_proxy_labels,
            loc="lower left",
            title=r"Database",
            title_fontproperties=bold_title,
            edgecolor="white",
        )

        # ax.add_artist(db_leg)

        # Save/show
        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)
            theta, jf, jgref, port = key
            fname = f"{x_param}_vs_{y_param}_theta{theta}_jf{jf}_jg{jgref}_{port}.png"
            fig.savefig(os.path.join(save_dir, fname))
            plt.close(fig)
        else:
            if show:
                plt.show()
            else:
                plt.close(fig)

    return