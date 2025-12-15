from .config import *

def plot_profiles2(cond, param, save_dir = '.', show=True, x_axis='vals', errorbars = False, 
                    const_to_plot = [90, 67.5, 45, 22.5, 0], include_complement = True, skip_1_comp = False,
                    fig_size=(4,4), fs = 10, title=True, label_str = '', legend_loc = 'best', xlabel_loc = 'center', include_const = False,
                    set_min = None, set_max = None, show_spines = True, xlabel_loc_coords = None, ylabel_loc_coords = None, cs=None, ms = None, ls = None) -> None:
    """Line plots of params
    
    **Args:**

        - ``param``: ``midas_dict`` parameter to plot. See :func:`~MARIGOLD.Condition.print_params` for options. Can be a single string, or a list of ``param`` strings.
        - ``save_dir``: directory in which to save the .png file. Will not save the file unless show = False. Defaults to '.'.
        - ``show``: display the figure (in an iPython notebook or have it pop up). Defaults to True.
        - ``x_axis``: the variable to put on the x-axis. Usually for vertical flow this is ``'rs'``, for horizontal, ``'vals'``. Also can be 'phis'. Defaults to 'vals'.
        - ``errorbars``: percentage errorbars to include. Can also specify ``'sigma'`` if you know what you're doing. Defaults to False.
        - ``const_to_plot``: a list of angles (if x-axis is 'vals' or 'rs') or rs (if x-axis is 'phis'). Defaults to [90, 67.5, 45, 22.5, 0].
        - ``include_complement``: includes the complementary angle of the const_to_plot (i.e. 270 with 90). Defaults to True.
        - ``skip_1_comp``: skip the r/R=1.0 value for the complementary angle. Useful for not having an ugly line interpolated down to 0 for r/R=1 if the data actually stops at something like r/R=-0.2. Defaults to False.
        - ``fig_size``: figure size tuple, in inches. Defaults to (4,4).
        - ``fs``: font size. Defaults to 10.
        - ``title``: option to display title. Defaults to True.
        - ``label_str``: param-axis label. Defaults to param name.
        - ``legend_loc``: option passed to ``plt.legend``. Defaults to 'best'.
        - ``xlabel_loc``: option passed to ``plt.legend``. Defaults to 'center'.
        - ``include_const``: Includes the constant angle in the plot legend. Defaults to False.
        - ``set_min``: minimum value for plot. If not specified, based on the data.
        - ``set_max``: maximum value for plot. If not specified, based on the data.
        - ``show_spines``: draw a box around the plot. Defaults to True.
        - ``xlabel_loc_coords``: coordinates to move the xlabel to. Defaults to None.
        - ``ylabel_loc_coords``: coordinates to move the ylabel to. Defaults to None.
        - ``cs``: colors to passed to ``plt.plot``. Can be a single value or list. If a list, will cycle through. Defaults to None.
        - ``ms``: marker style to passed to ``plt.plot``. Can be a single value or list. If a list, will cycle through. Defaults to None.
        - ``ls``: line style to passed to ``plt.plot``. Can be a single value or list. If a list, will cycle through. Defaults to None.
    """

    # TODO rewrite so it always loops over a list of params. If there's only one, just put it in a list at the begininng 

    plt.rcParams.update({'font.size': fs})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    fig, ax = plt.subplots(figsize=fig_size, dpi=300, layout='compressed')

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    # Tick marks facing in
    ax.tick_params(direction='in',which='both')

    if type(errorbars) is float and errorbars > 0:
        ax.plot([], [], ' ', label = f"{errorbars*100:0.1f}% error bars") # dummy to just get this text in the legend

    if not ms:
        ms = marker_cycle()
    else:
        ms = marker_cycle(marker_list=ms)

    if not ls:
        ls = line_cycle()
    else:
        ls = line_cycle(line_list=ls)

    if type(cs) == list:
        cs = color_cycle(color_list= cs)
    elif cs is None or type(param) == list:
        cs = color_cycle()
    elif cs == 'infer' and type(param) == str:
        cs = color_cycle(set_color = param)

    if set_min == None and type(param) == str:
        set_min = cond.min(param)
    elif set_min == None and type(param) == list:
        mins = []
        for specific_param in param:
            mins.append(cond.min(specific_param))

        set_min = min(mins)
    
    if set_max == None and type(param) == str:
        set_max = cond.max(param) *1.1

    elif set_min == None and type(param) == list:
        maxs = []
        for specific_param in param:
            maxs.append(cond.max(specific_param))

        set_max = max(maxs) *1.1

    if x_axis == 'vals' or x_axis == 'rs':
        for angle in const_to_plot:
            r_dict = cond.data[angle]
            rs = []
            vals = []
            errs = []
            for r, midas_output in r_dict.items():
                rs.append(r)

                if type(param) == str:
                    try:
                        vals.append(midas_output[param])
                    except:
                        if abs(r - 1) < 0.0001:
                            vals.append(0.0)
                        else:
                            vals.append(0.0)
                            print(f"Could not find {param} for φ = {angle}, r = {r}. Substituting 0")
                    
                    if errorbars == 'sigma':
                        if param == 'vr':
                            errs.append(midas_output['sigma_vr'])
                        elif param == 'vf':
                            errs.append(midas_output['sigma_vf'])
                        else:
                            print('issue with sigma, assuming 0 error')
                            errs.append(0)
                    elif type(errorbars) is float:
                            errs.append(abs(midas_output[param]*errorbars))
                    else:
                        errs.append(0)

                elif type(param) == list:
                    for i, specific_param in enumerate(param):
                        vals.append([])
                        try:
                            vals[i].append(midas_output[specific_param])
                        except:
                            if abs(r - 1) < 0.0001:
                                vals[i].append(0.0)
                            else:
                                vals[i].append(0.0)
                                print(f"Could not find {specific_param} for φ = {angle}, r = {r}. Substituting 0")

            if include_complement:
                if angle > 180:
                    print("Error: Cannot find complement to angle > 180. Skipping")
                else:
                    r_dict = cond.data[angle+180]
                    for r, midas_output in r_dict.items():
                        if skip_1_comp and r > 0.95:
                            # print('skipping')
                            continue
                        
                        rs.append(-r)
                        if type(param) == str:
                            try:
                                vals.append(midas_output[param])
                            except:
                                if abs(r - 1) < 0.0001:
                                    vals.append(0.0)
                                else:
                                    vals.append(0.0)
                                    print(f"Could not find {param} for φ = {angle}, r = {r}. Substituting 0")

                            if errorbars == 'sigma':
                                if param == 'vr':
                                    errs.append(midas_output['sigma_vr'])
                                elif param == 'vf':
                                    errs.append(midas_output['sigma_vf'])
                                else:
                                    print('issue with sigma, assuming 0 error')
                                    errs.append(0)
                            elif type(errorbars) is float:
                                errs.append(abs(midas_output[param]*errorbars))
                            else:
                                errs.append(0)
                        
                        elif type(param) == list:
                            for i, specific_param in enumerate(param):
                                vals.append([])
                                try:
                                    vals[i].append(midas_output[specific_param])
                                except:
                                    if abs(r - 1) < 0.0001:
                                        vals[i].append(0.0)
                                    else:
                                        vals[i].append(0.0)
                                        print(f"Could not find {specific_param} for φ = {angle}, r = {r}. Substituting 0")
            if type(param) == str:
                vals = [var for _, var in sorted(zip(rs, vals))]
                errs = [err for _, err in sorted(zip(rs, errs))]
                rs = sorted(rs)

                # errs = errorbars * np.abs(np.asarray(vals))

                if x_axis == 'vals':
                    if (type(errorbars) == float and errorbars == 0) or errorbars == False:
                        ax.plot(vals, rs, label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
                    else:
                        ax.errorbar(vals, rs, xerr = errs, capsize=3, ecolor = "black", label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
                elif x_axis == 'rs':
                    if type(errorbars) == float and errorbars == 0 or errorbars == False:
                        ax.plot(rs, vals, label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
                    else:
                        ax.errorbar(rs, vals, yerr = errs, capsize=3, ecolor = "black", label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
            
            elif type(param) == list:
                temp = []
                for vals_list in vals:
                    temp_list = [var for _, var in sorted(zip(rs, vals_list))]
                    temp.append(temp_list)
                rs = sorted(rs)

                for specific_param, val_list in zip(param, temp):
                    errs = errorbars * np.abs(np.asarray(val_list))
                    # print(val_list)
                    if cs == 'infer':
                        cs = color_cycle(set_color = specific_param)
                    
                    
                    legend_str = ''
                    if '_' not in specific_param:
                        if 'alpha' in specific_param:
                            specific_param = specific_param.replace('alpha', r'$\alpha$')
                        elif 'ai' in specific_param:
                            specific_param = specific_param.replace('ai', r'$a_{i}$')
                        elif 'ug1' in specific_param:
                            specific_param = specific_param.replace('ug1', r'$v_{g}$')
                        elif 'vf' in specific_param:
                            specific_param = specific_param.replace('vf', r'$v_{f}$')
                        
                        elif specific_param == 'vr':
                            legend_str = r'$v_{r, G1}$'
                        elif specific_param == 'vr2':
                            legend_str = r'$v_{r, G2}$'
                        elif specific_param == 'Dsm1':
                            legend_str = r'$D_{sm1}$'
                        elif specific_param == 'Dsm2':
                            legend_str = r'$D_{sm2}$'

                    else:
                        if 'ug1' in specific_param:
                            split_param = specific_param.split('_')
                            legend_str = '$' + 'v' + '_{g, ' + split_param[1] + '}$'
                        elif 'vf' in specific_param:
                            split_param = specific_param.split('_')
                            legend_str = '$' + 'v' + '_{f, ' + split_param[1] + '}$'
                        elif 'alpha' in specific_param:
                            split_param = specific_param.split('_')
                            legend_str = '$' + r'\alpha' + '_{' + split_param[1] + '}$'
                        else:
                            split_param = specific_param.split('_')
                            legend_str = '$' + split_param[0] + '_{' + split_param[1] + '}$'

                    if specific_param == 'vg_approx':
                        legend_str = r'$v_{g, model}$'

                    if legend_str == '':
                        legend_str = specific_param

                    if include_const:
                        legend_str += f", {angle}°"

                    if x_axis == 'vals':
                        if type(errorbars) == float and errorbars == 0 or errorbars == False:
                            ax.plot(val_list, rs, color=next(cs), marker=next(ms), linestyle = '--', label = legend_str)
                        else:
                            ax.errorbar(val_list, rs, xerr = errs, capsize=3, ecolor = "black", color=next(cs), marker=next(ms), linestyle = '--', label = legend_str)
                    elif x_axis == 'rs':
                        if type(errorbars) == float and errorbars == 0 or errorbars == False:
                            ax.plot(rs, val_list, color=next(cs), marker=next(ms), linestyle = '--', label = legend_str)
                        else:
                            ax.errorbar(rs, val_list, yerr = errs, capsize=3, ecolor = "black", color=next(cs), marker=next(ms), linestyle = '--', label = legend_str)
                
        
        if x_axis == 'vals':
            ax.set_ylim(-1, 1)
            ax.set_xlim(set_min, set_max)
        elif x_axis == 'rs':
            ax.set_xlim(-1, 1)
            ax.set_ylim(set_min, set_max)
        
    
    elif x_axis == 'phi':
        ax.plot([], [], label=r'$r/R$', color='white', linestyle = None)
        for rtarget in const_to_plot:
            phis = []
            vals = []
            for angle, r_dict in cond.data.items():
                for rstar, midas_output in r_dict.items():
                    if abs(rstar - rtarget) < 0.001:
                        phis.append(angle)
                        try:
                            vals.append(midas_output[param])
                        except:
                            if abs(rstar - 1) < 0.0001:
                                vals.append(0.0)
                            else:
                                vals.append(0.0)
                                print(f"Could not find {param} for φ = {angle}, r = {rstar}. Substituting 0")

            vals = [var for _, var in sorted(zip(phis, vals))]
            phis = sorted(phis)
            
            ax.plot(phis, vals, label=f'{rtarget:0.1f}', color=next(cs), marker=next(ms), linestyle = '--')
        
        ax.set_xlim(0, 360)
        ax.set_ylim(set_min, set_max)
            
    else:
        print(f"invalid axis for plot_profiles: {x_axis}. Current supported options are 'rs', 'vals' and 'phi'")
        return
    
    if label_str == '' and type(param) == str:
        if param == 'alpha':
            label_str = r'$\alpha\ [-]$'
        elif param == 'ai':
            label_str = r'$a_{i}\ [1/m]$'
        elif param == 'ug1' or param == 'ug2':
            label_str = r'$v_{g}\ [m/s]$'
        elif param == 'vf':
            label_str = r'$v_{f}\ [m/s]$'
        elif param == 'vr' or param == 'vr2':
            label_str = r'$v_{r}\ [m/s]$'
        elif param == 'Dsm1':
            label_str = r'$D_{sm1}\ [mm]$'
        else:
            label_str = param
    
    if x_axis == 'vals':
        ax.set_xlabel(label_str, loc = xlabel_loc)
        ax.set_ylabel(r'$r/R$ [-]')
        ax.set_yticks(np.arange(-1, 1.01, 0.2))
        
        ax.spines['bottom'].set_position(('data', 0))
        
        if set_min == 0 or set_max == 0:
            ax.spines['left'].set_position(('data', set_min))
            ax.spines['right'].set_position(('data', 0))
            # ax.yaxis.tick_right()
            # ax.yaxis.set_label_position("right")

        if not show_spines:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        else:
            ax2 = ax.twiny()
            ax2.get_xaxis().set_visible(False)


    elif x_axis == 'rs':
        ax.set_ylabel(label_str, loc = xlabel_loc)
        ax.set_xlabel(r'$r/R$ [-]')
        ax.set_xticks(np.arange(-1, 1.01, 0.2))

        ax.spines['bottom'].set_position(('data', max(0, set_min)))
        
        ax.spines['left'].set_position(('data', 0))
        
        if not show_spines:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        else:
            ax2 = ax.twiny()
            ax2.get_xaxis().set_visible(False)
        
        
    elif x_axis == 'phi':
        
        ax.set_ylabel(label_str, loc = xlabel_loc)
        ax.set_xlabel(r'$\varphi$ [-]')

        ax.set_xticks([0, 90, 180, 270, 360])
            
    if hasattr(title,'__len__'):
        # Set your own title! -DHK
        ax.set_title(title)
    elif title:
        ax.set_title(cond.name)

    
    

    leg = ax.legend(loc=legend_loc, edgecolor='white')

    if xlabel_loc_coords:
        ax.xaxis.set_label_coords(*xlabel_loc_coords)

    if ylabel_loc_coords:
        ax.yaxis.set_label_coords(*ylabel_loc_coords)


    
    ax.set_aspect('auto', adjustable='datalim', share=True)
    #fake_ax.set_box_aspect(1)
    
    if show:
        plt.show()
    else:
        if type(param) == str:
            plt.savefig(os.path.join(save_dir, f'{param}_profile_vs_{x_axis}_{cond.name}.png'))
        else:
            plt.savefig(os.path.join(save_dir, f'{"_".join(param)}_profile_vs_{x_axis}_{cond.name}.png'))
        plt.close()
    return

def plot_profiles(cond, param, save_dir = '.', show=True, x_axis='vals', 
                    const_to_plot = [90, 67.5, 45, 22.5, 0], include_complement = True, 
                    rotate=False, fig_size=(4,4), title=True, label_str = '', legend_loc = 'best', xlabel_loc = 'center',
                    set_min = None, set_max = None, show_spines = True, force_RH_y_axis = False, xlabel_loc_coords = None, cs=None) -> None:
    """Line plot of ``param``. Only a single ``param`` can be plotted at a time, but plot can be arbitrarily rotated. 

    **Args:**

        - ``param``: ``midas_dict`` parameter to plot. See :func:`~MARIGOLD.Condition.print_params` for options
        - ``save_dir``: _description_. Defaults to '.'.
        - ``show``: _description_. Defaults to True.
        - ``x_axis``: _description_. Defaults to 'vals'.
        - ``const_to_plot``: _description_. Defaults to [90, 67.5, 45, 22.5, 0].
        - ``include_complement``: _description_. Defaults to True.
        - ``rotate``: _description_. Defaults to False.
        - ``fig_size``: _description_. Defaults to (4,4).
        - ``title``: _description_. Defaults to True.
        - ``label_str``: _description_. Defaults to ''.
        - ``legend_loc``: _description_. Defaults to 'best'.
        - ``xlabel_loc``: _description_. Defaults to 'center'.
        - ``set_min``: _description_. Defaults to None.
        - ``set_max``: _description_. Defaults to None.
        - ``show_spines``: _description_. Defaults to True.
        - ``force_RH_y_axis``: _description_. Defaults to False.
        - ``xlabel_loc_coords``: _description_. Defaults to None.
        - ``cs``: _description_. Defaults to None.
    """

    cond.mirror()
    plt.rcParams.update({'font.size': 10})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    log_x = False # This breaks, so I removed it from the arguments to the function

    

    if rotate:

        # Set up the figure to be rotated by theta
        import matplotlib as mpl
        from matplotlib.transforms import Affine2D
        import mpl_toolkits.axisartist.floating_axes as floating_axes
        fig = plt.figure(figsize=fig_size)
        if x_axis == 'r':
            plot_extents = cond.min(param), cond.max(param)*1.1, -1, 1
            transform = Affine2D().scale(fig_size[0] / (cond.max(param)*1.1 - cond.min(param)), fig_size[1] / (1 - -1)).rotate_deg(cond.theta)
        else:
            plot_extents = cond.min(param), cond.max(param)*1.1, 0, 360
            transform = Affine2D().scale(fig_size[0] / (cond.max(param)*1.1 - cond.min(param)), fig_size[1] / (360-0)).rotate_deg(cond.theta)
        
        
        helper = floating_axes.GridHelperCurveLinear(transform, plot_extents)
        fake_ax = floating_axes.FloatingSubplot(fig, 111, grid_helper=helper)
        ax = fake_ax.get_aux_axes(transform)

    else:
        fig, ax = plt.subplots(figsize=fig_size, dpi=300, layout='compressed')
        fake_ax = ax

    # Only show ticks on the left and bottom spines
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    # Tick marks facing in
    ax.tick_params(direction='in',which='both')

    ms = marker_cycle()
    if cs is None:
        cs = color_cycle()
    elif cs == 'infer':
        cs = color_cycle(set_color = param)
    else:
        print("I hope cs is a generator that returns valid colors")

    if set_min == None:
        set_min = cond.min(param)
    
    if set_max == None:
        set_max = cond.max(param) *1.1

    if x_axis == 'vals' or x_axis == 'r':
        for angle in const_to_plot:
            r_dict = cond.data[angle]
            rs = []
            vals = []
            for r, midas_output in r_dict.items():
                rs.append(r)
                try:
                    vals.append(midas_output[param])
                except:
                    if abs(r - 1) < 0.0001:
                        vals.append(0.0)
                    else:
                        vals.append(0.0)
                        print(f"Could not find {param} for φ = {angle}, r = {r}. Substituting 0")

            if include_complement:
                if angle > 180:
                    print("Error: Cannot find complement to angle > 180. Skipping")
                else:
                    r_dict = cond.data[angle+180]
                    for r, midas_output in r_dict.items():
                        rs.append(-r)
                        try:
                            vals.append(midas_output[param])
                        except:
                            if abs(r - 1) < 0.0001:
                                vals.append(0.0)
                            else:
                                vals.append(0.0)
                                print(f"Could not find {param} for φ = {angle}, r = {r}. Substituting 0")

            vals = [var for _, var in sorted(zip(rs, vals))]
            rs = sorted(rs)
            if x_axis == 'vals':
                ax.plot(vals, rs, label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
            if x_axis == 'r':
                ax.plot(rs, vals, label=f'{angle}°', color=next(cs), marker=next(ms), linestyle = '--')
        
        if x_axis == 'vals':
            ax.set_ylim(-1, 1)
            ax.set_xlim(set_min, set_max)
        elif x_axis == 'r':
            ax.set_xlim(-1, 1)
            ax.set_ylim(set_min, set_max)
            fake_ax.set_xlim(-1, 1)
            fake_ax.set_ylim(set_min, set_max)
        
    
    elif x_axis == 'phi':
        ax.plot([], [], label=r'$r/R$', color='white', linestyle = None)
        for rtarget in const_to_plot:
            phis = []
            vals = []
            for angle, r_dict in cond.data.items():
                for rstar, midas_output in r_dict.items():
                    if abs(rstar - rtarget) < 0.001:
                        phis.append(angle)
                        try:
                            vals.append(midas_output[param])
                        except:
                            if abs(rstar - 1) < 0.0001:
                                vals.append(0.0)
                            else:
                                vals.append(0.0)
                                print(f"Could not find {param} for φ = {angle}, r = {rstar}. Substituting 0")

            vals = [var for _, var in sorted(zip(phis, vals))]
            phis = sorted(phis)
            if rotate:
                ax.plot(vals, phis, label=f'{rtarget:0.1f}', color=next(cs), marker=next(ms), linestyle = '--')
            else:
                ax.plot(phis, vals, label=f'{rtarget:0.1f}', color=next(cs), marker=next(ms), linestyle = '--')
        if rotate:
            ax.set_ylim(0, 360)
            ax.set_xlim(set_min, set_max)
            
        else:
            ax.set_xlim(0, 360)
            ax.set_ylim(set_min, set_max)
            
    else:
        print(f"invalid axis for plot_profiles: {x_axis}. Current supported options are 'r' and 'phi'")
        return
    
    if label_str == '':
        if param == 'alpha':
            label_str = r'$\alpha\ [-]$'
        elif param == 'ai':
            label_str = r'$a_{i}\ [1/m]$'
        elif param == 'ug1':
            label_str = r'$v_{g}\ [m/s]$'
        else:
            label_str = param
    
    if x_axis == 'vals':
        fake_ax.set_xlabel(label_str, loc = xlabel_loc)
        fake_ax.set_ylabel(r'$r/R$ [-]')
        fake_ax.set_yticks(np.arange(-1, 1.01, 0.2))
        #fake_ax.set_xticks(np.linspace(cond.min(param), cond.max(param), 7))

    elif x_axis == 'r':
        fake_ax.set_ylabel(label_str, loc = xlabel_loc)
        fake_ax.set_xlabel(r'$r/R$ [-]')
        fake_ax.set_xticks(np.arange(-1, 1.01, 0.2))

    elif x_axis == 'phi':
        if not rotate:
            fake_ax.set_ylabel(label_str, loc = xlabel_loc)
            fake_ax.set_xlabel(r'$\varphi$ [-]')

            fake_ax.set_xticks([0, 90, 180, 270, 360])
            #fake_ax.set_yticks(np.linspace(cond.min(param), cond.max(param), 7))
        else:
            fake_ax.set_xlabel(label_str, loc = xlabel_loc)
            fake_ax.set_ylabel(r'$\varphi$ [-]')
            
            fake_ax.set_yticks([0, 90, 180, 270, 360])
            # fake_ax.set_xticks(np.linspace(cond.min(param), cond.max(param), 7))

            #fake_ax.tick_params(axis='both', labelrotation=-cond.theta)
    
    if title:
        ax.set_title(cond.name)

    ax.spines['bottom'].set_position(('data', 0))


    if set_min == 0 or set_max == 0 or force_RH_y_axis:
        ax.spines['left'].set_position(('data', set_min))
        ax.spines['right'].set_position(('data', 0))
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")

    if not show_spines:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    else:
        ax2 = ax.twiny()
        ax2.get_xaxis().set_visible(False)


    if log_x:
        ax.set_xlim(cond.min(param, nonzero=True), cond.max(param)*1.2)
        ax.set_xscale('log')
        fake_ax.set_xscale('log')
    
    if rotate:
        fig.add_subplot(fake_ax)
    ax.legend(loc=legend_loc, edgecolor='white')

    if xlabel_loc_coords:
        ax.xaxis.set_label_coords(*xlabel_loc_coords)

    fake_ax.set_aspect('auto', adjustable='datalim', share=True)
    ax.set_aspect('auto', adjustable='datalim', share=True)
    #fake_ax.set_box_aspect(1)
    
    if show:
        plt.show()
    else:
        plt.savefig(os.path.join(save_dir, f'{param}_profile_vs_{x_axis}_{cond.name}.png'))
        plt.close()
    return   

def plot_isoline(cond, param:str, iso_axis:str, iso_val:float, fig_size=4, plot_res=100, 
                    save_dir = '.', show=True, extra_text = '') -> None:
    """_summary_
    
    **Args:**

        - ``param``: ``midas_dict`` parameter to plot. See :func:`~MARIGOLD.Condition.print_params` for options
        - ``iso_axis``: _description_
        - ``iso_val``: _description_
        - ``fig_size``: _description_. Defaults to 4.
        - ``plot_res``: _description_. Defaults to 100.
        - ``save_dir``: _description_. Defaults to '.'.
        - ``show``: _description_. Defaults to True.
        - ``extra_text``: _description_. Defaults to ''.
    """

    fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300, layout='compressed')

    plt.rcParams.update({'font.size': 10})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    cond.mirror()

    if iso_axis == 'x':
        lim = np.sqrt(1 - iso_val**2)
        ys = np.linspace(-lim, lim, plot_res)
        xs = iso_val * np.ones(ys.size)
        plt.plot(ys, cond(xs, ys, param, interp_method = 'linear_xy'), color = cond.marker_color, marker=None, linestyle= cond.line_style)

        plt.xlim(-1, 1)
        plt.xlabel(r'$y/R$ [-]')
        plt.ylabel(param)
        plt.show()

    elif iso_axis == 'y':
        lim = np.sqrt(1 - iso_val**2)
        xs = np.linspace(-lim, lim, plot_res)
        ys = iso_val * np.ones(xs.size)
        plt.plot(xs, cond(xs, ys, param, interp_method = 'linear_xy'), color = cond.marker_color, marker=None, linestyle= cond.line_style)

        plt.xlim(-1, 1)
        plt.xlabel(r'$x/R$ [-]')
        plt.ylabel(param)
        plt.show()

    if show:
        plt.show()
    else:
        plt.savefig( os.path.join(save_dir, f"{param}_iso_{iso_axis}={iso_val}_{cond.name + extra_text}.png") )
        plt.close()
    return

def plot_contour(cond, param:str, save_dir = '.', show=True, set_max = None, set_min = None, fig_size = 4, colorbar_label = None, suppress_colorbar = False,
                    rot_angle = 0, ngridr = 50, ngridphi = 50, colormap = 'plasma', num_levels = 0, level_step = 0.01, title = False, title_str = '', extra_save_text = '',
                    annotate_h = False, cartesian = False, h_star_kwargs = {'method': 'max_dsm', 'min_void': '0.05'}, plot_measured_points = False, font_size = 12) -> None:
    """Function to create a contour plot of a given param

    **Args:**

        - ``param``: ``midas_dict`` parameter to plot. See :func:`~MARIGOLD.Condition.print_params` for options
        - ``save_dir``: directory to save contour plot to. Defaults to '.'.
        - ``show``: whether or not to show the contour plot. Defaults to True.
        - ``set_max``: to specify the maximum value of the contour plot. If None, will caclulate based on data. Defaults to None.
        - ``set_min``: to specify the minimum value of the contour plot. If None, will caclulate based on data. Defaults to None.
        - ``fig_size``: size, in inches to make to figure square. Defaults to 4.
        - ``colorbar_label``: label to apply to the colorbar. Will default to the parameter name, with some common ones prettied up with LaTeX formatting. Defaults to None.
        - ``suppress_colorbar``: Don't include colorbar. Defaults to False.
        - ``rot_angle``: Rotate the contour plot by a specific angle (in degrees). Defaults to 0.
        - ``ngridr``: the number of interpolation points in the :math:`r` direction. Defaults to 50.
        - ``ngridphi``: the number of interpolation points in the :math:`\\varphi` direction. Defaults to 50.
        - ``colormap``: colormap to use for the contour plot. Defaults to 'plasma'.
        - ``num_levels``: number of levels for the contours. Defaults to 0.
        - ``level_step``: steps to define the number of contours. Not used if ``num_levels`` specified. Defaults to 0.01.
        - ``title``: title for the top of the plot. Default title is ``cond.name``. Defaults to False.
        - ``title_str``: string to use as the title, if specified title will be set to True. Defaults to ''.
        - ``extra_save_text``: extra text to include while saving. Defaults to ''.
        - ``annotate_h``: draw a line where :math:`h` is calculated to be. Set ``cartestian=True`` for best results. Defaults to False.
        - ``cartesian``: plot in Cartesian coordinates. Defaults to False.
        - ``h_star_kwargs``: for annotate_h. Defaults to {'method': 'max_dsm', 'min_void': '0.05'}.
        - ``plot_measured_points``: to plot red circles where the original data was measured (before mirroring). Defaults to False.
        - ``font_size``: font size. Defaults to 12.

    Returns:
        - ``ax``, the matplotlib axis the contour plot was made with
    """

    if cartesian:
        fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300)
    else:
        fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300, subplot_kw=dict(projection='polar'))
    plt.rcParams.update({'font.size': font_size})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"
    
    # cond.mirror()
    rs = []
    phis = []
    vals = []

    for phi_angle, r_dict in cond.data.items():
        for r, midas_output in r_dict.items():
            if r >= 0:
                rs.append(r)
                phis.append(phi_angle)
                try:
                    vals.append(midas_output[param])
                except:
                    if abs(r - 1) < 0.0001:
                        vals.append(0.0)
                    else:
                        vals.append(0.0)
                        print(f"Could not find {param} for φ = {phi_angle}, r = {r}. Substituting 0")

    rs = np.asarray(rs)
    phis = (np.asarray(phis) + rot_angle) * np.pi / 180 
    vals = np.asarray(vals)   

    ri = np.linspace(0, 1, ngridr)
    phii = np.linspace(rot_angle * np.pi / 180, 2*np.pi + rot_angle * np.pi / 180, ngridphi)

    triang = tri.Triangulation(phis, rs)
    interpolator = tri.LinearTriInterpolator( triang, vals )

    PHII, RI = np.meshgrid(phii, ri)
    XI = RI * np.cos(PHII)
    YI = RI * np.sin(PHII)
    parami = interpolator(PHII, RI)

    if debug and False: 
        print(cond.name, f'{param} to contour plot', file=debugFID)
        print(Xi, Yi, alphai, file=debugFID)
        print(f'max {param}', np.amax(alphai), file = debugFID)
    
    if param == 'alpha' and (np.nanmax(vals) > 1.0 or np.nanmax(vals) > 1.0):
        np.savetxt(f'contour_error_dump_{cond.name}.csv', vals, delimiter=',')
        #print(vals)
        raise(ValueError('Alpha exceeds 1.0!\nSaving problematic array'))
        print('Alpha exceeds 1.0!\nSaving problematic array')
        return
    
    extend_min = False
    extend_max = False
    
    if set_min == None:
        set_min = np.min(parami)

    if set_max == None:
        set_max = np.max(parami) + (np.max(parami) * 0.1)

    if set_min > np.min(parami):
        extend_min = True

    if set_max < np.max(parami):
        extend_max = True

    if extend_max and extend_min:
        extend_opt = 'both'
    elif extend_min and not extend_max and set_min > np.min(parami):
        extend_opt = 'min'
    elif extend_max and not extend_min:
        extend_opt = 'max'
    elif not extend_min and not extend_max:
        extend_opt = 'neither'

    if abs(set_max - set_min) < level_step:
        print(f'Warning: Level step {level_step} too larger for range {set_max}, {set_min}. Defaulting to 0.01*(set_max - set_min)')
        level_step = 0.01*(set_max - set_min)

    if num_levels:
        lvs = np.linspace(set_min, set_max, num_levels)
    else:
        if (abs(set_max) < 1e-8): # if the max is 0, start the counting from there
            lvs = np.arange(set_max, abs(set_min) + 1e-8, level_step)
            lvs = np.flip(lvs) * -1
        else:
            lvs = np.arange(set_min, set_max + 1e-8, level_step)
    
    if cartesian:
        mpbl = ax.contourf(XI, YI, parami, levels = lvs, vmin = set_min, vmax = set_max, cmap = colormap)

        x = np.linspace(-1, 1, 100)
        # Make a circle to look nice
        plt.plot(x, np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
        plt.plot(x, -np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
    else:
        mpbl = ax.contourf(PHII, RI, parami, levels = lvs, 
                            vmin = set_min, vmax = set_max, cmap = colormap, extend=extend_opt)
        if plot_measured_points: 
            # print(np.asarray(cond.original_mesh)[:,0]* np.pi/180 , np.asarray(cond.original_mesh)[:,1])
            measured_thetas = np.asarray(cond.original_mesh)[:,0]* np.pi/180
            measured_rs = np.asarray(cond.original_mesh)[:,1]

            for i, r in enumerate(np.asarray(cond.original_mesh)[:,1]):
                if r < 0:
                    measured_rs[i] = - measured_rs[i]
                    measured_thetas[i] = (measured_thetas[i] + np.pi) % (2*np.pi)

            ax.plot( measured_thetas, measured_rs, marker='o', fillstyle = 'none', mec = 'red', linewidth = 0, ms = 2)

    if annotate_h:
        if not cartesian:
            print("Warning: annotate_h assumes that we're plotting in Cartesian. Annotation may not be correct")
        hstar = cond.find_hstar_pos(**h_star_kwargs)
        ax.annotate('h', 
                (-1,1-hstar), 
                (1,1-hstar),
                arrowprops=dict(color='r', width=3, headwidth=3),
                color='r',
                verticalalignment='center'
    )
        
    ax.grid(False)
    if cartesian:
        plt.axis('square')
        plt.xlabel (r'$x/R$ [-]')
        plt.ylabel(r'$y/R$ [-]')
        ax.grid(False)
    else:
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    #plt.clim(vmin=set_min, vmax=set_max)
    
    if colorbar_label == None:
        colorbar_label = param
        if param == 'alpha':
            colorbar_label = r"$\alpha \ [-]$"
        elif param == 'ai':
            colorbar_label = r"$a_{i} \ [m^{-1}]$"
        elif param == 'Dsm1':
            colorbar_label = r"$D_{sm,1} \ [mm]$"
        elif param == 'ug1':
            colorbar_label = r"$v_{g} \ [m/s]$"

    if not suppress_colorbar:
        tx_step = round((set_max - set_min)/5,-int(np.floor(np.log10((set_max - set_min)/10))))
        tx = np.arange(set_min,set_max,tx_step)

        fig.colorbar(mpbl, label=colorbar_label, ticks=tx)

    if title_str != '':
        title = True
    
    if title:
        if title_str == '':
            plt.title(cond.name)
        else:
            plt.title(title_str)

    plt.tight_layout()
    
    #cb = plt.colorbar(ticks = [0, 0.05, 0.1, 0.15, 0.2])

    if show:
        plt.show()
    else:
        plt.savefig( os.path.join(save_dir, f"{param}_contours_{cond.name + extra_save_text}.png") )
        plt.close()
    return ax

def plot_surface(cond, param:str, save_dir = '.', show=True, set_max = None, set_min = None, rotate_gif=False, elev_angle = 145, 
                    azim_angle = 0, roll_angle = 180, title=True, ngridr = 50, ngridphi = 50, 
                    plot_surface_kwargs = None, solid_color = False, label_str = None, title_str = '', colormap = 'viridis') -> None:
    """Function to create a 3 dimensional surface plot of a given parameter
    
    **Args:**

        - ``param``: ``midas_dict`` parameter to plot. See :func:`~MARIGOLD.Condition.print_params` for options
        - ``save_dir``: directory to save surface plot image to. Defaults to '.'.
        - ``show``: option to show surface plot. Defaults to True.
        - ``set_max``: option to set the maximum value of the surface plot. Defaults to None.
        - ``set_min``: option to set the minimum value of the surface plot. Defaults to None.
        - ``rotate_gif``: option to produce a .gif file where the surface plot rotates around. Defaults to False.
        - ``elev_angle``: elevation angle to rotate image of surface plot. Defaults to 145.
        - ``azim_angle``: azimuthal angle to rotate image of surface plot. Defaults to 0.
        - ``roll_angle``: roll angle to rotate image of surface plot. Defaults to 180.
        - ``title``: title of plot, if ``title_str`` is not set, will use. Defaults to True.
        - ``ngridr``: the number of interpolation points in the :math:`r` direction. Defaults to 50.
        - ``ngridphi``: the number of interpolation points in the :math:`\\varphi` direction. Defaults to 50.
        - ``plot_surface_kwargs``: dictionary to pass additional parameters to ``plot_surface``. Defaults to None.
        - ``solid_color``: use a solid color instead of a colormap. Defaults to False.
        - ``label_str``: custom label for z-axis. Defaults to None.
        - ``title_str``: title of surface plot. Defaults to ''.
        - ``colormap``: colormap to use for surface plot. Defaults to 'viridis'.
    """

    if plot_surface_kwargs is None:
        plot_surface_kwargs = {}
    plt.rcParams.update({'font.size': 16})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    cond.mirror()
    rs = []
    phis = []
    vals = []

    for phi_angle, r_dict in cond.data.items():
        for r, midas_output in r_dict.items():
            if r >= 0:
                rs.append(r)
                phis.append(phi_angle)

                try:
                    vals.append(midas_output[param])
                except:
                    vals.append(np.NaN)
                    print(f"Could not find {param} for φ = {phi_angle}, r = {r}. Substituting NaN")

    rs = np.asarray(rs)
    phis = (np.asarray(phis)) * np.pi / 180 
    vals = np.asarray(vals)   

    ri = np.linspace(0, 1, ngridr)
    phii = np.linspace(0, 2*np.pi, ngridphi)

    triang = tri.Triangulation(phis, rs)
    interpolator = tri.LinearTriInterpolator( triang, vals )

    PHII, RI = np.meshgrid(phii, ri)
    parami = interpolator(PHII, RI)

    Xi = RI * np.cos(PHII)
    Yi = RI * np.sin(PHII)

    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw={"projection": "3d"})

    if 'vmin' not in plot_surface_kwargs.keys(): 
        plot_surface_kwargs.update({'vmin': cond.min(param)})

    if 'vmax' not in plot_surface_kwargs.keys(): 
        plot_surface_kwargs.update({'vmax': cond.max(param)})

    if 'cmap' not in plot_surface_kwargs.keys() and not solid_color: 
        plot_surface_kwargs.update({'cmap': colormap})

    surf = ax.plot_surface(Xi, Yi, parami, **plot_surface_kwargs)
    
    #plt.legend()
    ax.set_xlabel (r'$x/R$ [-]')
    ax.set_ylabel(r'$y/R$ [-]')

    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    
    ax.set_zlim([plot_surface_kwargs['vmin'], plot_surface_kwargs['vmax']])
    
    if set_min == None:
        set_min = np.min(parami)

    if set_max == None:
        set_max = np.max(parami) + (np.max(parami) * 0.1)

    tx_step = round((set_max - set_min)/5,-int(np.floor(np.log10((set_max - set_min)/10))))
    tx = np.arange(set_min,set_max,tx_step)
    
    if label_str:
        ax.set_zlabel(label_str)
        fig.colorbar(surf, label=label_str,ticks=tx)
    else:
        ax.set_zlabel(param)
        fig.colorbar(surf, label=param,ticks=tx)
    
    if title: 
        if title_str:
            plt.title(title_str)
        else:
            plt.title(cond.name)

    ax.view_init(azim=azim_angle, roll=roll_angle, elev=elev_angle)
            
    #cb = plt.colorbar(ticks = [0, 0.05, 0.1, 0.15, 0.2])
    if show:
        plt.show()
    else:
        plt.savefig( os.path.join(save_dir, f"{param}_surface_{cond.name}.png") )
    
    if rotate_gif:
        import matplotlib.animation as animation
        
        def rotate(angle):
            ax.view_init(azim=angle)
        
        rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
        rot_animation.save(os.path.join(save_dir, f'{param}_surface_rotation_{cond.name}.gif'), dpi=80)

    return

def plot_spline_contour(cond, param:str, save_dir = '.', show=True, set_max = None, set_min = None, fig_size = 4,
                    rot_angle = 0, ngridr = 50, ngridphi = 50, colormap = 'plasma', num_levels = 100, title = False,
                    annotate_h = False, cartesian = False, h_star_kwargs = {'method': 'max_dsm', 'min_void': '0.05'},
                    grad = 'None') -> None:
    """Plots a contour from a spline interpolation

    Will fit the spline if necessary. 
    
    """ 
    
    if param not in cond.spline_interp.keys():
        print(f"Warning: {param} not found in spline_interp dict, running fit_spline")
        cond.fit_spline(param)

    phii = np.linspace(0, 2*np.pi, 100)
    phii_arg = phii + rot_angle * np.pi/180
    ri = np.linspace(0, 1, 100)

    if cartesian:
        fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300)
    else:
        fig, ax = plt.subplots(figsize=(fig_size, fig_size), dpi=300, subplot_kw=dict(projection='polar'))
    plt.rcParams.update({'font.size': 10})
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["mathtext.fontset"] = "cm"

    PHII, RI = np.meshgrid(phii, ri)
    XI = RI * np.cos(PHII)
    YI = RI * np.sin(PHII)

    if grad == 'None':
        VALS = (cond.spline_interp[param](phii_arg, ri)).T
    elif grad == 'r':
        VALS = (cond.spline_interp[param](phii_arg, ri, dy=1)).T
    elif grad == 'r':
        VALS = (cond.spline_interp[param](phii_arg, ri, dy=1)).T
    elif grad == 'phi':
        VALS =  (1./RI *cond.spline_interp[param](phii_arg, ri, dx = 1)).T
    elif grad == 'phinor':
        VALS =  (cond.spline_interp[param](phii_arg, ri, dx = 1)).T
    elif grad == 'y':
        VALS = (cond.spline_interp['alpha'](phii_arg, ri, dx=1) * np.cos(phii_arg)*ri / (ri**2+1e-8) + cond.spline_interp['alpha'](phii_arg, ri, dy=1) * np.sin(phii_arg)).T
    else:
        print(f"Error: unrecognized grad type {grad}")


    if cartesian:
        plt.contourf(XI, YI, VALS, levels = num_levels, vmin = set_min, vmax = set_max, cmap = colormap)

        x = np.linspace(-1, 1, 100)
        # Make a circle to look nice
        plt.plot(x, np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
        plt.plot(x, -np.sqrt(1- x**2), marker= None, linestyle = '-', color = 'black', linewidth = 1)
    else:
        plt.contourf(PHII, RI, VALS, levels = num_levels, vmin = set_min, vmax = set_max, cmap = colormap)

    if annotate_h:
        if not cartesian:
            print("Warning: annotate_h assumes that we're plotting in Cartesian. Annotation may not be correct")
        hstar = cond.find_hstar_pos(**h_star_kwargs)
        ax.annotate('h', 
                (-1,1-hstar), 
                (1,1-hstar),
                arrowprops=dict(color='r', width=3, headwidth=3),
                color='r',
                verticalalignment='center'
    )
        
    ax.grid(False)
    if cartesian:
        plt.axis('square')
        plt.xlabel (r'$x/R$ [-]')
        plt.ylabel(r'$y/R$ [-]')
        ax.grid(True)
    else:
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    if grad == 'None':
        param_label = param
    else:
        param_label = "grad_" + param + "_" + grad 
    plt.colorbar(label=param_label)
    
    if title:
        plt.title(cond.name)

    plt.tight_layout()
    
    #cb = plt.colorbar(ticks = [0, 0.05, 0.1, 0.15, 0.2])

    if show:
        plt.show()
    else:
        plt.savefig( os.path.join(save_dir, f"{param_label}_spline_contours_{cond.name}.png") )
        plt.close()
    return

def color_cycle(set_color = None, color_list = []):
    """Custom generator for colors
    
    set_color can be:
     - None, for a basic cycle of blue, red, green, etc.
     - A single hexcolor code ('#000000' for black, etc)
     - 'alpha', 'ai', 'ug1', 'Dsm1' or 'vr', which have default built in colors
     - Otherwise, it assumes set_color is a list of colors to yield

    """
    
    if not color_list:
        if set_color == None:
            color_list = ['#FF0000',    # Red
                        '#0000FF',      # Blue
                        '#00FF00',      # Green
                        '#FF00FF',      # Magenta
                        '#00FFFF',      # Cyan
                        '#FFA500',      # Orange
                        '#000000',      # Black
                        '#7F00FF',      # Violet
                        '#007F7F',      # Teal
                        '#7F7F7F',      # Gray
                        '#008000',      # Dark Green
                        '#7FFF7F']      # Light Green
        elif re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', set_color):
            color_list = [set_color]
        # I want these to be able to override each other
        else:
            if 'alpha' in set_color:
                color_list = ['#000000']
            if 'alpha_G1' in set_color:
                color_list = ['#A0A0A0']
            if 'alpha_G2' in set_color:
                color_list = ['#606060']
            if 'ai' in set_color:
                color_list = ['#00FF00']
            if 'ai_G2' in set_color:
                color_list = ['#66FF66']
            if 'ai_G1' in set_color:
                color_list = ['#006600']
            if 'ug1' in set_color:
                color_list = ['#FF0000']
            if 'ug2' in set_color:
                color_list = ['#ff8080']
            if 'Dsm1' in set_color:
                color_list = ['#00FFFF']
            if 'Dsm2' in set_color:
                color_list = ['#6666FF']
            if 'vf' in set_color or 'vf_approx' in set_color:
                color_list = ['#0000FF']
            if 'vr' in set_color:
                color_list = ['#FF00FF']
            if 'vr2' in set_color:
                color_list = ['#FF80FF']
            
            elif not color_list:
                color_list = ['#000000']
        

    i = 0
    while True:
        yield color_list[ i % len(color_list)]
        i += 1

def marker_cycle(marker_list = []):
    """Custom generator for markers

    """
    if not marker_list:
        marker_list = ['o', '^', 's', 'v', 'D']
    i = 0
    while True:
        yield marker_list[ i % len(marker_list)]
        i += 1

def line_cycle(line_list = []):
    """Custom generator for linestyles

    """

    if not line_list:
        line_list = ['-', '--', '-.', ':']
    
    i = 0
    while True:
        yield line_list[ i % len(line_list)]
        i += 1