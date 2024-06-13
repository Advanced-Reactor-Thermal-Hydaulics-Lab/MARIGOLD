from .Condition import Condition
from .config import *
import subprocess

""" Functions for interfacing with CFX

"""

def write_CFX_BC(cond:Condition, save_dir = ".", z_loc = 'LoverD', only_90 = False, interp = False, csv_name = False, ngrid = 100):
    """ Write a csv file for CFX based on cond

    z_loc can be set by the user, or set to "LoverD" to use the cond L/D information

    vf information will be taken from condition, or approximated by approx_vf()

    only_90 option for writing data only down the 90 degree line
    
    """
    try:
        dummy = cond.run_ID
    except AttributeError:
        cond.run_ID = cond.database

    if not csv_name:
        if only_90:
            csv_name = f"{cond.run_ID}_{cond.theta}deg_jf{cond.jf:0.1f}_jg{cond.jgref}_{cond.port}_BC_90deg.csv"
        else:
            csv_name = f"{cond.run_ID}_{cond.theta}deg_jf{cond.jf:0.1f}_jg{cond.jgref}_{cond.port}_BC.csv"

    path_to_csv = os.path.join(save_dir, csv_name)

    if z_loc == 'LoverD':
        z_loc = cond.LoverD


    with open(path_to_csv, "w") as f:
        R = cond.Dh/2

        f.write("[Name],,,,,,,,,\n")
        f.write(f"{cond.port}data,,,,,,,,,\n")
        f.write("[Spatial Fields],,,,,,,,,\n")
        
        
        if not interp:
            f.write("radius,z,phi,,,,,,,\n")
            f.write("[Data],,,,,,,,,\n")
            f.write("radius [mm],z [m],phi [],Velocity u g[m s^-1],Velocity v g[m s^-1],Velocity w g[m s^-1],Volume Fraction [],Velocity u f [m s^-1],Velocity v f [m s^-1],Velocity w f [m s^-1],\n")

            for angle, r_dict in cond.phi.items():
                    for rstar, midas_output in r_dict.items():
                        if only_90:
                            if angle == 90 or angle == 270:
                                pass
                            else:
                                continue
                        r = rstar * 12.7 # r/R * R [mm]

                        f.write(f"{r},{z_loc},{angle * np.pi/180},{0},{0},{midas_output['ug1']},{midas_output['alpha']},{0},{0},{midas_output['vf']},\n")
        elif interp == 'xy':
            f.write("x,y,z,,,,,,,\n")
            f.write("[Data],,,,,,,,,\n")
            f.write("x [m],y [m],z [m],Velocity u g [m s^-1],Velocity v g [m s^-1],Velocity w g [m s^-1],Volume Fraction [],Velocity u f [m s^-1],Velocity v f [m s^-1],Velocity w f [m s^-1],\n")

            for x in np.linspace(-R / 2, R, ngrid):
                for y in np.linspace(-R, R, ngrid):
                    if np.sqrt(x**2 + y**2) <= R:
                        f.write(f"{x},{y},{z_loc},{0},{0},{ cond(x, y, 'ug1', interp_method='linear_xy') },{cond(x, y, 'alpha', interp_method='linear_xy')},{0},{0},{cond(x, y, 'vf', interp_method='linear_xy')},\n")

        else:
            f.write("radius,z,phi,,,,,,,\n")
            f.write("[Data],,,,,,,,,\n")
            f.write("radius [mm],z [m],phi [],Velocity u g [m s^-1],Velocity v g [m s^-1],Velocity w g [m s^-1],Volume Fraction [],Velocity u f [m s^-1],Velocity v f [m s^-1],Velocity w f [m s^-1],\n")

            for r in np.linspace(0, R, ngrid):
                for phi in np.linspace(0, 2*np.pi, ngrid):
                    f.write(f"{r*1000},{z_loc},{phi},{0},{0},{ cond(phi, r/R, 'ug1', interp_method='linear') },{cond(phi, r/R, 'alpha', interp_method='linear')},{0},{0},{cond(phi, r/R, 'vf', interp_method='linear')},\n")


    return

def read_CFX_export(csv_name, jf, jgref, theta, port, database, jgloc=None) -> Condition:
    """ Read CFX csv export into a MARIGOLD Condition object

    Must supply jf, jgref, theta, port, database, jgloc, etc. for Condition

    Returns Condition object
    
    """
    if jgloc is None:
        jgloc = jgref
    cond = Condition(jgref, jgloc, jf, theta, port, database)
    cond.run_ID = 'CFD'

    with open(csv_name) as fi:
        fi.readline()             # 
        fi.readline()             # [Name]
        fi.readline()             # port3
        fi.readline()             # 
        fi.readline()             # [data]
        variables = fi.readline() # variables

        variables = variables.split(",")

        vg_idx = [idx for idx, s in enumerate(variables) if 'gas.Velocity' in s][0]
        vf_idx = [idx for idx, s in enumerate(variables) if 'liquid.Velocity' in s][0]
        alpha_idx = [idx for idx, s in enumerate(variables) if 'gas.Volume' in s][0]
        x_idx = [idx for idx, s in enumerate(variables) if 'X [ m ]' in s][0]
        y_idx = [idx for idx, s in enumerate(variables) if 'Y [ m ]' in s][0]

        while True:
            try:
                data = fi.readline().split(",")
            except IOError:
                break
            
            if data == ['']:
                break

            try:
                x = float(data[x_idx])
                y = float(data[y_idx])
                vg = float(data[vg_idx])
                vf = float(data[vf_idx])
                alpha = float(data[alpha_idx])
            except Exception as e:
                print(e)
                print(variables)
                print("\nProblem data:")
                print(data)
                print(x_idx, y_idx, vg_idx, vf_idx, alpha_idx)

            data_dict = {'ug1': vg, 'vf': vf, 'alpha': alpha}
            
            roverR = np.sqrt(x**2 + y**2) / 0.0127
            phi_angle = int(np.arctan2(y, x) * 180/np.pi)

            try:
                cond.phi[phi_angle].update({roverR:data_dict})
            except:
                cond.phi.update({phi_angle:{}})
                cond.phi[phi_angle].update({roverR:data_dict})

    return cond


def make_ICEM_pipe_mesh(r_divs: int, theta_divs: int, z_divs: int, o_point: float, L: float, 
                        case_name: str, turb_model = 'ke', Ref = 20000, growth_ratio = 1.2, 
                        fluent_translator_path = "/apps/external/apps/ansys/2022r2/ansys_inc/v222/icemcfd/linux64_amd/icemcfd/output-interfaces/fluent6",
                        cleanup = True):
    """ Function for making an O-grid pipe mesh using ANSYS ICEM

    Adjusts first layer based on selected turbulence model, Ref

    o_point is the point where the mesh changes from the circular segments near the wall to the central square

    Writes a mesh replay file, runs ICEM using command check_call('module load ansys && icemcfd -x ./mesh_replay.rpl > auto_cfx_run.log', shell=True)
    
    """

    # Assuming water
    rho = 998
    mu = 0.001

    Uinf = Ref * mu / (rho * L)

    Cf = 0.026 / Ref**(1./7)
    tau_wall = Cf * rho * Uinf**2 / 2
    Ufric = (tau_wall / rho)**0.5

    if turb_model == 'ke': 
        yplus = 30
    else:
        yplus = 1

    first_layer = yplus*mu / (Ufric * rho)
    
    with open("mesh_replay.rpl", 'w') as fi:
        print(f'\
\
ic_set_global geo_cad 0 toptol_userset\n\
ic_set_global geo_cad 0.0 toler\n\
ic_undo_group_begin \n\
ic_geo_new_family GEOM\n\
ic_boco_set_part_color GEOM\n\
ic_empty_tetin \n\
ic_point {{}} GEOM pnt.00 0,0,0\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_point {{}} GEOM pnt.01 0,0,2.54\n\
ic_undo_group_end \n\
ic_set_global geo_cad 0 toptol_userset\n\
ic_set_global geo_cad 0.001 toler\n\
ic_set_global geo_cad 0.001 toler\n\
ic_set_global geo_cad 0.001 toler\n\
ic_undo_group_begin \n\
ic_surface cyl GEOM srf.00 {{pnt.00 pnt.01 0.0127 0.0127 1 1}}\n\
ic_set_global geo_cad 0.001 toler\n\
ic_set_dormant_pickable point 0 {{}}\n\
ic_set_dormant_pickable curve 0 {{}}\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_geo_new_family SOLID\n\
ic_boco_set_part_color SOLID\n\
ic_hex_initialize_blocking {{}} SOLID 0 101\n\
ic_hex_unblank_blocks \n\
ic_hex_multi_grid_level 0\n\
ic_hex_projection_limit 0\n\
ic_hex_default_bunching_law default 2.0\n\
ic_hex_floating_grid off\n\
ic_hex_transfinite_degree 1\n\
ic_hex_unstruct_face_type one_tri\n\
ic_hex_set_unstruct_face_method uniform_quad\n\
ic_hex_set_n_tetra_smoothing_steps 20\n\
ic_hex_error_messages off_minor\n\
ic_hex_set_piercing 0\n\
ic_undo_group_end \n\
ic_hex_find_comp_curve srf.00.C1\n\
ic_undo_group_begin \n\
ic_hex_undo_major_start set_edge_projection\n\
ic_hex_set_edge_projection 26 42 0 1 srf.00.C1\n\
ic_hex_set_edge_projection 38 42 0 1 srf.00.C1\n\
ic_hex_set_edge_projection 22 38 0 1 srf.00.C1\n\
ic_hex_set_edge_projection 22 26 0 1 srf.00.C1\n\
ic_hex_undo_major_end set_edge_projection\n\
ic_undo_group_end \n\
ic_hex_find_comp_curve srf.00.C0\n\
ic_undo_group_begin \n\
ic_hex_undo_major_start set_edge_projection\n\
ic_hex_set_edge_projection 25 41 0 1 srf.00.C0\n\
ic_hex_set_edge_projection 37 41 0 1 srf.00.C0\n\
ic_hex_set_edge_projection 21 37 0 1 srf.00.C0\n\
ic_hex_set_edge_projection 21 25 0 1 srf.00.C0\n\
ic_hex_undo_major_end set_edge_projection\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_project_to_surface GEOM SOLID\n\
ic_undo_group_end \n\
ic_hex_mark_blocks unmark\n\
ic_undo_group_begin \n\
ic_hex_mark_blocks superblock 13\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_mark_blocks face_neighbors corners {{ 22 38 26 42 }} {{ 21 37 25 41 }}\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_ogrid 1 m GEOM SOLID -version 50\n\
ic_hex_mark_blocks unmark\n\
ic_undo_group_end \n\
ic_hex_mark_blocks unmark\n\
\
ic_hex_set_node_location x {o_point} y {o_point} z 0 -csys global node_numbers {{  77  }}\n\
ic_hex_set_node_location x {-o_point} y {o_point} z 0 -csys global node_numbers {{  69  }}\n\
ic_hex_set_node_location x {-o_point} y {-o_point} z 0 -csys global node_numbers {{  65  }}\n\
ic_hex_set_node_location x {o_point} y {-o_point} z 0 -csys global node_numbers {{  73  }}\n\
\
ic_hex_set_node_location x {o_point} y {o_point} z 2.54 -csys global node_numbers {{  78  }}\n\
ic_hex_set_node_location x {-o_point} y {o_point} z 2.54 -csys global node_numbers {{  70  }}\n\
ic_hex_set_node_location x {-o_point} y {-o_point} z 2.54 -csys global node_numbers {{  66  }}\n\
ic_hex_set_node_location x {o_point} y {-o_point} z 2.54 -csys global node_numbers {{  74  }}\n\
\
ic_undo_group_begin \n\
ic_hex_set_mesh 26 42 n {theta_divs} h1 10000000000.0 h2 10000000000.0 r1 2 r2 2 lmax 1e+10 default unlocked\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_set_mesh 26 42 n {theta_divs} h1rel 501326007289.0 h2rel 501326007289.0 r1 2 r2 2 lmax 1e+10 default copy_to_parallel unlocked\n\
ic_undo_group_begin \n\
ic_undo_group_end \n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_set_mesh 38 42 n {theta_divs} h1rel 501326007289.0 h2rel 501326007289.0 r1 2 r2 2 lmax 1e+10 default copy_to_parallel unlocked\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_set_mesh 38 42 n {theta_divs} h1rel 501326007289.0 h2rel 501326007289.0 r1 2 r2 2 lmax 1e+10 default copy_to_parallel unlocked\n\
ic_undo_group_begin \n\
ic_undo_group_end \n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_set_mesh 26 70 n {r_divs} h1rel 0.0 h2rel 0.0 r1 2 r2 2 lmax 0 default copy_to_parallel unlocked\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_set_mesh 26 70 n {r_divs} h1rel {first_layer} h2rel 1.31233595801 r1 {growth_ratio} r2 2 lmax 0 default copy_to_parallel unlocked\n\
ic_undo_group_begin \n\
ic_undo_group_end \n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_set_mesh 25 26 n {z_divs} h1rel 3937007874.02 h2rel 3937007874.02 r1 2 r2 2 lmax 1e+10 default copy_to_parallel unlocked\n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_hex_set_mesh 25 26 n {z_divs} h1rel 3937007874.02 h2rel 3937007874.02 r1 2 r2 2 lmax 1e+10 default copy_to_parallel unlocked\n\
ic_undo_group_begin \n\
ic_undo_group_end \n\
ic_undo_group_end \n\
ic_hex_create_mesh GEOM SOLID proj 2 dim_to_mesh 3\n\
ic_hex_write_file ./hex.uns GEOM SOLID proj 2 dim_to_mesh 3 no_boco\n\
ic_uns_load ./hex.uns 3 0 {{}} 1\n\
ic_uns_update_family_type visible {{GEOM ORFN SOLID}} {{!LINE_2 QUAD_4 !HEXA_8}} update 0\n\
ic_boco_solver \n\
ic_boco_clear_icons \n\
ic_uns_update_family_type visible {{ORFN SOLID}} {{!LINE_2 QUAD_4 !HEXA_8}} update 0\n\
ic_uns_update_family_type visible {{GEOM ORFN SOLID}} {{!LINE_2 QUAD_4 !HEXA_8}} update 0\n\
ic_uns_create_selection_subset 0\n\
ic_uns_create_selection_edgelist 1\n\
ic_uns_subset_configure uns_sel_0 -draw_nodes 1\n\
ic_uns_update_family_type visible {{GEOM ORFN SOLID}} {{!LINE_2 !QUAD_4 !HEXA_8}} update 0\n\
ic_uns_subset_configure uns_sel_0 -width 0\n\
ic_uns_subset_configure uns_sel_0 -draw_nodes 0\n\
ic_uns_uniqify uns_sel_0\n\
ic_uns_subset_visible uns_sel_0 0\n\
ic_uns_create_selection_edgelist 0\n\
ic_uns_subset_delete uns_sel_0\n\
ic_undo_group_begin \n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_geo_set_part surface srf.00.S2 OUTLET 1\n\
ic_uns_update_family_type visible {{GEOM OUTLET ORFN SOLID}} {{!LINE_2 !QUAD_4 !HEXA_8}} update 0\n\
ic_delete_empty_parts \n\
ic_undo_group_end \n\
ic_undo_group_begin \n\
ic_geo_set_part surface srf.00.S1 INLET 1\n\
ic_uns_update_family_type visible {{INLET GEOM OUTLET ORFN SOLID}} {{!LINE_2 !QUAD_4 !HEXA_8}} update 0\n\
ic_delete_empty_parts \n\
ic_undo_group_end \n\
ic_uns_update_family_type visible {{GEOM OUTLET ORFN SOLID}} {{!LINE_2 !QUAD_4 !HEXA_8}} update 0\n\
ic_uns_update_family_type visible {{INLET GEOM OUTLET ORFN SOLID}} {{!LINE_2 !QUAD_4 !HEXA_8}} update 0\n\
ic_uns_update_family_type visible {{INLET OUTLET ORFN SOLID}} {{!LINE_2 !QUAD_4 !HEXA_8}} update 0\n\
ic_uns_update_family_type visible {{INLET GEOM OUTLET ORFN SOLID}} {{!LINE_2 !QUAD_4 !HEXA_8}} update 0\n\
ic_geo_rename_family GEOM WALLS 0\n\
ic_geo_rename_family GEOM WALLS 1\n\
ic_hex_write_file ./hex.uns WALLS SOLID OUTLET INLET proj 2 dim_to_mesh 3 no_boco\n\
ic_unload_mesh \n\
ic_delete_empty_parts \n\
ic_uns_load ./hex.uns 3 0 {{}} 1\n\
ic_uns_update_family_type visible {{INLET OUTLET ORFN WALLS SOLID}} {{!LINE_2 QUAD_4 !HEXA_8}} update 0\n\
ic_boco_solver \n\
ic_boco_clear_icons \n\
ic_boco_solver \n\
ic_boco_solver {{Ansys Fluent}}\n\
ic_solution_set_solver {{Ansys Fluent}} 1\n\
ic_boco_save ./ansys.fbc\n\
ic_boco_save_atr ./ansys.atr\n\
ic_delete_empty_parts \n\
ic_delete_empty_parts \n\
ic_save_tetin project1.tin 0 0 {{}} {{}} 0 0 1\n\
ic_uns_check_duplicate_numbers \n\
ic_save_unstruct project1.uns 1 {{}} {{}} {{}}\n\
ic_uns_set_modified 1\n\
ic_hex_save_blocking project1.blk\n\
ic_boco_solver \n\
ic_boco_solver {{Ansys Fluent}}\n\
ic_solution_set_solver {{Ansys Fluent}} 1\n\
ic_boco_save project1.fbc\n\
ic_boco_save_atr project1.atr\n\
ic_save_project_file ./project1.prj {{array\ set\ file_name\ \{{ {{    catia_dir .}} {{    parts_dir .}} {{    domain_loaded 0}} {{    cart_file_loaded 0}} {{    cart_file {{}}}} {{    domain_saved project1.uns}} {{    archive {{}}}} {{    med_replay {{}}}} {{    topology_dir .}} {{    ugparts_dir .}} {{    icons {{{{$env(ICEM_ACN)/lib/ai_env/icons}} {{$env(ICEM_ACN)/lib/va/EZCAD/icons}} {{$env(ICEM_ACN)/lib/icons}} {{$env(ICEM_ACN)/lib/va/CABIN/icons}}}}}} {{    tetin project1.tin}} {{    family_boco project1.fbc}} {{    iges_dir .}} {{    solver_params_loaded 0}} {{    attributes_loaded 0}} {{    project_lock {{}}}} {{    attributes project1.atr}} {{    domain project1.uns}} {{    domains_dir .}} {{    settings_loaded 0}} {{    settings project1.prj}} {{    blocking project1.blk}} {{    hexa_replay {{}}}} {{    transfer_dir .}} {{    mesh_dir .}} {{    family_topo {{}}}} {{    gemsparts_dir .}} {{    family_boco_loaded 0}} {{    tetin_loaded 0}} {{    project_dir .}} {{    topo_mulcad_out {{}}}} {{    solver_params {{}}}} \}} array\ set\ options\ \{{ {{    expert 1}} {{    remote_path {{}}}} {{    tree_disp_quad 2}} {{    tree_disp_pyra 0}} {{    evaluate_diagnostic 0}} {{    histo_show_default 1}} {{    select_toggle_corners 0}} {{    remove_all 0}} {{    keep_existing_file_names 0}} {{    record_journal 0}} {{    edit_wait 0}} {{    face_mode 1}} {{    select_mode all}} {{    med_save_emergency_tetin 1}} {{    user_name adix}} {{    diag_which all}} {{    uns_warn_if_display 500000}} {{    bubble_delay 1000}} {{    external_num 1}} {{    tree_disp_tri 2}} {{    apply_all 0}} {{    default_solver {{Ansys Fluent}}}} {{    temporary_directory {{}}}} {{    flood_select_angle 0}} {{    home_after_load 1}} {{    project_active 0}} {{    histo_color_by_quality_default 1}} {{    undo_logging 1}} {{    tree_disp_hexa 0}} {{    histo_solid_default 1}} {{    host_name brown-a016.rcac.purdue.edu}} {{    xhidden_full 1}} {{    replay_internal_editor 1}} {{    editor {{}}}} {{    mouse_color orange}} {{    clear_undo 1}} {{    remote_acn {{}}}} {{    remote_sh csh}} {{    tree_disp_penta 0}} {{    n_processors 1}} {{    remote_host {{}}}} {{    save_to_new 0}} {{    quality_info Quality}} {{    tree_disp_node 0}} {{    med_save_emergency_mesh 1}} {{    redtext_color red}} {{    tree_disp_line 0}} {{    select_edge_mode 0}} {{    use_dlremote 0}} {{    max_mesh_map_size 1024}} {{    show_tris 1}} {{    remote_user {{}}}} {{    enable_idle 0}} {{    auto_save_views 1}} {{    max_cad_map_size 512}} {{    display_origin 0}} {{    uns_warn_user_if_display 1000000}} {{    detail_info 0}} {{    win_java_help 0}} {{    show_factor 1}} {{    boundary_mode all}} {{    clean_up_tmp_files 1}} {{    auto_fix_uncovered_faces 1}} {{    med_save_emergency_blocking 1}} {{    max_binary_tetin 0}} {{    tree_disp_tetra 0}} \}} array\ set\ disp_options\ \{{ {{    uns_dualmesh 0}} {{    uns_warn_if_display 500000}} {{    uns_normals_colored 0}} {{    uns_icons 0}} {{    uns_locked_elements 0}} {{    uns_shrink_npos 0}} {{    uns_node_type None}} {{    uns_icons_normals_vol 0}} {{    uns_bcfield 0}} {{    backup Solid/wire}} {{    uns_nodes 0}} {{    uns_only_edges 0}} {{    uns_surf_bounds 0}} {{    uns_wide_lines 0}} {{    uns_vol_bounds 0}} {{    uns_displ_orient Triad}} {{    uns_orientation 0}} {{    uns_directions 0}} {{    uns_thickness 0}} {{    uns_shell_diagnostic 0}} {{    uns_normals 0}} {{    uns_couplings 0}} {{    uns_periodicity 0}} {{    uns_single_surfaces 0}} {{    uns_midside_nodes 1}} {{    uns_shrink 100}} {{    uns_multiple_surfaces 0}} {{    uns_no_inner 0}} {{    uns_enums 0}} {{    uns_disp Wire}} {{    uns_bcfield_name {{}}}} {{    uns_color_by_quality 0}} {{    uns_changes 0}} {{    uns_cut_delay_count 1000}} \}} {{set icon_size1 24}} {{set icon_size2 35}} {{set thickness_defined 0}} {{set solver_type 1}} {{set solver_setup -1}} array\ set\ prism_values\ \{{ {{    n_triangle_smoothing_steps 5}} {{    min_smoothing_steps 6}} {{    first_layer_smoothing_steps 1}} {{    new_volume {{}}}} {{    height {{}}}} {{    prism_height_limit {{}}}} {{    interpolate_heights 0}} {{    n_tetra_smoothing_steps 10}} {{    do_checks {{}}}} {{    delete_standalone 1}} {{    ortho_weight 0.50}} {{    max_aspect_ratio {{}}}} {{    ratio_max {{}}}} {{    incremental_write 0}} {{    total_height {{}}}} {{    use_prism_v10 0}} {{    intermediate_write 1}} {{    delete_base_triangles {{}}}} {{    ratio_multiplier {{}}}} {{    verbosity_level 1}} {{    refine_prism_boundary 1}} {{    max_size_ratio {{}}}} {{    triangle_quality {{}}}} {{    max_prism_angle 180}} {{    tetra_smooth_limit 0.3}} {{    max_jump_factor 5}} {{    use_existing_quad_layers 0}} {{    layers 3}} {{    fillet 0.10}} {{    into_orphan 0}} {{    init_dir_from_prev {{}}}} {{    blayer_2d 0}} {{    do_not_allow_sticking {{}}}} {{    top_family {{}}}} {{    law exponential}} {{    min_smoothing_val 0.1}} {{    auto_reduction 0}} {{    stop_columns 1}} {{    stair_step 1}} {{    smoothing_steps 12}} {{    side_family {{}}}} {{    min_prism_quality 0.01}} {{    ratio 1.2}} \}} {{set aie_current_flavor {{}}}} array\ set\ vid_options\ \{{ {{    auxiliary 0}} {{    show_name 0}} {{    inherit 1}} {{    default_part GEOM}} {{    new_srf_topo 1}} {{    DelPerFlag 0}} {{    show_item_name 0}} {{    composite_tolerance 1.0}} {{    replace 0}} {{    same_pnt_tol 1e-4}} {{    tdv_axes 1}} {{    vid_mode 0}} {{    DelBlkPerFlag 0}} \}} {{set savedTreeVisibility {{geomNode 1 geom_subsetNode 0 geomPointNode 0 geomCurveNode 0 geomSurfNode 2 meshNode 1 mesh_subsetNode 0 meshLineNode 0 meshShellNode 2 meshQuadNode 2 meshVolumeNode 0 meshHexaNode 0 blockingNode 1 block_subsetNode 2 block_vertNode 0 block_edgeNode 0 block_faceNode 0 block_blockNode 0 block_meshNode 0 topoNode 2 topo-root 2 partNode 2 part-INLET 2 part-OUTLET 2 part-SOLID 2 part-VORFN 0 part-WALLS 2}}}} {{set last_view {{rot {{3.19905930365e-05 0.999933204826 0.00429823005183 0.0107289367044}} scale {{12258.8222189 12258.8222189 12258.8222189}} center {{0 0 1.2700001001400001}} pos {{-93.0 19.0 0}}}}}} array\ set\ cut_info\ \{{ {{    active 0}} {{    whole 1}} \}} array\ set\ hex_option\ \{{ {{    default_bunching_ratio 2.0}} {{    floating_grid 0}} {{    project_to_topo 0}} {{    n_tetra_smoothing_steps 20}} {{    sketching_mode 0}} {{    trfDeg 1}} {{    wr_hexa7 0}} {{    hexa_projection_mode 0}} {{    smooth_ogrid 0}} {{    find_worst 1-3}} {{    hexa_verbose_mode 0}} {{    old_eparams 0}} {{    uns_face_mesh_method uniform_quad}} {{    multigrid_level 0}} {{    uns_face_mesh one_tri}} {{    check_blck 0}} {{    proj_limit 0}} {{    check_inv 0}} {{    project_bspline 0}} {{    hexa_update_mode 1}} {{    default_bunching_law BiGeometric}} {{    worse_criterion Quality}} \}} array\ set\ saved_views\ \{{ {{    views {{}}}} \}}}} {{ICEM CFD}}\n\
ic_exec {fluent_translator_path} -dom ./hex.uns -b project1.fbc ./{case_name}.msh\n\
ic_uns_num_couplings \n\
ic_undo_group_begin \n\
ic_uns_create_diagnostic_edgelist 1 \n\
ic_uns_diagnostic subset all diag_type uncovered fix_fam FIX_UNCOVERED diag_verb {{Uncovered faces}} fams {{}} busy_off 1 quiet 1 \n\
ic_uns_create_diagnostic_edgelist 0 \n\
ic_undo_group_end \n\
ic_uns_min_metric Quality {{}} {{}} \n\
ic_exit\n', file=fi)
    try:
        subprocess.check_call("ml ansys", shell=True)
    except subprocess.CalledProcessError as e:
        print(e)
        print("Continuing...")
    comp_process = subprocess.call('icemcfd -x ./mesh_replay.rpl > auto_cfx_run.log', shell=True)
    print(comp_process)
    if cleanup:
        subprocess.check_call('rm mesh_replay.rpl', shell=True)
        subprocess.check_call('rm project1.*', shell=True)
        subprocess.check_call('rm ansys.*', shell=True)
        subprocess.check_call('rm hex.uns', shell=True)
    return


def calculate_CL_Ryan(jf, jg, theta=0.0, Db = 0.0018, testing = False):

    """ Calculates the lift coefficient based on Ryan et al. (2023) method (i.e. DFM at various angles)

    Also applies Sharma correction
    
    """
    j = jg + jf

    if theta == 90:
        C0 = 1.05
        Vgj = 0.112
    elif theta == 80:
        C0 = 1.042
        Vgj = -0.134
    elif theta == 60:
        C0 = 1.008
        Vgj = -0.19
    elif theta == 30:
        C0 = 0.988
        Vgj = -0.159
    elif theta == 0:
        C0 = 0.986
        Vgj = -0.248
    else:
        print('Did not recognize %f, defaulting to vertical' % theta)
        C0 = 1.05
        Vgj = 0.112
    
    alpha90 = jg / (1.05*j + 0.112)
    S90 = (jg / alpha90) / (jf / (1-alpha90))

    alpha = jg / (C0*j + Vgj)
    
    S = (jg / alpha) / (jf / (1-alpha) )

    I = (S - 1) / (S90 - 1)
    #print(I, theta, alpha, alpha90)

    if testing:
        return I

    ### Tomiyama lift coefficient ###
    vr = jf / (1 - alpha) * (S-1)

    g = 9.81
    # Hardcoded air/water properties
    rhof = 998
    rhog = 1.225
    sigma = 0.072
    muf = 0.001

    Eo = g * (rhof - rhog) * Db**2 / sigma
    Reb = rhof * Db * np.abs(vr) / muf

    f = 0.00105*Eo**3-0.0159*Eo**2-0.0204*Eo+0.474
    CL_T = min(0.288 * np.tanh(0.1221 * Reb), f)

    ### Sharma correction ###
    Ref = rhof * jf * 0.0254 / muf # Hardcoded diameter!!!
    f_sharma = np.exp( (max(Ref, 5e4) - 5e4) / -1.22e5 )

    CL = CL_T * I * f_sharma
    return CL

def write_CCL(mom_source = 'normal_drag_mom_source', ccl_name = 'auto_setup.ccl',
                          inDataFile= '/home/adix/CFD/exp_BCs/in_H_0.csv', 
                          outDataFile = '/home/adix/CFD/exp_BCs/out_H_0.csv', 
                          Db=0.0018, CL = 0.25, CTD = 0.25, theta = 0, mdot = 2.02277,
                          Kf = 0.083, Kw = 0.98, CD = 0.44, jf = 4, jg = 0.11, p_out = 20,
                          num_iter = 10000, resid_target = 1e-6):
    """ Function to write the CCL file that CFX reads

    Can use different momentum source terms, current options are
    - "drag_model", drag force calculated with 
    - "IS", 
    - "normal_drag_mom_source", regular drag formulation

    Other important simulation parameters
    - Db, bubble diameter
    - mdot, mass flow rate for BC
    - CL, lift coefficient. Can set to a constant, 'tomiyama', or 'ryan_DFM'. For the latter cases, also need jf and jg. They both also use Sharma's lift modification
    - CTD, turbulent dispersion coefficient. Using Lopez-de-Bertodano model
    - CD, drag coefficient
    - Kf, Kw, for relative velocity model
    - jf, jg, needed for lift coefficient models
    - num_iter, number of iterations
    - resid_target, residual convergence criteria
    
    """
    
    with open(inDataFile, 'r') as infi:
        trash = infi.readline()
        InletData = infi.readline().split(',')[0] # Name of function in csv

    with open(outDataFile, 'r') as infi:
        trash = infi.readline()
        OutletData = infi.readline().split(',')[0] # Name of function in csv
    CD_CFX = 0
    if mom_source == 'IS':
        mom_source_to_write = ['ISfx', 'ISfy', 'ISfz', 'ISgx', 'ISgy', 'ISgz']
    elif mom_source == 'drag_model':
        mom_source_to_write = ['FDLx', 'FDLy', 'FDLz', 'FDGx', 'FDGy', 'FDGz']
    elif mom_source == 'normal_drag_mom_source':
        mom_source_to_write = ['FDLx', 'FDLy', 'FDLz', 'FDGx', 'FDGy', 'FDGz']
        Kf = 0
        Kw = 0
    else:
        mom_source_to_write = ['0 [kg m^-2 s^-2]', '0 [kg m^-2 s^-2]', '0 [kg m^-2 s^-2]', '0 [kg m^-2 s^-2]', '0 [kg m^-2 s^-2]', '0 [kg m^-2 s^-2]']
        CD_CFX = CD

    if CL == 'tomiyama':
        lift_string = f"Option = Tomiyama \n"
    elif CL == 'legendre':
        lift_string = f"Option = Legendre \n"
    else:

        if CL == 'ryan_DFM':
            CL = calculate_CL_Ryan(jf=jf, jg=jg, theta=theta, Db=Db)
        elif CL == 'tomi-sharma':
            CL = 'CL'

        lift_string = f"Lift Coefficient = {CL} \nOption = Lift Coefficient \n"

    with open(ccl_name, 'w') as fi:

        strToWrite = f"\
\
LIBRARY: \n\
CEL: \n\
EXPRESSIONS: \n\
CD = 0.44 \n\
FDGx = - 3/4 * CD / (gas.Mean Particle Diameter) * (gas.Volume Fraction) * liquid.density * (vrNorm) * (gas.u - liquid.u) \n\
FDGy = - 3/4 * CD / (gas.Mean Particle Diameter) * (gas.Volume Fraction) * liquid.density* (vrNorm) * (gas.v - liquid.v) \n\
FDGz = - 3/4 * CD / (gas.Mean Particle Diameter) * (gas.Volume Fraction) * liquid.density* (vrNorm) * (gas.w - liquidWEff) \n\
FDLx = 3/4 * CD / (gas.Mean Particle Diameter) * (gas.Volume Fraction) * liquid.density * (vrNorm) * (gas.u - liquid.u) \n\
FDLy = 3/4 * CD / (gas.Mean Particle Diameter) * (gas.Volume Fraction) * liquid.density * (vrNorm) * (gas.v - liquid.v) \n\
FDLz = 3/4 * CD / (gas.Mean Particle Diameter) * (gas.Volume Fraction) * liquid.density * (vrNorm) * (gas.w- liquidWEff) \n\
Kf = {Kf} \n\
Kw = {Kw} \n\
facilitytheta = {theta} \n\
gravy = -9.81 [m s^-2]*cos(facilitytheta* pi / 180) \n\
gravz = -9.81 [m s^-2]*sin(facilitytheta* pi / 180) \n\
liquidWEff = (1 - Kf - Kw * gas.Volume Fraction * CD^(1/3) ) * liquid.w \n\
phi = atan2(y, x) \n\
radius = sqrt(x^2 + y^2) \n\
vrNorm = sqrt( (gas.u - liquid.u)^2+ (gas.v - liquid.v)^2+ (gas.w - liquidWEff)^2 ) \n\
Eo = gravz * (liquid.density - gas.density) * ( gas.Mean Particle Diameter )^2 / (gas | liquid.Surface Tension Coefficient) \n\
Eop = gravz * (liquid.density - gas.density) * ( gas.Mean Particle Diameter * (1+0.136*Eo^0.757)^(1/3.) )^2 / (gas | liquid.Surface Tension Coefficient) \n\
Rep = liquid.density * abs(gas.w - liquidWEff) * gas.Mean Particle Diameter / liquid.viscosity \n\
f = 0.00105*Eop^3-0.0159*Eop^2-0.0204*Eop + 0.474 \n\
CLtomiyama = min( 0.288*tanh(0.121*Rep) , f )\n\
Ref = liquid.density * areaAve(liquid.w)@inlet * 0.0254 [m] / liquid.viscosity\n\
sharmaFactor = exp( (max(Ref, 5e4) - 5e4) / (-1.22e5) )\n\
CL = CLtomiyama * sharmaFactor\n\
END \n\
FUNCTION: {InletData} \n\
Argument Units = [mm], [m], [] \n\
Option = Profile Data \n\
Reference Coord Frame = Coord 0 \n\
Spatial Fields = radius, z, phi \n\
DATA FIELD: Velocity u g\n\
Field Name = Velocity u g\n\
Parameter List = U,Velocity r Component,Wall U,Wall Velocity r Component \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity u f \n\
Field Name = Velocity u f \n\
Parameter List = U,Velocity r Component,Wall U,Wall Velocity r Component \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity v g\n\
Field Name = Velocity v g\n\
Parameter List = V,Velocity Theta Component,Wall V,Wall Velocity Theta Component \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity v f \n\
Field Name = Velocity v f \n\
Parameter List = V,Velocity Theta Component,Wall V,Wall Velocity Theta Component \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity w g\n\
Field Name = Velocity w g\n\
Parameter List = Velocity Axial Component,W,Wall Velocity Axial Component,Wall W \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity w f \n\
Field Name = Velocity w f \n\
Parameter List = Velocity Axial Component,W,Wall Velocity Axial Component,Wall W \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Volume Fraction \n\
Field Name = Volume Fraction \n\
Parameter List = Volume Fraction \n\
Result Units = [] \n\
END \n\
DATA SOURCE: \n\
File Name = {inDataFile} \n\
Option = From File \n\
END \n\
END \n\
FUNCTION: {OutletData} \n\
Argument Units = [mm], [m], [] \n\
Option = Profile Data \n\
Reference Coord Frame = Coord 0 \n\
Spatial Fields = radius, z, phi \n\
DATA FIELD: Velocity u g\n\
Field Name = Velocity u g\n\
Parameter List = U,Velocity r Component,Wall U,Wall Velocity r Component \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity u f \n\
Field Name = Velocity u f \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity v g\n\
Field Name = Velocity v g\n\
Parameter List = V,Velocity Theta Component,Wall V,Wall Velocity Theta Component \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity v f \n\
Field Name = Velocity v f \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity w g\n\
Field Name = Velocity w g\n\
Parameter List = Velocity Axial Component,W,Wall Velocity Axial Component,Wall W \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Velocity w f \n\
Field Name = Velocity w f \n\
Result Units = [m s^-1] \n\
END \n\
DATA FIELD: Volume Fraction \n\
Field Name = Volume Fraction \n\
Parameter List = Volume Fraction \n\
Result Units = [] \n\
END \n\
DATA SOURCE: \n\
File Name = {outDataFile} \n\
Option = From File \n\
END \n\
END \n\
END \n\
MATERIAL: Air Ideal Gas \n\
Material Description = Air Ideal Gas (constant Cp) \n\
Material Group = Air Data, Calorically Perfect Ideal Gases \n\
Option = Pure Substance \n\
Thermodynamic State = Gas \n\
PROPERTIES: \n\
Option = General Material \n\
EQUATION OF STATE: \n\
Molar Mass = 28.96 [kg kmol^-1] \n\
Option = Ideal Gas \n\
END \n\
SPECIFIC HEAT CAPACITY: \n\
Option = Value \n\
Specific Heat Capacity = 1.0044E+03 [J kg^-1 K^-1] \n\
Specific Heat Type = Constant Pressure \n\
END \n\
REFERENCE STATE: \n\
Option = Specified Point \n\
Reference Pressure = 1 [atm] \n\
Reference Specific Enthalpy = 0. [J/kg] \n\
Reference Specific Entropy = 0. [J/kg/K] \n\
Reference Temperature = 25 [C] \n\
END \n\
DYNAMIC VISCOSITY: \n\
Dynamic Viscosity = 1.831E-05 [kg m^-1 s^-1] \n\
Option = Value \n\
END \n\
THERMAL CONDUCTIVITY: \n\
Option = Value \n\
Thermal Conductivity = 2.61E-2 [W m^-1 K^-1] \n\
END \n\
ABSORPTION COEFFICIENT: \n\
Absorption Coefficient = 0.01 [m^-1] \n\
Option = Value \n\
END \n\
SCATTERING COEFFICIENT: \n\
Option = Value \n\
Scattering Coefficient = 0.0 [m^-1] \n\
END \n\
REFRACTIVE INDEX: \n\
Option = Value \n\
Refractive Index = 1.0 [m m^-1] \n\
END \n\
END \n\
END \n\
MATERIAL: Water \n\
Material Description = Water (liquid) \n\
Material Group = Water Data, Constant Property Liquids \n\
Option = Pure Substance \n\
Thermodynamic State = Liquid \n\
PROPERTIES: \n\
Option = General Material \n\
EQUATION OF STATE: \n\
Density = 997.0 [kg m^-3] \n\
Molar Mass = 18.02 [kg kmol^-1] \n\
Option = Value \n\
END \n\
SPECIFIC HEAT CAPACITY: \n\
Option = Value \n\
Specific Heat Capacity = 4181.7 [J kg^-1 K^-1] \n\
Specific Heat Type = Constant Pressure \n\
END \n\
REFERENCE STATE: \n\
Option = Specified Point \n\
Reference Pressure = 1 [atm] \n\
Reference Specific Enthalpy = 0.0 [J/kg] \n\
Reference Specific Entropy = 0.0 [J/kg/K] \n\
Reference Temperature = 25 [C] \n\
END \n\
DYNAMIC VISCOSITY: \n\
Dynamic Viscosity = 8.899E-4 [kg m^-1 s^-1] \n\
Option = Value \n\
END \n\
THERMAL CONDUCTIVITY: \n\
Option = Value \n\
Thermal Conductivity = 0.6069 [W m^-1 K^-1] \n\
END \n\
ABSORPTION COEFFICIENT: \n\
Absorption Coefficient = 1.0 [m^-1] \n\
Option = Value \n\
END \n\
SCATTERING COEFFICIENT: \n\
Option = Value \n\
Scattering Coefficient = 0.0 [m^-1] \n\
END \n\
REFRACTIVE INDEX: \n\
Option = Value \n\
Refractive Index = 1.0 [m m^-1] \n\
END \n\
THERMAL EXPANSIVITY: \n\
Option = Value \n\
Thermal Expansivity = 2.57E-04 [K^-1] \n\
END \n\
END \n\
END \n\
END \n\
FLOW: Flow Analysis 1 \n\
SOLUTION UNITS: \n\
Angle Units = [rad] \n\
Length Units = [m] \n\
Mass Units = [kg] \n\
Solid Angle Units = [sr] \n\
Temperature Units = [K] \n\
Time Units = [s] \n\
END \n\
ANALYSIS TYPE: \n\
Option = Steady State \n\
EXTERNAL SOLVER COUPLING: \n\
Option = None \n\
END \n\
END \n\
DOMAIN: Default Domain \n\
Coord Frame = Coord 0 \n\
Domain Type = Fluid \n\
Location = SOLID \n\
BOUNDARY: inlet \n\
Boundary Type = INLET \n\
Location = INLET \n\
Use Profile Data = False \n\
BOUNDARY CONDITIONS: \n\
FLOW REGIME: \n\
Option = Subsonic \n\
END \n\
MASS AND MOMENTUM: \n\
Option = Fluid Velocity \n\
END \n\
TURBULENCE: \n\
Option = Medium Intensity and Eddy Viscosity Ratio \n\
END \n\
END \n\
FLUID: gas \n\
BOUNDARY CONDITIONS: \n\
VELOCITY: \n\
Option = Cartesian Velocity Components \n\
U = 0 [m s^-1] \n\
V = 0 [m s^-1] \n\
W = {InletData}.Velocity w g(radius,z,phi) \n\
END \n\
VOLUME FRACTION: \n\
Option = Value \n\
Volume Fraction = {InletData}.Volume Fraction(radius,z,phi) \n\
END \n\
END \n\
END \n\
FLUID: liquid \n\
BOUNDARY CONDITIONS: \n\
FLOW DIRECTION: \n\
Option = Normal to Boundary Condition \n\
END \n\
VELOCITY: \n\
Mass Flow Rate = {mdot} [kg s^-1] \n\
Option = Mass Flow Rate \n\
END \n\
VOLUME FRACTION: \n\
Option = Value \n\
Volume Fraction = 1-{InletData}.Volume Fraction(radius,z,phi) \n\
END \n\
END \n\
END \n\
END \n\
BOUNDARY: outlet\n\
Boundary Type = OPENING\n\
Location = OUTLET\n\
Use Profile Data = False\n\
BOUNDARY CONDITIONS:\n\
FLOW DIRECTION:\n\
Option = Normal to Boundary Condition\n\
END\n\
FLOW REGIME:\n\
Option = Subsonic\n\
END\n\
MASS AND MOMENTUM:\n\
Option = Opening Pressure and Direction\n\
Relative Pressure = {p_out} [psi]\n\
END\n\
TURBULENCE:\n\
Option = Medium Intensity and Eddy Viscosity Ratio\n\
END\n\
END\n\
FLUID: gas\n\
BOUNDARY CONDITIONS:\n\
VOLUME FRACTION:\n\
Option = Zero Gradient\n\
END\n\
END\n\
END\n\
FLUID: liquid\n\
BOUNDARY CONDITIONS:\n\
VOLUME FRACTION:\n\
Option = Zero Gradient\n\
END\n\
END\n\
END\n\
END\n\
BOUNDARY: walls \n\
Boundary Type = WALL \n\
Location = WALLS \n\
Use Profile Data = False \n\
BOUNDARY CONDITIONS: \n\
MASS AND MOMENTUM: \n\
Option = Fluid Dependent \n\
END \n\
WALL CONTACT MODEL: \n\
Option = Use Volume Fraction \n\
END \n\
WALL ROUGHNESS: \n\
Option = Smooth Wall \n\
END \n\
END \n\
FLUID: gas \n\
BOUNDARY CONDITIONS: \n\
MASS AND MOMENTUM: \n\
Option = No Slip Wall \n\
END \n\
END \n\
END \n\
FLUID: liquid \n\
BOUNDARY CONDITIONS: \n\
MASS AND MOMENTUM: \n\
Option = No Slip Wall \n\
END \n\
END \n\
END \n\
END \n\
DOMAIN MODELS: \n\
BUOYANCY MODEL: \n\
Buoyancy Reference Density = 998 [kg m^-3] \n\
Gravity X Component = 0 [m s^-2] \n\
Gravity Y Component = gravy \n\
Gravity Z Component = gravz \n\
Option = Buoyant \n\
BUOYANCY REFERENCE LOCATION: \n\
Option = Automatic \n\
END \n\
END \n\
DOMAIN MOTION: \n\
Option = Stationary \n\
END \n\
MESH DEFORMATION: \n\
Option = None \n\
END \n\
REFERENCE PRESSURE: \n\
Reference Pressure = 1 [atm] \n\
END \n\
END \n\
FLUID DEFINITION: gas \n\
Material = Air Ideal Gas \n\
Option = Material Library \n\
MORPHOLOGY: \n\
Mean Diameter = {Db} [m] \n\
Minimum Volume Fraction = 0.000001 \n\
Option = Dispersed Fluid \n\
END \n\
END \n\
FLUID DEFINITION: liquid \n\
Material = Water \n\
Option = Material Library \n\
MORPHOLOGY: \n\
Minimum Volume Fraction = 0.000001 \n\
Option = Continuous Fluid \n\
END \n\
END \n\
FLUID MODELS: \n\
COMBUSTION MODEL: \n\
Option = None \n\
END \n\
FLUID: gas \n\
FLUID BUOYANCY MODEL: \n\
Option = Density Difference \n\
END \n\
TURBULENCE MODEL: \n\
Option = Dispersed Phase Zero Equation \n\
END \n\
END \n\
FLUID: liquid \n\
FLUID BUOYANCY MODEL: \n\
Option = Density Difference \n\
END \n\
TURBULENCE MODEL: \n\
Option = k epsilon \n\
BUOYANCY TURBULENCE: \n\
Option = None \n\
END \n\
END \n\
TURBULENT WALL FUNCTIONS: \n\
Option = Scalable \n\
END \n\
END \n\
HEAT TRANSFER MODEL: \n\
Fluid Temperature = 25 [C] \n\
Homogeneous Model = False \n\
Option = Isothermal \n\
END \n\
THERMAL RADIATION MODEL: \n\
Option = None \n\
END \n\
TURBULENCE MODEL: \n\
Homogeneous Model = False \n\
Option = Fluid Dependent \n\
END \n\
END \n\
FLUID PAIR: gas | liquid \n\
Surface Tension Coefficient = 0.072 [N m^-1] \n\
INTERPHASE TRANSFER MODEL: \n\
Minimum Volume Fraction for Area Density = 0.000001 \n\
Option = Particle Model \n\
END \n\
MASS TRANSFER: \n\
Option = None \n\
END \n\
MOMENTUM TRANSFER: \n\
DRAG FORCE: \n\
Drag Coefficient = {CD_CFX} \n\
Option = Drag Coefficient \n\
END \n\
LIFT FORCE: \n\
{lift_string}\
END \n\
TURBULENT DISPERSION FORCE: \n\
Option = Lopez de Bertodano \n\
Turbulent Dispersion Coefficient = {CTD} \n\
END \n\
VIRTUAL MASS FORCE: \n\
Option = None \n\
END \n\
WALL LUBRICATION FORCE: \n\
Lubrication Coefficient C1 = -0.01 \n\
Lubrication Coefficient C2 = 0.05 \n\
Option = Antal \n\
END \n\
END \n\
TURBULENCE TRANSFER: \n\
ENHANCED TURBULENCE PRODUCTION MODEL: \n\
Option = Sato Enhanced Eddy Viscosity \n\
END \n\
END \n\
END \n\
INITIALISATION:\n\
Option = Automatic\n\
FLUID: gas\n\
INITIAL CONDITIONS:\n\
Velocity Type = Cartesian\n\
CARTESIAN VELOCITY COMPONENTS:\n\
Option = Automatic with Value\n\
U = 0 [m s^-1]\n\
V = 0 [m s^-1]\n\
W = 3.5 [m s^-1]\n\
END\n\
VOLUME FRACTION:\n\
Option = Automatic\n\
END\n\
END\n\
END\n\
FLUID: liquid\n\
INITIAL CONDITIONS:\n\
Velocity Type = Cartesian\n\
CARTESIAN VELOCITY COMPONENTS:\n\
Option = Automatic with Value\n\
U = 0 [m s^-1]\n\
V = 0 [m s^-1]\n\
W = 4 [m s^-1]\n\
END\n\
TURBULENCE INITIAL CONDITIONS:\n\
Option = Medium Intensity and Eddy Viscosity Ratio\n\
END\n\
VOLUME FRACTION:\n\
Option = Automatic\n\
END\n\
END\n\
END\n\
INITIAL CONDITIONS: \n\
STATIC PRESSURE: \n\
Option = Automatic \n\
END \n\
END \n\
END \n\
MULTIPHASE MODELS: \n\
Homogeneous Model = False \n\
FREE SURFACE MODEL: \n\
Option = None \n\
END \n\
END \n\
SUBDOMAIN: fluid \n\
Coord Frame = Coord 0 \n\
Location = SOLID \n\
FLUID: gas \n\
SOURCES: \n\
MOMENTUM SOURCE: \n\
GENERAL MOMENTUM SOURCE: \n\
Momentum Source Coefficient = 1000000 [kg m^-3 s^-1] \n\
Momentum Source X Component = {mom_source_to_write[3]} \n\
Momentum Source Y Component = {mom_source_to_write[4]} \n\
Momentum Source Z Component = {mom_source_to_write[5]} \n\
Option = Cartesian Components \n\
END \n\
END \n\
END \n\
END \n\
FLUID: liquid \n\
SOURCES: \n\
MOMENTUM SOURCE: \n\
GENERAL MOMENTUM SOURCE: \n\
Momentum Source Coefficient = 1000000 [kg m^-3 s^-1] \n\
Momentum Source X Component = {mom_source_to_write[0]} \n\
Momentum Source Y Component = {mom_source_to_write[1]} \n\
Momentum Source Z Component = {mom_source_to_write[2]} \n\
Option = Cartesian Components \n\
END \n\
END \n\
END \n\
END \n\
END \n\
END \n\
OUTPUT CONTROL: \n\
RESULTS: \n\
File Compression Level = Default \n\
Option = Standard \n\
END \n\
END \n\
SOLVER CONTROL: \n\
Optimization Level = 3 \n\
Turbulence Numerics = High Resolution \n\
ADVECTION SCHEME: \n\
Option = High Resolution \n\
END \n\
CONVERGENCE CONTROL: \n\
Length Scale Option = Conservative \n\
Maximum Number of Iterations = {num_iter} \n\
Minimum Number of Iterations = 1 \n\
Timescale Control = Auto Timescale \n\
Timescale Factor = 0.5 \n\
END \n\
CONVERGENCE CRITERIA: \n\
Conservation Target = 0.01 \n\
Residual Target = {resid_target} \n\
Residual Type = RMS \n\
END \n\
DYNAMIC MODEL CONTROL: \n\
Global Dynamic Model Control = On \n\
END \n\
INTERRUPT CONTROL: \n\
Option = Any Interrupt \n\
CONVERGENCE CONDITIONS: \n\
Option = Default Conditions \n\
END \n\
END \n\
MULTIPHASE CONTROL: \n\
Initial Volume Fraction Smoothing = Volume-Weighted \n\
Optimization Level = 3 \n\
Volume Fraction Coupling = Coupled \n\
END \n\
END \n\
EXPERT PARAMETERS: \n\
diffusion coef averaging type = 3 \n\
stress coef averaging type = 3 \n\
END \n\
END \n\
COMMAND FILE: \n\
Version = 22.2 \n\
END"

        print(strToWrite, file = fi)
        return

def make_CFX_case(case_name, mesh_name = None, ccl_name='auto_setup'):
    """ Makes CFX case based on CCL file

    Takes "setup_ccl" and writes a case file, "casename.def"
    
    """

    if mesh_name is None:
        mesh_name = case_name

    with open('CFXPre_Commands.pre', 'w') as fi:
        print(f"\
# CFX-Pre commands\n\
\n\
COMMAND FILE:\n\
  CFX Pre Version = 22.1\n\
END\n\
\n\
>load mode=new\n\
> update\n\
\n\
> gtmImport filename={mesh_name}.msh, type=Fluent, units=m, genOpt= -n, nameStrategy=Assembly\n\
> update\n\
\n\
>importccl filename={ccl_name}.ccl, mode=replace, autoLoadLibrary=none\n\
> update\n\
\n\
>writeCaseFile filename=./{case_name}.def, operation=write def file\n\
> update\n\
\n\
> quit ", file = fi)
    print('$\ncfx5pre -s CFXPre_Commands.pre')
    
    subprocess.run("ml ansys", shell=True)
    
    comp_process = subprocess.check_call('cfx5pre -s CFXPre_Commands.pre -line > auto_cfx_run.log', shell=True)
    print(comp_process)
    return

def run_CFX_case(case_name, parallel=True, npart = 4, init_fi = None, interactive = False):
    """ Runs CFX case case_name.def

    Must be in the same directory, or a full path (without the .def extension) supplied

    Can run in parallel, specify the number of cores with npart
    
    """
    comp_process = subprocess.run("ml ansys", shell=True)
    print(comp_process)

    run_string = f"cfx5solve -def {case_name}.def -monitor {case_name}_001.out -double"
    
    if parallel:
        run_string += f" -par -par-local -part {npart}"

    if init_fi:
        run_string += f" -initial-file {init_fi}"

    if interactive:
        run_string += f" -interactive"
    
    run_string += " > auto_cfx_run.log"

    print(f"${run_string}")
    comp_process = subprocess.run(run_string, shell=True)
    if comp_process.returncode != 0:
        print(comp_process)
    return


def post_process_CFX_case(case_name, template = 'CFX'):
    """ Post-process a CFX case, case_name.res

    Must have a template available, "CFX" is the default. Might need to set it up beforehand based on a .cst file

    TODO update to export csv file that can be read into MARIGOLD Condition
    
    """

    with open("CFXPost_Commands.cse", 'w') as fi:
        print(f'\
\
COMMAND FILE:\n\
  CFX Post Version = 22.1\n\
END\n\
\n\
DATA READER:\n\
  Clear All Objects = false\n\
  Append Results = false\n\
  Edit Case Names = false\n\
  Multi Configuration File Load Option = Last Case\n\
  Open in New View = true\n\
  Keep Camera Position = true\n\
  Load Particle Tracks = true\n\
  Multi Configuration File Load Option = Last Case\n\
  Construct Variables From Fourier Coefficients = true\n\
  Open to Compare = false\n\
  Files to Compare =\n\
END\n\
\n\
DATA READER:\n\
  Domains to Load=\n\
END\n\
\n\
>load filename={case_name}_001.res, force_reload=true\n\
\n\
VIEW:View 1\n\
  Camera Mode = User Specified\n\
  CAMERA:\n\
    Option = Pivot Point and Quaternion\n\
    Pivot Point = 0, 0, 1.27\n\
    Scale = 1.70237\n\
    Pan = 0, 0\n\
    Rotation Quaternion = 0.279848, -0.364705, -0.115917, 0.880476\n\
  END\n\
\n\
END\n\
\n\
> update\n\
>report generatemode=APPEND, loadtemplate={template}\n\
\n\
REPORT:\n\
  PUBLISH:\n\
    Generate 3D Viewer Files = Off\n\
    Report Format = HTML\n\
    Report Path = {case_name}_Results.html\n\
    Save Images In Separate Folder = On\n\
\n\
  END\n\
END\n\
>report save\n\
EXPORT:\n\
  ANSYS Export Data = Element Heat Flux\n\
  ANSYS File Format = ANSYS\n\
  ANSYS Reference Temperature = 0.0 [K]\n\
  ANSYS Specify Reference Temperature = Off\n\
  ANSYS Supplemental HTC = 0.0 [W m^-2 K^-1]\n\
  Additional Variable List =\n\
  BC Profile Type = Inlet Velocity\n\
  CSV Type = CSV\n\
  Case Name = Case {case_name}_001\n\
  Export Connectivity = Off\n\
  Export Coord Frame = Global\n\
  Export File = {case_name}_port3_phi90.csv\n\
  Export Geometry = Off\n\
  Export Location Aliases =\n\
  Export Node Numbers = Off\n\
  Export Null Data = On\n\
  Export Type = Generic\n\
  Export Units System = Current\n\
  Export Variable Type = Current\n\
  External Export Data = None\n\
  Include File Information = Off\n\
  Include Header = On\n\
  Location = inlet\n\
  Location List = /LINE:port 3 phi 90 rake\n\
  Null Token = null\n\
  Overwrite = On\n\
  Precision = 8\n\
  Separator = ", "\n\
  Spatial Variables = X,Y,Z\n\
  Variable List = X, Y, Z, gas.Velocity w, gas.Volume Fraction, liquid.Velocity w\n\
  Vector Brackets = ()\n\
  Vector Display = Scalar\n\
END\n\
>export\n\
EXPORT:\n\
  ANSYS Export Data = Element Heat Flux\n\
  ANSYS File Format = ANSYS\n\
  ANSYS Reference Temperature = 0.0 [K]\n\
  ANSYS Specify Reference Temperature = Off\n\
  ANSYS Supplemental HTC = 0.0 [W m^-2 K^-1]\n\
  Additional Variable List =\n\
  BC Profile Type = Inlet Velocity\n\
  CSV Type = CSV\n\
  Case Name = Case {case_name}_001\n\
  Export Connectivity = Off\n\
  Export Coord Frame = Global\n\
  Export File = {case_name}_port3.csv\n\
  Export Geometry = Off\n\
  Export Location Aliases =\n\
  Export Node Numbers = Off\n\
  Export Null Data = On\n\
  Export Type = Generic\n\
  Export Units System = Current\n\
  Export Variable Type = Current\n\
  External Export Data = None\n\
  Include File Information = Off\n\
  Include Header = On\n\
  Location = inlet\n\
  Location List = /PLANE:port3\n\
  Null Token = null\n\
  Overwrite = On\n\
  Precision = 8\n\
  Separator = ", "\n\
  Spatial Variables = X,Y,Z \n\
  Variable List = gas.Velocity w, gas.Volume Fraction, liquid.Velocity w, X, Y, Z\n\
  Vector Brackets = ()\n\
  Vector Display = Scalar\n\
END\n\
>export\n\
\n\
\n\
>quit     \n\
        \
        ', file=fi)

    subprocess.run("ml ansys", shell=True)

    print(os.getcwd())
    subprocess.check_call(f'rm -rf ./{case_name}_Results', shell=True)
    comp_process = subprocess.run(f'cfx5post -play CFXPost_Commands.cse -line > auto_cfx_run.log', shell=True)
    print(comp_process)




