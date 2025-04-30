from .config import *

def color_cycle(set_color = None, color_list = []):
    """Custom generator for colors
    set_color can be 
     * None, for a basic cycle of blue, red, green, etc.
     * A single hexcolor code ('#000000' for black, etc)
     * 'alpha', 'ai', 'ug1', 'Dsm1' or 'vr', which have default built in colors
     * Otherwise, it assumes set_color is a list of colors to yield

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