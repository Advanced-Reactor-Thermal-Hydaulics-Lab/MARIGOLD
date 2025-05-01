Example Analysis
================

----------

MARIGOLD includes some example files (under MARIGOLD/examples) to help you get started with two-phase flow analysis. Here we will discuss and explain the example files in more detail.

example_analysis.ipynb
######################

First thing to note is that this is a interactive python, or Jupytr, notebook. This allows each cell to be run individually. Typically we use VSCode, with the Jupytr notebook extension installed, to run these.

Lets take a look at the first cell::

   # Imports and iPython magic setup
   import numpy as np
   import matplotlib as mpl
   import matplotlib.pyplot as plt
   import os
   import MARIGOLD as mgd
   from MARIGOLD import Condition
   %matplotlib widget
   %load_ext autoreload
   %autoreload 2

This cell just imports the important packages, and sets up some iPython settings that are useful. Usually ``MARIGOLD`` is imported as ``mgd`` for easy reference later, while the entire ``Condition`` object is imported for better type hinting. ``%matplolib widget`` makes it possible to interact with ``matplotlib`` figures directly, such as panning the camera, zooming in, etc. ``load_ext autoreload`` and ``autoreload 2`` are useful if you're changing MARIGOLD functions and running them in iPython, as you don't need to re-import them after making changes

The next cell::

   # Load Data from spreadsheets
   
   cwd = os.path.abspath('')
   dat_name = 'example_database.dat'
   path_to_dump = os.path.join(cwd, "example_sheets") # Assumes the "example_sheets" directory is present in the cwd

   # If the .dat file already exists, we could just load it immediatiately by setting this to False
   refetch = True
   if not os.path.isfile(os.path.join(cwd, dat_name)) or refetch:
       mgd.extractLocalDataFromDir(path_to_dump, dump_file=dat_name, sheet_type='adix_template4')
   
   database : list[Condition] = mgd.loadData(dat_name)

This cell is responsible for loading in the experimental data. ``extractLocalDataFromDir`` is the star of the show here. This function is responsible for looking through a directory and finding any .xlsx files that have probe data. It then saves them to a .dat file, in this example, ``example_database.dat``. This allows easy transfer of data, without having to move the .xlsx files around. The .dat file is loaded by the ``loadData`` function.

The database is now filled with ``Condition`` objects. These are the primary way we will be interacting with the data. Data at a specific  :math:`(r, \varphi)` can be accessed via calling the condition, ``cond(phi, r)`` where :math:`varphi` is in radians and :math:`r` is nondimensional. There are also lots of built in methods for interacting and manipulating the data.

For instance, a common task in data analysis is looking at area-averaged parameters. The next cell is an example of making a table of these. ::
      
   database.sort(key=lambda cond: cond.jgloc+cond.jf) # Sort the database by jf and jg, in ascending order
   
   # This will print out the area averages
   print("jf\tjgloc\t⟨α⟩\t⟨α vg⟩\t⟨(1-α) vf⟩\tε_jg")
   for cond in database:
       cond.mirror(method='axisym')
       print(f"{cond.jf:0.2f}\t{cond.jgloc:0.3f}\t{cond.area_avg('alpha'):0.3f}\t{cond.area_avg('alpha_ug1'):0.3f}\t{cond.area_avg('jf'):0.3f}\t\t{ (cond.area_avg('alpha_ug1') - cond.jgloc)/(cond.jgloc) *100 :0.1f}")
   

Before area-averaging, it is critical that each ``Condition`` object is mirrored. This ensures there is data around the entire cross section. For 90 degree data, the ``method`` we use is ``axisym``, for horizontal it would be ``sym90``. The area-averaging is performed by the ``.area_avg()`` method of the ``Condition`` object. Any ``param`` can be area-averaged, such as ``alpha``, ``alpha_ug1``, etc. A full list of params can be found by calling ``mgd.print_tab_keys()``. 

The next cell is an example agreement plot for the gas benchmark. ::
   
   # Gas benchmark
   fig, ax = plt.subplots(figsize = (5, 5))
   plt.rcParams.update({'font.size': 12})
   plt.rcParams["font.family"] = "Times New Roman"
   plt.rcParams["mathtext.fontset"] = "cm"
   for cond in database:
   
       # mpbl = ax.scatter(cond.jgloc, cond.area_avg('alpha_ug1'), marker=cond.marker_type, color = cond.marker_color)
       mpbl = ax.scatter(cond.jgloc, cond.area_avg('alpha_ug1'), marker='o', color = 'r', edgecolors='black', s=30)

   ax.set_xscale('log')
   ax.set_yscale('log')
   plt.plot([-10, 10], [-11, 11], '--', color='grey')
   plt.plot([-10, 10], [-9, 9], '--', color='grey', label = "± 10%")
   plt.plot([-10, 10], [-10, 10], linestyle='dotted', color='black')
   
   ax.scatter([-10], [-10], linestyle='None', marker='o', color = 'r', edgecolor='black', s=30, label = "Data")
   
   plt.xlim(0.075, 0.35)
   plt.ylim(0.075, 0.35)
   
   ax.set_xticks([0.075, 0.1, 0.2, 0.3])
   ax.set_yticks([0.075, 0.1, 0.2, 0.3])
   ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
   ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
   
   ax.set_aspect('equal')
   
   plt.xlabel(r'$\langle j_{g} \rangle _{rotameters}\ m/s$')
   plt.ylabel(r'$\langle \alpha v_{g} \rangle_{probe}\  m/s$')
   
   # plt.colorbar(mpbl, label=r'$\langle j_{g} \rangle \ m/s$', ticks=np.arange(0, 0.51, 0.1), boundaries = np.arange(0, 0.51, 0.01), values = np.arange(0, 0.5, 0.01))
   
   plt.legend()
   plt.tight_layout()
   plt.savefig(r".\gas_benchmark.png", dpi=500)
   plt.show()
   
This is honestly more ``matplotlib`` stuff than MARIGOLD, so no additional comments will be made.

Plotting is an important feature in MARIGOLD. Below is an example of how to do line plots. Again, this is a method on a ``Condition`` object, so it has some similarities with ``.avea_avg()``. ::
   
   for cond in database:
       for param in ['alpha', 'ai', 'ug1', 'Dsm1', 'vf', 'vr']:
           # Usually the minimum of the graph is set by the minimum of the param value. But for Dsm, the diameter doesn't 
           # exactly go to zero at the wall, so when we plot it makes more sense to not go to 0
           if param == 'Dsm1':
               set_min = np.floor(cond.min('Dsm1', nonzero = True)) - 1
           else:
               set_min = cond.min(param)
           
           # This function will plot the data down the 90° line, with the r data plotted on the x axis. 
           # The colors are set by the function based on the param type. Black for void, red for vg, etc.
           # 6.35 x 3 is the best size to fit 4 graphs all on the same ppt slide
           cond.plot_profiles2(param, x_axis = 'r', const_to_plot = [90], title = False, fig_size = (6.35, 3), show = True, cs = 'infer', set_min = set_min)

Hopefully this gave you a start on performing two-phase flow analysis with MARIGOLD. Some things that weren't covered in this tutorial include contour plotting with ``plot_contour`` and interpolation. 