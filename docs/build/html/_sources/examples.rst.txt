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

This cell just imports the important packages, and sets up some iPython settings that are useful. ``%matplolib widget`` makes it possible to interact with ``matplotlib`` figures directly, such as panning the camera, zooming in, etc.
