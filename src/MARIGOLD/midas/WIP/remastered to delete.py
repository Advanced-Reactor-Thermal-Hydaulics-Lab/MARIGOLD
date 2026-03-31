from ..config import *



from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import numpy as np
from openpyxl import load_workbook
from openpyxl.utils import coordinate_to_tuple, get_column_letter
from pathlib import Path
import re
import shutil
import subprocess
import xlwings as xw


# Root directory
BASE_DIR = Path(r"D:\NEUP")

# Angles: folders whose names end with "deg"
ANGLE_SUFFIX = "deg"

# Subdirectory that contains conductivity data under each ANGLE folder
CONDUCTIVITY_SUBDIR = "Conductivity"

# Name of the input file in each AZIMUTH directory
INP_FILENAME = "Input.inp"

# Sample frequency (Hz)
SAMPLE_FREQ = 50000

# Order of radial positions (TAB files) to compile into each AZIMUTH.TAB
RADIUS_ORDER = [
    "r9o",
    "r85o",
    "r8o",
    "r7o",
    "r6o",
    "r5o",
    "r4o",
    "r3o",
    "r2o",
    "r1o",
    "r0o",  # or r0i or r0
    "r1i",
    "r2i",
    "r3i",
    "r4i",
    "r5i",
    "r6i",
    "r7i",
    "r8i",
    "r85i",
    "r9i",
]

# Azimuth folder names
AZIMUTHS = ["90", "67.5", "45", "22.5", "00"]

# Mapping from AZIMUTH string to top-left Excel cell for pasting
AZIMUTH_CELL_MAP = {
    "90": "B9",
    "67.5": "B56",
    "45": "B105",
    "22.5": "B152",
    "00": "B201",
}

