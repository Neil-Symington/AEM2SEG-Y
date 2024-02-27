"""
This script converts airborne electromagnetic (AEM) conductivity section data from a column based ascii format
 ( such as ASEG-GDF) to SEG-Y format on a line-by-line basis.
"""

import numpy as np
import os
import glob
import segyio

# Column indices (Note: These are 1-based indices, not zero-based as is typically used in python)
EASTING_COL = 7  # Column for Easting coordinates (1-based index)
NORTHING_COL = 8  # Column for Northing coordinates (1-based index)
ELEVATION_COL = 11  # Column for Elevation data (1-based index)
LINE_COL = 5  # Column for line numbers (1-based index)
CONDUCTIVITY_COLS = [26, 55]  # Columns for conductivity data (1-based index range)
THICKNESS_COLS = None # Adjust if using thickness of layers directly
#THICKNESS_COLS = [56, 85]  # Columns for layer thickness data (1-based index range)
LAYER_TOP_ELEVATION_COLS = None  # Adjust if using top elevation of layers directly
#LAYER_TOP_DEPTH_COLS = None  # Adjust if using top depth of layers directly
LAYER_TOP_DEPTH_COLS = [86, 115]
Y_RES = 2.0  # Vertical resolution in the output SEG-Y file
MAX_DEPTH = 400.0  # Maximum depth to be represented in the output SEG-Y file
SCALING = 1.0  # scale data byt his factor depending on desired units (i.e. mS/m to S/m or ohm.m to S/m)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INDIR = os.path.join(SCRIPT_DIR, '..', '..', 'data', 'lines') # or define absolute path here
FILE_EXTENSION = "asc" # file extension for data files
OUTDIR = r"..\..\data\segy"
#INFILE = r"..\..\data\inversion.output"
INFILE = None # if INFILE = None then the script assume multiple files in INDIR with a common extension


def thickness_to_depth(thickness):
    """
    Calculates depth top from a thickness array.

    :param thickness: An array of thicknesses.
    :return: A flat array of depth.
    """
    layer_top_depth = np.zeros(shape=thickness.shape, dtype=float)
    layer_top_depth[:, 1:] = np.cumsum(thickness[:, :-1], axis=1)
    return layer_top_depth


def line2segy(line_number, lines, easting, northing, elevation, conductivity, layer_top_elevations, outdir, yres, max_depth):
    """
    Converts line data to SEG-Y format.

    :param line_number: The number of the line being processed.
    :param line: an array of line numbers
    :param easting: Easting coordinates.
    :param northing: Northing coordinates.
    :param elevation: Elevation data.
    :param conductivity: Conductivity data.
    :param layer_top_elevations: Top elevations of layers.
    :param outdir: Output directory.
    :param yres: Vertical resolution.
    :param max_depth: Maximum depth.
    """

    line_inds = np.where(lines == line_number)[0]

    # Define top datum
    ymax = np.ceil(np.max(elevation) / 10.0) * 10.0
    ymin = np.floor(np.min(elevation) / 10.0) * 10.0 - max_depth
    grid_elevations = np.arange(ymin, ymax + yres, yres)[::-1]

    interpolated = -1 * np.ones(shape=(len(grid_elevations), len(line_inds)), dtype=float)

    for i in range(len(line_inds)):
        layer_top_elev = layer_top_elevations[i, :]
        sigma_trace = conductivity[i]
        for j in range(len(layer_top_elev) - 1):
            mask = (grid_elevations < layer_top_elev[j]) & (grid_elevations >= layer_top_elev[j + 1])
            interpolated[mask, i] = sigma_trace[j]

    # Writing to SEG-Y
    path = os.path.join(outdir, f"{line_number}.segy")
    segyio.tools.from_array2D(path, interpolated.T, dt=int(yres * 1e3), delrt=-1 * int(ymax))

    with segyio.open(path, mode='r+', ignore_geometry=True) as f:
        for i, x in enumerate(f.header[:]):
            x.update({segyio.TraceField.GroupX: int(easting[i]),
                      segyio.TraceField.SourceX: int(easting[i]),
                      segyio.TraceField.GroupY: int(northing[i]),
                      segyio.TraceField.SourceY: int(northing[i]),
                      segyio.TraceField.ElevationScalar: int(elevation[i]),
                      segyio.TraceField.ReceiverGroupElevation: int(elevation[i])})

def ascii2segy(file):

    """
    Parse the text file and convert to segy.

    :param file: path to asci file.
    """
    data = np.loadtxt(file, skiprows=1)
    X = data[:, EASTING_COL - 1]
    Y = data[:, NORTHING_COL - 1]
    Z = data[:, ELEVATION_COL - 1]
    lines = data[:, LINE_COL - 1].astype(int)
    conductivity = data[:,CONDUCTIVITY_COLS[0] - 1: CONDUCTIVITY_COLS[1]] * SCALING

    if THICKNESS_COLS is not None:
        thickness = data[:,THICKNESS_COLS[0] - 1:THICKNESS_COLS[1]]
        layer_top_depth = thickness_to_depth(thickness)
        layer_top_elevations = Z[:, np.newaxis] - layer_top_depth

    elif LAYER_TOP_ELEVATION_COLS is not None:
        layer_top_elevations = data[:,LAYER_TOP_ELEVATION_COLS[0] - 1:LAYER_TOP_ELEVATION_COLS[1]]
        print(layer_top_elevations)

    elif LAYER_TOP_DEPTH_COLS is not None:
        layer_top_depth = data[:,LAYER_TOP_DEPTH_COLS[0] - 1: LAYER_TOP_DEPTH_COLS[1]]
        layer_top_elevations = Z[:, np.newaxis] - layer_top_depth

    for line_number in np.unique(lines):
        line2segy(line_number, lines, X, Y, Z, conductivity, layer_top_elevations, OUTDIR, Y_RES, MAX_DEPTH)

# Main execution
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)

if not INFILE is None:
    ascii2segy(INFILE)
else:
    for file in glob.glob(os.path.join(INDIR, "*.{}".format(FILE_EXTENSION))):
        print(file)
        ascii2segy(file)


