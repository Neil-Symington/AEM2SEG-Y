import numpy as np
import segyio
import os

# choose the vertical resampling interval
yres = 4.0;
max_depth = 500.;

##
infile = r"C:\Users\u77932\Documents\NGS\data\AEM\EarheedyDesertStrip\GA_layer_earth_inversion\AusAEM_20_Earaheedy&DesertStrip_GA_layer_earth_inversion\Earaheedy& DesertStrip_AusAEM2020_GA_vsum_inversion.dat"
outdir = r"C:\Users\u77932\Documents\NGS\data\AEM\EarheedyDesertStrip\GA_layer_earth_inversion\AusAEM_20_Earaheedy&DesertStrip_GA_layer_earth_inversion"

# define column numbers
easting_col = 6
northing_col = 7
elevation_col = 8
conductivity_cols = [22,51] #ranger for 30 layer model;
thickness_cols = [52,81]
line_col = 4

X = np.loadtxt(infile, usecols=[easting_col])
Y = np.loadtxt(infile, usecols=[northing_col])
Z = np.loadtxt(infile, usecols=[elevation_col])
conductivity = np.loadtxt(infile, usecols=np.arange(conductivity_cols[0], conductivity_cols[1]+1))
thicknesses = np.loadtxt(infile, usecols=np.arange(thickness_cols[0],thickness_cols[1]+1))
lines = np.loadtxt(infile, usecols=[line_col]).astype(int)

##

def get_line_indices(line):
    inds = np.where(lines == line)[0]
    return inds


def line2segy(line_number):
    # find our actual line indices
    line_inds = get_line_indices(line_number)

    easting = X[line_inds]
    northing = Y[line_inds]
    elevation = Z[line_inds]
    thickness = thicknesses[line_inds]
    sigma = conductivity[line_inds]

    # our first step is to define our datum as 20 m above the highest point in survey

    ymax = np.ceil(np.max(elevation)/10.) *10.

    # create an elevation grid

    ymin = np.floor(np.min(elevation)/10.)*10. - max_depth

    grid_elevations = np.arange(ymin, ymax, yres)[::-1]

    # get conductivity

    ## create our array of interpolated AEM conductivity, negative indicates free air
    interpolated = -1*np.ones(shape = (len(grid_elevations), len(line_inds)), dtype = float)

    # iterate
    for i in range(len(line_inds)):

        layer_top_elev = elevation[i] * np.ones(shape = thickness.shape[1])

        layer_top_elev[1:] = layer_top_elev[0] - np.cumsum(thickness[i,:-1])

        #layer_top_elev = layer_top_elevation[i,:]
        sigma_trace = sigma[i]
        # iterate through our layers
        for j in range(len(layer_top_elev) - 1):
            mask = (grid_elevations < layer_top_elev[j]) & (grid_elevations >= layer_top_elev[j + 1])
            interpolated[mask,i] = sigma_trace[j]

    path = os.path.join(outdir, "{}.segy".format(line_number))
    segyio.tools.from_array2D(path, interpolated.T, dt  = int(yres * 1e3), delrt=-1*int(ymax))

    with segyio.open(path, mode = 'r+', ignore_geometry=True) as f:
        for i, x in enumerate(f.header[:]):
            x.update({segyio.TraceField.GroupX: int(easting[i])})
            x.update({segyio.TraceField.SourceX: int(easting[i])})
            x.update({segyio.TraceField.GroupY: int(northing[i])})
            x.update({segyio.TraceField.SourceY: int(northing[i])})
            x.update({segyio.TraceField.ElevationScalar: int(elevation[i])})
            x.update({segyio.TraceField.ReceiverGroupElevation: int(elevation[i])})


for l in np.unique(lines):

    if not str(l).startswith('9'): # avoid high-altitude lines
        line2segy(l)

