import numpy as np
import segyio
import os

# choose the vertical resampling interval
yres = 4.0
max_depth = 400.
scaling = 0.001 # convert from mS/m to S/m
##
infile = r"C:\Users\u77932\Documents\GAB\data\AEM\2017_survey\lci\Injune_WB_MGA55.dat"
outdir = r"C:\Users\u77932\Documents\GAB\data\AEM\2017_survey\lci\segy"

# define column numbers with python indexing (i.e. columns 1 has index of 0)
easting_col = 5
northing_col = 6
elevation_col =7
line_col = 4
conductivity_cols = [42,71] #range for 30 layer model;

# depending on the inversion, the data may have thicknesses, layer_top_depths or layer_top_elevations
# pick one and make the others = None
thickness_cols = None
layer_top_elevation_cols = [12,41]
layer_top_depth_cols = None

X = np.loadtxt(infile, usecols=[easting_col])
Y = np.loadtxt(infile, usecols=[northing_col])
Z = np.loadtxt(infile, usecols=[elevation_col])
conductivity = np.loadtxt(infile, usecols=np.arange(conductivity_cols[0], conductivity_cols[1]+1)) * scaling

if thickness_cols is not None:
    thicknesses = np.loadtxt(infile, usecols=np.arange(thickness_cols[0],thickness_cols[1]+1))
    layer_top_elevations = np.zeros(shape = conductivity.shape, dtype = float)
    layer_top_elevations[:,0] = Z
    layer_top_elevations[:,1:] = layer_top_elevations[:,0] - np.cumsum(thicknesses[:,-1])

elif layer_top_elevation_cols is not None:
    layer_top_elevations = np.loadtxt(infile, usecols=np.arange(layer_top_elevation_cols[0],
                                                                layer_top_elevation_cols[1]+1))
elif layer_top_depth_cols is not None:
    layer_top_depths = np.loadtxt(infile, usecols=np.arange(layer_top_depth_cols[0],
                                                            layer_top_depth_cols[1] + 1))
    layer_top_elevations = Z - layer_top_depths
else:
    print("please define column numbers for one of thickness, layer top elevation or layer top depths")


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
    layer_top_elevation = layer_top_elevations[line_inds]
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

        layer_top_elev = layer_top_elevation[i,:]

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

