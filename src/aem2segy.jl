using SegyIO
using Statistics
using PyPlot
using DelimitedFiles

##
# choose the vertical resampling interval
yres = 4.0;
max_depth = 400.;

##
infile = raw"C:\Users\u77932\Documents\MORPH\data\AEM\AusAEM_2020\earaheedy_desertStrip\run.01\ref.01\inversion.output.dat";
outdir = raw"C:\Users\u77932\Documents\MORPH\data\AEM\AusAEM_2020\segy";

inversion_output = readdlm(infile);

# define column numbers
easting_col = 7;
northing_col = 8;
elevation_col = 9;
conductivity_cols = [23,52] #ranger for 30 layer model;
thickness_cols = [53,82];
line_col = 5;

lines = round.(Int, inversion_output[:,line_col]);

##

function get_line_indices(line)
    inds = findall(x -> x in line, lines)
    return inds
end

function line2segy(line_number)
    # find our actual line indices
    line_inds = get_line_indices(line_number)

    # our first step is to define our datum as 20 m above the highest point in survey

    elevations = inversion_output[line_inds,elevation_col]
    ymax = ceil(maximum(elevations)/10.) *10.

    # create an elevation grid

    ymin = floor(minimum(elevations)/10.)*10. - max_depth

    grid_elevations = collect(range(ymax, ymin, length = Int(ceil((ymax - ymin)/yres))))


    # get conductivity
    sigma = transpose(inversion_output[line_inds,conductivity_cols[1]:conductivity_cols[2]])
    thickness = transpose(inversion_output[line_inds,thickness_cols[1]:thickness_cols[2]])

    ## get coordinates

    easting = inversion_output[line_inds,easting_col]
    northing = inversion_output[line_inds,northing_col]
    ## create our array of interpolated AEM conductivity, negative indicates free air
    interpolated = -1*ones(Float32, length(grid_elevations), length(line_inds))

    # iterate
    for i in 1:length(line_inds)

        layer_top_elev = elevations[i] * ones(size(thickness)[1])
        layer_top_elev[2:end] = layer_top_elev[2:end] - cumsum(thickness[1:end-1,i])
        sigma_trace = sigma[:,i]
        # iterate through our layers
        for j in 1:(length(layer_top_elev) - 1)
            interpolated[(grid_elevations .< layer_top_elev[j]) .& (grid_elevations .>= layer_top_elev[j + 1]),i] .= sigma_trace[j]
        end
    end

    ## try writing these out in depth space

    block = SeisBlock(interpolated)

    ##

    set_header!(block, "SourceX", round.(Int, easting))
    set_header!(block, "SourceY", round.(Int, northing))
    set_header!(block, "GroupX", round.(Int, easting))
    set_header!(block, "GroupY", round.(Int, northing))
    set_header!(block, "ElevationScalar", round.(Int, elevations))
    set_header!(block, "TraceNumWithinLine", Array(1:length(line_inds)))
    set_header!(block, "TraceNumWithinFile", Array(1:length(line_inds)))
    #set_header!(block, "TraceIDCode", round.(Int, fiducial * 10))
    set_header!(block, "RecGroupElevation", round.(Int, elevations))
    set_header!(block, "DelayRecordingTime", Int(ymax) * -1)
    set_header!(block, "dt", Int(yres * 1e3))

    ##
    segy_write(joinpath(outdir,string(line_number) * ".segy"), block)
end
##
for l in unique(lines)
    line2segy(l)
end
