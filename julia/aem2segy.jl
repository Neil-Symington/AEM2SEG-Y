using SegyIO
using NetCDF
using Statistics
using PyPlot
## define variables

#lines = [6001002,6001003,6001004]
yres = 4.0

##
ncfile = "/home/nsymington/Documents/GA/AEM/AusAEM/AusAEM_02_GA_layer_earth_inversion/vectorsum/AusAEM_WA_2xxx_regional.nc"

lines = ncread(ncfile, "line")

ncinfo(ncfile)

##

function get_line_indices(line)
    nclines = ncread(ncfile, "line")
    inds = findall(x -> x in line, nclines) .- 1
    linds = findall(x -> x in inds, ncread(ncfile, "line_index"))
    return linds[2:end]
end
function line2segy(line_number)
    # find our actual line indices
    line_inds = get_line_indices(line_number)

    # our first step is to define our datum as 20 m above the highest point in survey

    elevations = ncread(ncfile, "elevation")[line_inds]
    ymax = ceil(maximum(elevations)/10.) *10.

    # create an elevation grid

    ymin = floor(minimum(elevations)/10.)*10. - 600.

    grid_elevations = collect(range(ymax, ymin, length = Int(ceil((ymax - ymin)/yres))))

    ##

    # get conductivity
    σ = ncread(ncfile, "conductivity")[:,line_inds]
    thickness = ncread(ncfile, "thickness")[:,line_inds]

    ## get coordinates

    easting = ncread(ncfile, "easting")[line_inds]
    northing = ncread(ncfile, "northing")[line_inds]
    fiducial = ncread(ncfile, "fiducial")[line_inds]
    ## create our array of interpolated AEM conductivity, negative indicates free air
    interpolated = -1*ones(Float32, length(grid_elevations), length(line_inds))

    # iterate
    for i in 1:length(line_inds)
        layer_top_elev = elevations[i] * ones(size(thickness)[1])
        layer_top_elev[2:end] = layer_top_elev[2:end] - cumsum(thickness[1:end-1,i])
        σ_trace = σ[:,i]
        # iterate through our layers
        for j in 1:(length(layer_top_elev) - 1)
            interpolated[(grid_elevations .< layer_top_elev[j]) .& (grid_elevations .>= layer_top_elev[j + 1]),i] .= σ_trace[j]
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
    segy_write(string(line_number) * ".segy", block)
end
##
for l in lines
    line2segy(l)
end
