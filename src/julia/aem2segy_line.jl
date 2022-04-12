using SegyIO
using DelimitedFiles

function line2segy(infile)
    print(infile);
    inversion_output = readdlm(infile);
    line_number = unique(round.(Int, inversion_output[:,line_col]))[1];

    # our first step is to define our datum as 20 m above the highest point in survey

    elevations = inversion_output[:,elevation_col]
    ymax = ceil(maximum(elevations)/10.) *10.

    # create an elevation grid

    ymin = floor(minimum(elevations)/10.)*10. - max_depth

    grid_elevations = collect(range(ymax, ymin, length = Int(ceil((ymax - ymin)/yres))))

    # get conductivity
    sigma = transpose(inversion_output[:,conductivity_cols[1]:conductivity_cols[2]])
    thickness = transpose(inversion_output[:,thickness_cols[1]:thickness_cols[2]])
    #layer_top_elevation = inversion_output[:,layer_top_elevation_cols[1]:layer_top_elevation_cols[2]]

    ## get coordinates

    easting = inversion_output[:,easting_col]
    northing = inversion_output[:,northing_col]
    ## create our array of interpolated AEM conductivity, negative indicates free air
    interpolated = -1*ones(Float32, length(grid_elevations), length(easting))

    # iterate
    for i in 1:length(easting)

        layer_top_elev = elevations[i] * ones(size(thickness)[1])
        layer_top_elev[2:end] = layer_top_elev[2:end] - cumsum(thickness[1:end-1,i])
        #layer_top_elev = layer_top_elevation[i,:]
        sigma_trace = sigma[:,i]
        # iterate through our layers
        for j in 1:(length(layer_top_elev) - 1)
            interpolated[(grid_elevations .< layer_top_elev[j]) .& (grid_elevations .>= layer_top_elev[j + 1]),i] .= sigma_trace[j]
        end
    end

    ## try writing these out in depth space

    block = SeisBlock(interpolated)

    set_header!(block, "SourceX", round.(Int, easting))
    set_header!(block, "SourceY", round.(Int, northing))
    set_header!(block, "GroupX", round.(Int, easting))
    set_header!(block, "GroupY", round.(Int, northing))
    set_header!(block, "ElevationScalar", round.(Int, elevations))
    set_header!(block, "TraceNumWithinLine", Array(1:length(easting)))
    set_header!(block, "TraceNumWithinFile", Array(1:length(easting)))
    #set_header!(block, "TraceIDCode", round.(Int, fiducial * 10))
    set_header!(block, "RecGroupElevation", round.(Int, elevations))
    set_header!(block, "DelayRecordingTime", Int(ymax) * -1)
    set_header!(block, "dt", Int(yres * 1e3))

    ##
    segy_write(joinpath(outdir,string(line_number) * ".segy"), block)
end

##
# choose the vertical resampling interval
yres = 4.0;
max_depth = 500.;

##
indir = raw"C:\Users\u77932\Documents\LEB\data\AEM\Frome\galei\lines2";
outdir = raw"C:\Users\u77932\Documents\LEB\data\AEM\Frome\galei\segy";

cd(indir);
files = readdir();


# define column numbers
easting_col = 7;
northing_col = 8;
elevation_col = 9;
conductivity_cols = [24,53] #ranger for 30 layer model;
thickness_cols = [54,83];
#layer_top_elevation_cols = [14 43]
line_col = 5;



##
for fname in files

    line2segy(fname)
end
