#!/usr/bin/env python

#===============================================================================
#    Copyright 2017 Geoscience Australia
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#===============================================================================

'''
Created on 8/5/2019
@author: Neil Symington

Workflow for converting the aseg gdf data to seg-y
'''

import sys, os, glob
import math
import numpy as np
import aem2segy
# Needs to be obspy 0.10
import obspy

# Get the control file

control_file = sys.argv[1]


# Parse the control file to get the variables
var_dict = aem2segy.parse_control_file(control_file)


try:
    header_file = var_dict['header_file']
except KeyError:
    header_file = None

# Set the variables

# Set directory with AEM line data
AEM_dir = var_dict['AEM_dir']


# Set directory to save the outfiles
outdir = var_dict['segy_dir']

# Create the out directory if it doesn't exit
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Select resampling interval in metres
interval = np.float(var_dict['vertical_interval'])

# Define the max depth

max_depth = np.float(var_dict['max_depth'])

# Define the height datum on which all the data will be hung. This should be
# at least as high as the highest topography in the survey area

datum = np.float(var_dict['datum'])

# Define the new elevation profile onto which the conductivity data will be interpolated
z_new = np.arange(-1 * datum + (interval/2.), max_depth + (interval/2.), interval)

# Set the job number

job = int(var_dict['job_id'])

# If we are running a single line then search for it in the var_dict

try:
    line_no = var_dict['AEM_line'] + '.asc'
except KeyError:
    # Batch mode
    line_no = '*.asc'



for file in glob.glob(os.path.join(AEM_dir, line_no)):
    
    print file
    
    # Get the line name
    
    line = file.split('\\')[-1].split('.')[0]

    
    # Only run the process if the files don't already exist
    if not os.path.isfile(outdir + line + '.segy'):
        
        # Retrieve the data from the asci file

        data_dict = aem2segy.parse_AEM(file, var_dict)
        
        # Define numpy array for continuous, interpolated conductivity data
        
        stream_cond = np.ones((np.shape(data_dict['data'])[0],
                               math.ceil((datum + max_depth)/interval)))

        ### TODO remove for loop and vectorise this subroutine
        for i in range(len(stream_cond)):

            interp_dat = aem2segy.interpolate_layer_data(data_dict['depth_top'][i],  z_new,
                                                         data_dict['data'][i], data_dict['elevation'][i],
                                                         max_depth, datum)


            # If flagged then mask all data below the doi
            if aem2segy.to_bool(var_dict['doi_mask']) == True:
                try:
                    interp_dat = aem2segy.remove_below_doi(interp_dat, z_new,
                                                           data_dict['depth_of_investigation'][i],
                                                           data_dict['elevation'][i])
                except KeyError:
                    print "Invalild input for depth_of_investigation"

            # Assert that the datum is greater than the max elevation

            assert np.max(data_dict['elevation']) < datum

            stream_cond[i] = interp_dat

        # Now write each fiducial into the trace. For this we use the obspy package
        traces = []


        for row in stream_cond:
            trace = obspy.core.trace.Trace(row)
            # This is a hack to try and produce a segy with an acceptably low
            # sampling frequency. I highly recommend to include sample depth information
            # in the header
            trace.stats.delta = interval * 0.001
            traces.append(trace)
        
        # Write into a stream
        stream = obspy.core.stream.Stream(traces)
        
        
        for tr in stream:
            tr.data = np.require(tr.data, dtype=np.float32)

        # Write the stream into a temporary segy
        
        tempoutfile = os.path.join(outdir, line + '_temp.segy')
        
        stream.write(tempoutfile, format = 'SEGY', data_encoding = 1)
        
        # Reopen the segy
        
        st = obspy.segy.core.readSEGY(tempoutfile)
        
        os.remove(tempoutfile)
        
        # Now write in some of the important header information
        
        st.stats.binary_file_header['line_number'] = int(line)
        st.stats.binary_file_header['job_identification_number'] = job
        
        
        # Now write trace header information


        for i in range(len(st)):
            st[i].stats.segy.trace_header.trace_sequence_number_within_line = i
            st[i].stats.segy.trace_header.trace_number_within_the_original_field_record = data_dict['fiducial'][i]
            st[i].stats.segy.trace_header.trace_number_within_the_ensemble = data_dict['fiducial'][i]
            st[i].stats.segy.trace_header.ensemble_number = data_dict['fiducial'][i]
            st[i].stats.segy.trace_header.source_coordinate_x = data_dict['easting'][i]
            st[i].stats.segy.trace_header.source_coordinate_y = data_dict['northing'][i]
            st[i].stats.segy.trace_header.group_coordinate_x = data_dict['easting'][i]
            st[i].stats.segy.trace_header.group_coordinate_y = data_dict['northing'][i]
            st[i].stats.segy.trace_header.receiver_group_elevation = data_dict['elevation'][i]
            st[i].stats.segy.trace_header.datum_elevation_at_receiver_group = 0
            st[i].stats.segy.trace_header.datum_elevation_at_source = 0
            st[i].stats.segy.trace_header.delay_recording_time = -1 * datum
        
        
        outfile = os.path.join(outdir, line + '.segy')
        
        st.write(outfile, format = 'SEGY', data_encoding = 1)
