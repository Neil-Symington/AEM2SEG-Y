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

Functions for converting the aseg gdf data to seg-y
'''


import numpy as np
import ast
from scipy import interpolate


# Define a function for parsing the control file

# From  https://stackoverflow.com/questions/715417/converting-from-a-string-to-boolean-in-python
def to_bool(value):
    """
    Converts 'something' to boolean. Raises exception for invalid formats
    Possible True  values: 1, True, "1", "TRue", "yes", "y", "t"
    Possible False values: 0, False, None, [], {}, "", "0", "faLse", "no", "n", "f", 0.0, ...
    :param value: string

    :return: boolean: True or False
    """
    if str(value).lower() in ("yes", "y", "true",  "t", "1"): return True
    if str(value).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "[]", "{}"): return False
    raise Exception('Invalid value for boolean conversion: ' + str(value))

def RepresentsInt(s):
    """
    check if a string can be represented by an interger
    :param s: string
    :return: Boolean
    """
    try:
        int(s)
        return True
    except ValueError:
        return False

def parse_control_file(infile):
    """
    A function for parsing the control file
    :param infile: path to control file
    :return:
    dictionary with key infomation needed to convert
    """
    # Create diciotnary
    var_dict = {}

    # OPen the file
    with open(infile, 'r') as f:

        # Iterate through the lines in the file
        for line in f:
            s = line.strip()
            
            # Pass for empty lines
            if len(s.strip()) == 0:
                pass
            
            # Pass for comments
            elif s.strip()[0] == '#':
                pass
            # Otherwise split the string on the equals and add to the dictionary with the key word as the key
            else:
                l = s.split('=')
                var_dict[l[0].strip()] = l[1].strip()
    return var_dict


def listify_data_columns(string):
    """
    Take a string representing a range of integer values (e.g. 43-72) and create a pythonic range
    :param string:
    :return:
    """
    d1 = int(string.split('-')[0])
    d2 = int(string.split('-')[1])
    
    return range(d1, d2+1)

def check_range_string(s):
    """
    Check if a string is a valid range of type "34-67"
    :param s: string
    :return: boolean
    """

    L = s.split("-")

    if (len(L) != 2):
        return False
    elif (RepresentsInt(L[0])) & (RepresentsInt(L[1])):
        return True
    else:
        return False

def parse_AEM(AEM_file, var_dict):
    """
    This function parses the AEM asci file
    :param AEM_file: path to file
    :param var_dict: dictionary with column information
    :return: dictionary with numpy arrays for numerical data with keyword from var_dict as the key
    """

    # Dictionary for column numbers
    col_dict = {}
    data_dict = {}

    # Find the column indices and rewrite them into a dictionary as lists
    # WE will use this dictionary to extract the data
    
    columns = ['easting', 'northing', 'elevation', 'fiducial', 'depth_of_investigation', 'data', 'depth_top']

    # Flags
    depth_top_in_file = True

    # Iterate through columns
    for item in columns[:-2]:

        try:
            # Extract entry from dictionary
            entry = var_dict[item]

            # check the data type is a string, integer or another
            if type(entry) == int:

                col_dict[item] = entry

            elif type(entry) == str:
                # get from header file
                if RepresentsInt(entry):
                    col_dict[item] = entry
                else:
                    print "Invalid string entry for ", item

        except KeyError:
            if item == 'depth_of_investigation':
                data_dict[item] = None
            else:
                print "Please create a valid control file entry for variable ", item
                return None

    for item in columns[-2:]:



        try:
            # Extract entry from dictionary
            entry = var_dict[item]

            # check the data type is a string with a range or mapped to the .hdr file
            if type(entry) == str:

                # Checck if it is a valid range

                if check_range_string(entry):

                    col_dict[item] = listify_data_columns(entry)
                else:

                    # Raise flag
                    depth_top_in_file = False
                    col_dict[item] = np.array(ast.literal_eval(entry))

        except KeyError:
            print "Please create a valid control file entry for variable ", item

    # Convert to pythonic indexing
    
    first_col = 1

    # Search for first_col keyword in case it has been included
    
    if 'first_col' in var_dict.keys():
        first_col = int(var_dict['first_col'])
    
    t = (2 - first_col)
    
    # Now subtract the value to all list elements in the cols directory
    
    for item in columns[:-2]:

        # Get the columns

        cols = int(col_dict[item]) - t

        # Extract as a numpy array

        data_dict[item] = np.loadtxt(AEM_file, usecols= cols)

    # Get the data cols
    cols = [int(x) - t for x in col_dict['data']]

    data_dict['data'] = np.loadtxt(AEM_file, usecols=cols)
    
    # If data is resistivity, convert to conductivity
    if to_bool(var_dict['resistivity']):

        data_dict['data'] = 1./data_dict['data']
    
    # Multiply data by the scaling factor 

        data_dict['data'] = data_dict['data'] * np.float(var_dict['scaling_factor'])

    # If the depth tops are in the file extract

    if depth_top_in_file:

        cols = [x - t for x in col_dict['depth_top']]


        # Extract and tile
        data_dict['depth_top'] = np.loadtxt(AEM_file, usecols = cols)

    else:
        # Otherwise extract the parsed list and tile it to fit the data

        data_dict['depth_top'] = np.tile(np.array(col_dict['depth_top']),
                                         (data_dict['data'].shape[0],1))

    # Assert that the depth top and data array are the same shape

    assert data_dict['depth_top'].shape == data_dict['data'].shape

    return data_dict

# Function for nulling all values below the doi
def remove_below_doi(interpolated_data, z_new, doi, elevation):

    """

    :param interpolated_data: numpy array with interpolated data
    :param z_new: new elevation intervals for segy trace
    :param doi: float with fiducial depth of investigation
    :param elevation: float fiducial with elevation
    :return:
    interpolated_data with below doi values changed to -1.
    """

    doi_elevation = -1 * (elevation - doi)
    
    # Find the indices that are below the depth of investigation
    interpolated_data[np.where(z_new > doi_elevation)] = -1

    return interpolated_data

# Interpolate so that we have a continuously spaced data

def interpolate_layer_data(depth_top, z_new, dat, elev, max_depth, datum):
    # First find layer bottom (by adding a small delta d)

    depth_bottom = depth_top[1:] - 0.01

    # Now add the layer tops and bottoms into a single array and produce a
    # corresponding conductivity array

    # The aim is to book end each layer
    z = []
    new_dat = []
    for i in range(len(depth_bottom)):
        z.append(depth_top[i])
        z.append(depth_bottom[i])
        new_dat.append(dat[i])
        new_dat.append(dat[i])

    # Convert the depth to elevation (where negative values are above msl)
    z = [x - elev for x in z]

    # Finally bookend the air and give it a conductivity of 0

    z.insert(0, z[0] - 0.01)
    z.insert(0, datum * -1)

    new_dat.insert(0, -1)
    new_dat.insert(0, -1)

    # Now bookend the bottom half-space to the max depth
    z.append(z[-1] + 0.01)
    z.append(-1 * max_depth * -1)

    new_dat.append(dat[-1])
    new_dat.append(dat[-1])

    f = interpolate.interp1d(z, new_dat)

    interpolated_dat = f(z_new)

    return interpolated_dat

