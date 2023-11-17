#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:41:49 2023

@author: amt
"""

from scipy.io import loadmat

# Load the MATLAB .mat file
test = loadmat('test.mat')

loadTide={}
for ii in range(16):
    loadTide[ii] = {
        'Station': {'Name': test['loadTides'][0][ii][0][0][0][0][0],
                    'Lat': test['loadTides'][0][ii][0][0][0][1][0][0],
                    'Lon': test['loadTides'][0][ii][0][0][0][2][0][0],
                    'Height': test['loadTides'][0][ii][0][0][0][3][0][0]},
        'Type': test['loadTides']['Type'][0][ii][0],
        'Phase': test['loadTides']['Phase'][0][ii][0],
        'Info': test['loadTides']['Info'][0][ii][0][0][0][0],
        'Darwin': test['loadTides']['Darwin'][0][ii][0],
        'DoodsonNum': test['loadTides']['DoodsonNum'][0][ii][0],
        'gravLoadTide': test['loadTides']['gravLoadTide'][0][ii][0],
        'potentialHeight': test['loadTides']['potentialHeight'][0][ii][0],
        'displacement': test['loadTides']['displacement'][0][ii][0],
        'tilt': test['loadTides']['tilt'][0][ii][0],
        'strain': test['loadTides']['strain'][0][ii][0],
    }
    
def combine_tide_struct(*args,flag='s'):
    
    """
    Combines multiple load tide structures based on the Darwin Symbol and type.

    Call: output = combineTideStruct(loadStruct1, loadStruct2, ...)
    
    Concatenates all tide structures into a single tide structure array.
    After concatenation, it looks by type/Darwin pairs to combine (sum) tide
    structures.
    
    Call: output = combineTideStruct(flag, loadStruct1, ...)
    
    Sets the combine flag. If flag is 's' (default), then all types are kept
    separate. If the flag is 'c' (combine), then it will look for load and
    earth tides, and sum them together changing the type to combined. The
    separate earth and load are still kept in their original forms, but all
    previous combined tides are thrown out.
    
    Output: A single (Nx1) load structure, one element for each unique Darwin
    Symbol/type pair.
    
    Written 1/2012 by Charlie Sievers @ UNAVCO
    """
    
    if len(args) < 1:
        raise ValueError('combine_tide_struct expects an input.')
    
    if not isinstance(flag, str) or len(flag) != 1 or sum(c in 'sc' for c in flag) == 0:
        print(flag)
        raise ValueError('The flag must be either s (separate) or c (combine).')
    
    if len(args) < 1:
        raise ValueError('No tidal structure passed to combineTideStruct')
    
    # Check to make sure there aren't any body tides
    good_types = ['EarthTide', 'LoadTide', 'CombinedTide']
    for arg in args:
        for jj in range(len(arg)):
            if arg[jj]['Type'] not in good_types:
                raise ValueError('CombineTideStruct is only valid for Earth, Load, and Combined Tides')
    
    # Error checking inputs
    for ii in range(len(args)):
        if not isinstance(args[ii], dict):
            print(args[ii])
            raise ValueError('All inputs (other than the flag) must be dictionaries.')
        if ii == 0:
            fnames = list(args[ii][0].keys())
        else:
            fnames2 = list(args[ii][0].keys())
            if set(fnames) != set(fnames2):
                raise ValueError('All input structures must be the same type.')
                
    # Combine all inputs into one array structure
    input_struct = args[0]
    print(len(input_struct))
    for struct in args[1:]:
        input_struct.update(struct)
    print(len(input_struct))
    
    # Rest of your translation, continuing the logic with if-else blocks, loops, and function calls
    # ...
    # Your code here
    # ...

    return

out=combine_tide_struct(loadTide,loadTide,loadTide,flag='c')
