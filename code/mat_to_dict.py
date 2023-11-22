#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 11:41:49 2023

@author: amt
"""

from scipy.io import loadmat
import numpy as np

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
earthTide={}
for ii in range(16):
    earthTide[ii] = {
        'Station': {'Name': test['loadTides'][0][ii][0][0][0][0][0],
                    'Lat': test['loadTides'][0][ii][0][0][0][1][0][0],
                    'Lon': test['loadTides'][0][ii][0][0][0][2][0][0],
                    'Height': test['loadTides'][0][ii][0][0][0][3][0][0]},
        'Type': 'EarthTide',
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

    
def check_mandatory_fields(dictionary, mandatory_fields):
    """
    Check if a dictionary contains all mandatory fields.

    Args:
    - dictionary (dict): The dictionary to be checked for mandatory fields.
    - mandatory_fields (list): A list of strings representing mandatory fields.

    Returns:
    - bool: True if all mandatory fields are present in the dictionary, False otherwise.
    """
    missing_fields = [field for field in mandatory_fields if field not in dictionary]
    if missing_fields:
        print(f"Missing fields: {', '.join(missing_fields)}")
        return False
    return True

def check_other_fields(mandatory_fields, dictionary):
    """
    Check for missing fields in a dictionary against a list of mandatory fields.

    Args:
    - mandatory_fields (list): A list of strings representing mandatory fields.
    - dictionary (dict): The dictionary to be checked for missing mandatory fields.

    Returns:
    - list: A list of missing fields not found in the dictionary.
    """
    missing_fields = [field for field in mandatory_fields if field not in dictionary]
    return missing_fields

def compare_fields(dict1, dict2, fields_to_compare):
    """
    Compare specific fields between two dictionaries to check if their values are equal.

    Args:
    - dict1 (dict): First dictionary for comparison.
    - dict2 (dict): Second dictionary for comparison.
    - fields_to_compare (list): List of field names to compare.

    Returns:
    - bool: True if the values of specified fields in both dictionaries are equal, False otherwise.
    """

    for field in fields_to_compare:
        if dict1.get(field) != dict2.get(field):
            return False
    return True

def merge_lists(nested_list):
    """
    Combines lists within a nested list that share common elements.

    Args:
    nested_list (list): A nested list of number pairs.

    Returns:
    list: A modified nested list where lists with common elements are merged.
    """

    def merge_helper(list1, list2):
        return list(set(list1) | set(list2))

    merged = True

    while merged:
        merged = False
        i = 0

        while i < len(nested_list):
            j = i + 1
            while j < len(nested_list):
                if set(nested_list[i]) & set(nested_list[j]):
                    nested_list[i] = merge_helper(nested_list[i], nested_list[j])
                    del nested_list[j]
                    merged = True
                else:
                    j += 1
            i += 1

    return nested_list

def tide_phase_rotation(tide, phase):
    """
    Rotates the phase of all the components in the tidal structure.

    Args:
    - tide (list): List containing tidal structure.
    - phase (str or int): If phase is an integer, it rotates the phase by that value in degrees.
                          If phase is 'g' or 'l', it converts the phase to Greenwich or Local respectively.

    Returns:
    - list: Modified tidal structure with updated phases.
    """
    if len(tide) == 0 or phase is None:
        raise ValueError("Bad Elements")

    long = tide[0]['Station']['Long']

    if isinstance(phase, str):
        if phase.lower() == 'l':
            angle = -long
            for item in tide:
                if item['Phase'].lower() == 'greenwich':
                    rot = np.exp(1j * item['DoodsonNum'][0] * angle * np.pi / 180)
                    item = apply_rotation(item, rot)
                    item['Phase'] = 'Local'
        elif phase.lower() == 'g':
            angle = long
            for item in tide:
                if item['Phase'].lower() == 'local':
                    rot = np.exp(1j * item['DoodsonNum'][0] * angle * np.pi / 180)
                    item = apply_rotation(item, rot)
                    item['Phase'] = 'Greenwich'
        else:
            raise ValueError(f"Bad phase: {phase}")
    else:
        if isinstance(phase, int):
            angle = np.array(phase)
        else:
            angle = np.array(phase)
            if angle.size != len(long):
                raise ValueError("Bad Phase")

        for item in tide:
            rot = np.exp(1j * angle * np.pi / 180)
            item = apply_rotation(item, rot)
    return tide


def apply_rotation(tide, rot):
    """
    Applies rotation to specified fields in the tidal structure.

    Args:
    - tide (dict): Tidal structure dictionary.
    - rot (complex): Rotation value.

    Returns:
    - dict: Modified tidal structure with rotated values.
    """
    fields_ignore = ['Station', 'Type', 'Phase', 'Info', 'Darwin', 'DoodsonNum']
    fields = list(tide.keys())
    fields = [field for field in fields if field not in fields_ignore]

    for field in fields:
        value = tide[field]
        if isinstance(value, (int, float)):
            tide[field] = value * rot
    return tide

def compare_fields(dict1, dict2, fields_to_compare):
    """
    Compare specific fields between two dictionaries to check if their values are equal.

    Args:
    - dict1 (dict): First dictionary for comparison.
    - dict2 (dict): Second dictionary for comparison.
    - fields_to_compare (list): List of field names to compare.

    Returns:
    - bool: True if the values of specified fields in both dictionaries are equal, False otherwise.
    """
    c=0
    for field in fields_to_compare:
        if np.sum(dict1.get(field) == dict2.get(field)):
            c+=1
    if c==len(fields_to_compare):
        return 1
    else:
        return 0
    
def merge_dicts(d1, d2):
    """
    Merges two dictionaries by adding elements from d2 to d1 starting from the length of d1.

    Args:
    d1 (dict): The base dictionary to which elements from d2 will be added.
    d2 (dict): The dictionary whose elements will be added to d1.

    Returns:
    dict: A merged dictionary containing elements from both d1 and d2.
    """
    if isinstance(d1,dict) and isinstance(d2,dict):
        merged = d1.copy()
    
        # Enumerate through the range of length from d1 to d1 + length of d2
        # and add elements from d2 to d1
        for cc, ii in enumerate(np.arange(len(merged), len(merged) + len(d2))):
            merged[ii] = d2[cc]
        return merged
    if isinstance(d1,dict):
        return d1
    if isinstance(d2,dict):
        return d2

def combine_tide_struct_sub(tidestruct,obligatory_fields):
    """
    Combines entries within tidestruct based on identical obligatory fields.
    
    Args:
    tidestruct (list): A list of dictionaries containing tide structure information.
    obligatory_fields (set): A set of obligatory field names that must be identical for merging.
    
    Returns:
    dict: A dictionary containing combined tide structures based on identical obligatory fields.
    """
    result=np.zeros((len(tidestruct),len(tidestruct)))
    for ii in range(len(tidestruct)):
        for jj in range(ii+1,len(tidestruct)):
            # compares fields and returns 1 where the obligatory fields are identical
            result[ii,jj]=compare_fields(tidestruct[ii],tidestruct[jj],obligatory_fields)
    # condense entries that have identical obligatory fields
    ii,jj=np.where(result==1)
    # makes list of index pairs of entries with the same obligatory fields
    ind_pairs = [[ii[kk],jj[kk]] for kk in range(len(ii))]
    # combines all indexes that have similar structure
    combos=merge_lists(ind_pairs)
    combinedStruct={}
    for kk in range(len(combos)):
        # have you already been combined
        for field in obligatory_fields:
            print(field)
            combinedStruct[kk][field] = tidestruct[combos[kk][0]][field]
        # if info field exists, add it
        if 'Info' in tidestruct[combos[kk][0]].keys():
            combinedStruct[kk]['Info']=tidestruct[combos[kk][0]]['Info']
        # other fields common to the two dicts
        other_fields=list(set(tidestruct[ii[kk]]).intersection(tidestruct[jj[kk]])-set(obligatory_fields)-set({'Info'}))
        for field in other_fields:
            init=True
            for inds in combos[kk]: 
                if init:
                    combinedStruct[kk][field]=tidestruct[inds][field] 
                    init=False
                else:
                    combinedStruct[kk][field]+=tidestruct[inds][field]
    return combinedStruct
    
def combine_tide_struct(*args,flag='s'):
    
    """
    Combines multiple load tide structures based on a set of obligatory_fields.

    Call: output = combineTideStruct(loadStruct1, loadStruct2, ...)
    
    Concatenates all tide structures into a single tide structure array.
    After concatenation, it looks by type/Darwin pairs to combine (sum) tide
    structures.
    
    Sets the combine flag. If flag is 's' (default), then all types are kept
    separate. If the flag is 'c' (combine), then it will look for load and
    earth tides, and sum them together changing the type to combined. The
    separate earth and load are still kept in their original forms, but all
    previous combined tides are thrown out.
    
    Output: A single (Nx1) load structure, one element for each unique Darwin
    Symbol/type pair.
    
    Written 1/2012 by Charlie Sievers @ UNAVCO
    """
    
    # Check that there is an input
    if len(args) < 1:
        raise ValueError('combine_tide_struct expects an input.')
    
    # check that the flag is a string and exists and is either s or c
    if not isinstance(flag, str) or len(flag) != 1 or sum(c in 'sc' for c in flag) == 0:
        print(flag)
        raise ValueError('The flag must be either s (separate) or c (combine).')
    
    # Check that at least one tide strucutre is passed
    if len(args) < 1:
        raise ValueError('No tidal structure passed to combineTideStruct')
    
    # Check to make sure there are only Earth, Load, or Combinded tide types
    good_types = ['EarthTide', 'LoadTide', 'CombinedTide']
    for arg in args:
        for jj in range(len(arg)):
            if arg[jj]['Type'] not in good_types:
                raise ValueError('combine_tide_struct is only valid for Earth, Load, and Combined Tides')
    
    # Error checking inputs
    for ii in range(len(args)):
        # check all args are dicts
        if not isinstance(args[ii], dict):
            print(args[ii])
            raise ValueError('All inputs (other than the flag) must be dictionaries.')
        # check that all entries in tide struct have the same keys
        if ii == 0:
            fnames = list(args[ii][0].keys())
        else:
            fnames2 = list(args[ii][0].keys())
            if set(fnames) != set(fnames2):
                raise ValueError('All input structures must be the same type.')
                
    # Combine all inputs into one array structure
    input_struct = {}
    cc=0
    for struct in args:
        for jj in range(len(struct)):
            input_struct[cc]=struct[jj]
            cc+=1
            
    # Remove possible empty load structures
    for idx in range(len(input_struct)):
        if not input_struct[idx]['Station']:
            input_struct.pop(idx)  
    
    # Fields that throw an error if not found
    obligatory_fields = {'Station', 'Type', 'Phase', 'Darwin', 'DoodsonNum'}
    
    # Check if mandatory fields are present 
    # Exit if they arent present
    # Determine remaining nonmandantory fields
    for ii in range(len(input_struct)):
        result = check_mandatory_fields(input_struct[0], obligatory_fields)
        if not result:
            return
        else:
            fnames=check_other_fields(input_struct[0], obligatory_fields)
    
    # Sets the combine flag.  If flag is 's' (default), then all types are kept 
    # separate.  If the flag is 'c' (combine) then it will look for load and
    # earth tides, and sum them together changing the type to combined.  The
    # seperate earth and load are still kept in their original forms, but all
    # previous combined tides are thrown out.
    types = [input_struct[ii]['Type'] for ii in range(len(input_struct))]
    ind_load = [idx for idx, t in enumerate(types) if t == 'LoadTide']
    ind_earth = [idx for idx, t in enumerate(types) if t == 'EarthTide']
    ind_combined = [idx for idx, t in enumerate(types) if t == 'CombinedTide']  
    # how many of each type are there
    n_ind = [len(ind_load), len(ind_earth), len(ind_combined)]
    load_tide, earth_tide, combined_tide=[],[],[]
    # for the load tides
    if n_ind[0] > 0:
        load_tide_dict={}
        for ii,idx in enumerate(ind_load):
            load_tide_dict[ii]=input_struct[idx] 
        len(load_tide_dict)
        load_tide = combine_tide_struct_sub(load_tide_dict, obligatory_fields)
    # for the earth tides
    if n_ind[1] > 0:
        earth_tide_dict = {}
        for ii,idx in enumerate(ind_earth):
            earth_tide_dict[ii]=input_struct[idx]             
        earth_tide = combine_tide_struct_sub(earth_tide_dict, obligatory_fields)
    # for the combined tides
    if n_ind[2] > 0:
        combined_tide_dict = {}
        for ii,idx in enumerate(ind_combined):
            combined_tide_dict[ii]=input_struct[idx] 
        combined_tide = combine_tide_struct_sub(combined_tide_dict, obligatory_fields)
            
    # if you want to keep the types separate
    if flag == 's':
        # merged dictionaries with different types
        output=merge_dict(load_tide,merge_dict(earth_tide,combined_tide))
    # else:  # (flag == 'c') # if you want to combine them

    return output    

out=combine_tide_struct(loadTide)