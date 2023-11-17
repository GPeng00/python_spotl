#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Tue Oct 31 12:53:05 2023

A python version of the matlab version of SPOTL

@author: amt

"""

import numpy as np
import os
import pickle

def nloadf(station, tideFile, greenFile, modo, **kwargs):
    
    """
    Calculate load tides at a specified station using tide and Greens function data.
    
    :param station: Station structure with fields Lat, Long, and Height. 
                    See/use makeStation for creating this structure.
    :param tideFile: This is either the tide model filename or the loaded tide
                    model structure.  See read_tide_model for more information
                    on the requirements for the text file and the structure
                    information.  It is suggested that you preload the tide
                    model before running nloadf if you are going to call it
                    multiple times using readTideModel.
    :param greenFile: Either the file for the Greens function information or a
                    structure containing that information.  See/use 
                    read_greens_function_file for more information.
    :param modo: A single character flag to set the phase to be local (l) or
                Greenwich (g).  The original spotl code had a map option(m)
                but that option is currently unimplemented.
    :param args: Flags and additional arguments.
    :param kwargs: Key-value pairs for various options.
    
    Key/Flag:   

    Polygon     Key/Value  Value is the polygon file name or polygon
                structure.  The default include/exclude is used.
    
    PolygonInclude      Key/Value.  Same as Polygon, but sets default
                        to include.
    
    PolygonExclude      Key/Value.  Same as Polygon, but sets default
                        to exclude.
    
    ShowPlot        Flag. Tells the program to plot the points used to
                    compute the load.
    
    OverrideLndsea  Flag. Tells the program to not use the landsea database
                          to determine if a point is over land and should
                          be ignored.  Normally the landsea database is
                          used if a fine scale lookup is used
    
    OverrideCourse  Flag. Tells the program to use a fine scale lookup for
                          tidal information regardless of the greens 
                          function default.
    
    GridFactor/AngleGridFactor/RangeGridFactor   Key/Value   These three
                Key/Value pairs change the resolution factor for computing
                the tidal loads.  This is a factor, so a value of two for 
                RangeGridFactor will double the range sample rate.  
                Likewise, setting AngleGridFactor to 2 will double the
                angle sample rate from the default ceil(360*sin(Range)).
                GridFactor is shorthand for setting both the Angle and
                Range grid factors to the same value.  I would be wary of
                the product of the two gridfactors being higher than say 
                100 without using SeekPolygon
    
    SeekPolygon    Flag.   This is a special run mode that seeks out points
                   in a polygon and allows for higher gridfactors without
                   the memory crashing.  The program normally find all the
                   points on the sphere over all ranges and angles, and
                   then does a check to see if they are inside/outside the
                   polygon.  For a high grid factor, this can be
                   prohibitively memory intensive.  SeekPolygon first looks
                   for the range and angle boundaries for the polygon and
                   only samples that area.  Because the points sampled will
                   be different, using SeekPolygon will not produce
                   identical results as without.  Furthermore, an error is
                   thrown if SeekPolygon is used without using a polygon
                   include.
    
    PointSource    Key/Value.  This is a special run mode where a set of
                   lat/long/extents are the only sources used to compute
                   the load tides.  The value is a structure giving those
                   points:
    
                   Value.Station:    A Station structure that includes
                                the fields Lat/Long/Height
                   Value.Area:  (sq km) The effective area weighting in square
                                kilometers. 
                   Value.Depth  (mm) This field gives the over
                                    rides saying to ignore the tidalModel
                                    and uses these value instead.  The
                                    value to this field needs to be an Nx1
                                    structure array with the fields Darwin
                                    and Tide.
    
                   When using the Overrides, the tidal model is still used
                   to determine the current harmonic, matched using the
                   Darwin symbol.  If no Darwin symbol matched is found,
                   the tidal model value is used and a warning is thrown.
                   If multiple point sources are used, then the input
                   structure should become a structure array where each
                   element in the array corresponds to a point source.
    
    Output:        A Nx1 structure where each element in the structure
                   corresponds to a specific frequncy constituent.  The
                   fields are:
 

        
    :return: Output structure containing load tide information.
    :rtype: dict
    
    This function calculates load tides at the specified station using tide model and Greens function data. The input
    parameters allow for various options and flags to customize the calculation. The output is a dictionary with load
    tide information.
    
    Default run parameters:
    - polygon: Specifies a polygon file or structure for area weighting.
    - defaultPoly: Specifies the default polygon mode (include/exclude).
    - overrideLndSea: Flag to override land/sea determination.
    - overrideCourse: Flag to use fine-scale lookup for tidal information.
    - specialRunMode: Special run mode (e.g., SeekPolygon).
    - specialRunParameters: Parameters for special run modes.
    - specialRunOverride: Flag to override special run modes.
    - phiResolutionFactor: Resolution factor for computing tidal loads.
    - delResolutionFactor: Resolution factor for computing tidal loads.
    - showPlot: Flag to plot the points used to compute the load.
    
    Additional functions and implementation details are needed to complete the code.
    
    Example usage:
    >>> station = {'Lat': 45.0, 'Long': -120.0, 'Height': 1000.0}
    >>> tideFile = 'tide_model.txt'
    >>> greenFile = 'greens_function.txt'
    >>> modo = 'g'
    >>> output = nloadf(station, tideFile, greenFile, modo)
    
    Written 1/2012 by Charlie Sievers @ UNAVCO.  Modified 11/2023 by Amanda M. Thomas at the University of Oregon.
    
    """
    
    # Read tide model and Greens function
    tideModel = read_tide_model_file(tideFile)
    green = read_greens_function_file(greenFile)
    
    # Check station
    if (not isinstance(station, dict) or any(field not in station for field in ['Lat', 'Long', 'Height'])):
        print(station)
        raise ValueError('Station must have the fields Lat, Long, and Height.')
    
    if modo.lower() not in {'l', 'g'}:
        raise ValueError('Invalid Input: The modo must be local time (l) or Greenwich time (g).')
    
    # Default run parameters
    argStruct = {
        'polygon': '',
        'defaultPoly': '',
        'overrideLndSea': 0,
        'overrideCourse': 0,
        'specialRunMode': 'none',
        'specialRunParameters': {},
        'specialRunOverride': 0,
        'phiResolutionFactor': 1,
        'delResolutionFactor': 1,
        'showPlot': 0
    }
        
    # I don't think this block is needed since function requires 4 args
    # if Nargs == 1 and isinstance(args[0], dict):
    #     # Make sure the structure is of the right type
    #     fnAS = argStruct.keys()
    #     fnVA = args[0].keys()
    #     if len(fnAS) != len(fnVA) or len(fnAS) != sum(key in fnAS for key in fnVA):
    #         raise ValueError('The attempted override for the argument structure is wrong. The structure is not of the right type.')
    
    #     argStruct = args[0]
    #     if len(args) > 1:
    #         print("Attempts to modify an argument structure with key/value and flags are not allowed. All modifications are ignored.")
    # else:
    #     iArg = 0
    #     while iArg < Nargs:
    #         if not isinstance(args[iArg], str):
    #             key = args[iArg]
    #             raise ValueError(f'Invalid key: {key}. Keys must be strings.')
          
    if 'polygon' in kwargs.keys():
        # If polygon input is string, then read the polygon file
        if isinstance(kwargs['polygon'], str):
            argStruct['polygon'] = read_polygon_file(kwargs['polygon'])
        else:
            argStruct['polygon'] = kwargs['polygon']
        kwargs.pop('polygon')
    if 'iore' in kwargs.keys():
        if kwargs['iore'].lower() == 'polygoninclude':
            argStruct['defaultPoly'] = '+'
        elif kwargs['iore'].lower() == 'polygonexclude':
            argStruct['defaultPoly'] = '-'
        kwargs.pop('iore')
    if 'showplot' in kwargs.keys():
        argStruct['showPlot'] = kwargs['showplot']
        kwargs.pop('showplot')
    if 'overridelndsea' in kwargs.keys():
        argStruct['overrideLndSea'] = kwargs['overridelndsea']
        kwargs.pop('overridelndsea')
    if 'overridecourse' in kwargs.keys():
        argStruct['overridecourse'] = kwargs['overridecourse']
        kwargs.pop('overridecourse')
    if 'specialrunmode' in kwargs.keys():
        argStruct['specialrunmode'] = kwargs['specialrunmode']
        kwargs.pop('specialrunmode')
    if 'rangegridfactor' in kwargs.keys():
        argStruct['rangegridfactor'] = kwargs['rangegridfactor']
        kwargs.pop('rangegridfactor')
    if 'anglegridfactor' in kwargs.keys():
        argStruct['phiResolutionFactor'] = kwargs['anglegridfactor'] 
        kwargs.pop('anglegridfactor')
    if 'gridfactor' in kwargs.keys():
        argStruct['phiResolutionFactor'] = kwargs['gridfactor'] 
        argStruct['delResolutionFactor'] = kwargs['gridfactor'] 
    if len(kwargs) > 0:
        raise ValueError("Unknown keyword arguments ", list(kwargs.keys()))


    # include/exclude in function call supercedes that in the polygon file
    if argStruct['defaultPoly']:
        argStruct['polygon']['default'] = argStruct['defaultPoly']
        
    # Input consistency checks
    if argStruct['specialRunMode'] == 'seekPolygon' and argStruct['polygon']['default'].lower() != '+':
        raise ValueError('nloadf attempted to seek points within a polygon, but either no polygon was given, or the default was to exclude.')

    print(argStruct)

    # # TODO: what is specialrunmode
    # # Checks the input structure for point source to make sure it makes sense.
    # if argStruct['specialRunMode'] == 'pointSource':
    #     if len(station['Lat']) > 1:
    #         # Use recursion for point source for multiple stations.
    #         # It would be too confusing otherwise.
    #         N = len(station['Lat'])
    #         out = {}
    #         for ii in range(N):
    #             singleStation = make_station(station, ii)
    #             out[ii] = nloadf(singleStation, tideFile, greenFile, modo, argStruct)
    #         output = out[0]
    #         fields = ['gravLoadTide', 'potentialHeight', 'displacement', 'tilt', 'strain']
    #         for jj in range(len(fields)):
    #             output[fields[jj]] = np.concatenate([out[ii][fields[jj]] for ii in range(N)])
    #         output['Station'] = station
    #         return  # End Recursions
    
    #     sRP = argStruct['specialRunParameters']
    #     if not isinstance(sRP, dict) or not sRP or len(set(['Station', 'Area']) & set(sRP.keys())) != 2:
    #         print(sRP)
    #         # Replace 'keyboard' with an appropriate action
    #         raise ValueError('The Point Source structure appears to be of the wrong form.')
    
    #     try:
    #         sensorTrial = np.concatenate([sRP['Station']])
    #         areaTrial = np.concatenate([sRP['Area']])
    #         latTrial = np.concatenate([sensorTrial['Lat']])
    #         longTrial = np.concatenate([sensorTrial['Long']])
    #         if len(areaTrial) != len(latTrial) or len(areaTrial) != len(longTrial):
    #             raise ValueError('Numel Error')
    #     except Exception as E:
    #         print(E)
    #         raise ValueError('Errors trying to evaluate the point source structure.')
    
    #     if 'Overrides' in sRP:
    #         argStruct['specialRunOverride'] = 1

    # Ok, the program uses recursion to loop all the tide models 
    # (IE K1, O2, etc.) 
    if len(tideModel) > 1:
        outputSub = []
        for ii in range(len(tideModel)):
            outputSub.append(nloadf(station, tideModel[ii], green, modo, argStruct))
        output = combineTideStruct(outputSub)
        return  # End Tidal Model Recursion
    
    '''
    # Parse the green function data
    nGreen = sum([g['ngr'] for g in green])
    greenDel = np.zeros(nGreen)
    greenWidth = np.zeros(nGreen)
    greenData = np.zeros((nGreen, green[0]['ngreen']))
    greenFine = np.zeros(nGreen)

    iGr = 0
    for ii in range(len(green)):
        ind = np.arange(1, green[ii]['ngr'] + 1) + iGr
        greenDel[ind] = np.arange(green[ii]['rin'], green[ii]['rout'] + 1, green[ii]['spc'])
        greenWidth[ind] = green[ii]['spc']
        greenData[ind, :] = green[ii]['data'] / (green[ii]['spc'] * np.pi / 180)
        if green[ii]['fingrd'].lower() == 'f':
            greenFine[ind] = 1
        else:
            greenFine[ind] = 0
        iGr = ind[-1]

    # Sort the arrays
    order = np.argsort(greenDel)
    greenDel = greenDel[order]
    greenWidth = greenWidth[order]
    greenData = greenData[order, :]
    greenFine = greenFine[order]
    
    # Force fine grid lookup
    if argStruct['overrideCourse']:
        greenFine[:] = 1
    
    # Resample the range values in the greens function
    # by a multiplicative factor. Interpolate and Spline
    # the data to match
    if argStruct['delResolutionFactor'] != 1 and argStruct['specialRunMode'] != 'pointSource':
        gLower = greenDel - 0.5 * greenWidth
        gUpper = greenDel + 0.5 * greenWidth
        gBound = np.union1d(gLower, gUpper)
        gBound = np.unique(np.round(1e8 * gBound) / 1e8)
    
        xBoundL = np.linspace(0, 1, len(gBound))
        xBoundU = np.linspace(0, 1, len(gBound) * argStruct['delResolutionFactor'])
        gBound = np.interp(xBoundU, xBoundL, gBound)
    
        greenWidth = np.diff(gBound)
        tmpGreenDel = (gBound[:-1] + gBound[1:]) / 2
        greenData = np.transpose(np.interp(tmpGreenDel, greenDel, np.transpose(greenData)))
        greenFine = np.round(np.interp(tmpGreenDel, greenDel, greenFine))
        greenDel = tmpGreenDel
    
    lat = np.array(station['Lat'])
    long = np.array(station['Long'])
    height = np.array(station['Height'])
    rLat = lat * np.pi / 180
    rLong = long * np.pi / 180

    boundaries_az = [0, 2 * np.pi]
    if argStruct['specialRunMode'].lower() == 'seekpolygon':
        # Find a polygon and explicitly ignore all points outside
        # the polygon. This allows oversampling that region without wasting
        # too much time on the outside points.
        bDel, bPhi = rThetaBoundry(lat, long, argStruct['polygon'])
        boundaries_az = [np.min(bPhi[:, 0]), np.max(bPhi[:, 1])] * np.pi / 180
        boundaries_del = [np.min(bDel[:, 0]), np.max(bDel[:, 1])]
        # Sanity Check
        if np.diff(boundaries_az) > 2 * np.pi:
            boundaries_az[1] = boundaries_az[0] * 2 * np.pi
        # Remove unneeded Ranges
        ind = np.where((greenDel >= boundaries_del[0]) + (greenDel <= boundaries_del[1]) == 2)[0]
        greenDel = greenDel[ind]
        greenWidth = greenWidth[ind]
        greenData = greenData[ind, :]
        greenFine = greenFine[ind]
    
    N = len(height)
    grav0 = np.zeros(N, dtype=complex)
    # Stations
    ct = np.cos((90 - lat) * np.pi / 180)
    st = np.sin((90 - lat) * np.pi / 180)
    rho = 1.025e3
    gParam = constructGFParams(height, ct)
    
    # Point Sources.
    if argStruct['specialRunMode'] == 'pointSource':
        pSS = argStruct['specialRunParameters']['Station']
        psLat = pSS['Lat'] * np.pi / 180
        psLong = pSS['Long'] * np.pi / 180
    
        if len(lat) != 1:
            raise ValueError('The point source run mode does not work with multiple sensors')
    
        # Find the range and bearing from the Sensor to the point sources
        deltaLat = psLat - rLat
        deltaLong = psLong - rLong
        a = np.sin(deltaLat / 2) ** 2 + np.cos(rLat) * np.cos(psLat) * (np.sin(deltaLong / 2) ** 2)
        range2PSR = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        range2PSD = range2PSR * 180 / np.pi
        bearing2PSR = np.angle(1j * np.sin(deltaLong) * np.cos(psLat) +
                               np.cos(rLat) * np.sin(psLat) - np.sin(rLat) * np.cos(psLat) * np.cos(deltaLong))
    
        # Incorporate the area weighting equally between range and azimuth
        solidAngle = argStruct['specialRunParameters']['Area'] / 5.10072e8
        angleWeighting = np.sqrt(solidAngle / (range2PSR / 2 + range2PSR ** 3 / 48 + range2PSR ** 5 / 3840))
    
        greenData = np.transpose(np.interp(range2PSD, greenDel, np.transpose(greenData)))
        greenDel = range2PSD
        greenWidth = angleWeighting * 180 / np.pi
        greenFine = np.ones_like(greenDel)
    
    # Range variables
    NGreen = len(greenDel)
    del_val = greenDel * (np.pi / 180)
    width_val = greenWidth * (np.pi / 180)
    cd = np.cos(del_val)
    sd = np.sin(del_val)
    
    useLndSea = 1 - argStruct['overrideLndSea']
    
    # Use the zero point if it's below water.
    sel = np.where(height < 0)
    if len(sel) > 0 and argStruct['specialRunMode'] != 'pointSource':
        amp = rho * tideModelLookup(tideModel, lat[sel], long[sel],
                                   fingrd=green.fingrd, polygon=polygon, uselndsea=useLndSea)
        grav0[sel] = 2 * np.pi * amp * integratedGreensFunction(del_val[0], 2 * del_val[0], gParam)
    
    g = np.zeros((N, NGreen), dtype=complex)  # Grav
    t = np.zeros((N, NGreen), dtype=complex)  # Tilt
    pot = np.zeros((N, NGreen), dtype=complex)  # Potential Height.
    
    # Process the greens functions with the sensor location corrections
    for ii in range(NGreen):
        g[:, ii], t[:, ii], pot[:, ii] = integratedGreensFunction(del_val[ii], width_val[ii], gParam)


    g = g + np.tile(np.transpose(greenData[:, 2] * width_val), (N, 1))
    t = t + np.tile(np.transpose(greenData[:, 3] * width_val), (N, 1))
    pot = pot + np.tile(np.transpose(greenData[:, 6] * width_val), (N, 1))
    u = greenData[:, 0] * width_val
    v = greenData[:, 1] * width_val
    ett = greenData[:, 4] * width_val
    ell = greenData[:, 5] * width_val
    elt = np.zeros_like(ett)
    
    # Amanda edit to add cells for larger greens functions
    if green[0]['ngreen'] == 9:
        ezz = greenData[:, 7] * width_val
        etz = greenData[:, 8] * width_val
        elz = np.zeros_like(greenData[:, 8])
    
    # Compute the lat/long for points in the concentric rings for each
    # greens function radius. I am keeping to the original spotl code
    # so I'm using a different number of points for each ring, but
    # Python works efficiently with loops, so I'm vectoring the
    # azimuth computations, even if the number of azimuthal points varies.
    
    # How many azimuthal points for each ring
    nAzWhole = np.floor(360 * sd)
    sel = np.where((del_val < np.pi / 2) + (nAzWhole < 150) == 2)
    nAzWhole[sel] = 150
    
    # Rescale for higher resolution
    nAzWhole = np.ceil(nAzWhole * argStruct['phiResolutionFactor'])
    nAzWhole[nAzWhole == 0] = 1  # Make sure the point on the opposite
    
    # side of the globe doesn't add NaN
    
    # For zooming in.
    nAzNum = np.ceil(nAzWhole * (np.diff(boundaries_az) / (2 * np.pi)))
    maxNAZ = np.max(nAzNum)
    
    # Override for the point source. I'm doing this a little backward
    # so I can feed it into my program more easily.
    if argStruct['specialRunMode'] == 'pointSource':
        stepAzRad = angleWeighting
        azimuthRad = bearing2PSR
        logicalIndexZero = np.zeros_like(azimuthRad, dtype=bool)
        nAzNum = np.ones_like(azimuthRad, dtype=int)
        maxNAZ = 1
    else:  # This is for the general case.
        # Compute the azimuthal values for each ring
        stepAzRad = 2 * np.pi / nAzWhole
        azimuthRad = np.full((NGreen, maxNAZ), np.nan)
        for ii in range(NGreen):
            tmpAz = boundaries_az[0] + ((np.arange(nAzNum[ii]) + 0.5) * stepAzRad[ii])
            azimuthRad[ii, :len(tmpAz)] = tmpAz
        logicalIndexZero = np.isnan(azimuthRad)
        azimuthRad[logicalIndexZero] = 0
    
    # The original spotl had an iterative way of propagating sin and cos.
    # There was no reason to keep it, though, since I took it out of the loop.
    caz = np.cos(azimuthRad)
    saz = np.sin(azimuthRad)

    NpointsCall = NGreen * N
    cb = np.zeros((NpointsCall, maxNAZ))
    rlong = np.zeros((NpointsCall, maxNAZ))
    greenMapping = np.zeros(NpointsCall)
    
    iPointsCall = 0
    for iGreen in range(1, NGreen + 1):
        for jPoint in range(1, N + 1):
            iPointsCall += 1
            greenMapping[iPointsCall - 1] = iGreen
            tmpcb = cd[iGreen - 1] * ct[jPoint - 1] + sd[iGreen - 1] * st[jPoint - 1] * caz[iGreen - 1, 0:nAzNum[iGreen - 1]]
            sel = np.abs(tmpcb) > 1
            tmpcb[sel] = np.sign(tmpcb[sel])
            cb[iPointsCall - 1, 0:tmpcb.size] = tmpcb
            sb = np.sqrt(1 - tmpcb**2)
            sel = sb > 1e-3
            sg = sd[iGreen - 1] * saz[iGreen - 1, sel] / sb[sel]
            cg = (st[jPoint - 1] * cd[iGreen - 1] - sd[iGreen - 1] * ct[jPoint - 1] * caz[iGreen - 1, sel]) / sb[sel]
            rlong[iPointsCall - 1, sel] = long[jPoint - 1] + np.arctan2(sg, cg) * 180 / np.pi
    
    rlat = 90 - np.arccos(cb) * 180 / np.pi
    rlong = np.mod(rlong, 360)
    
    ampLookup = np.zeros(rlat.shape, dtype=complex)
    
    fineCourse = np.zeros(rlat.shape)
    fineCourse[np.where(greenFine[greenMapping])] = 1
    
    selFine = fineCourse == 1
    selCourse = fineCourse == 0
    
    if not (argStruct.specialRunMode == 'pointSource'):
        if np.sum(selFine) > 0:
            ampLookup[selFine] = rho * tideModelLookup(tideModel, rlat[selFine], rlong[selFine], 'fingrd', 'f', 'polygon', argStruct.polygon, 'uselndsea', useLndSea)
        if np.sum(selCourse) > 0:
            ampLookup[selCourse] = rho * tideModelLookup(tideModel, rlat[selCourse], rlong[selCourse], 'fingrd', 'c', 'polygon', argStruct.polygon)
    else:
        ampLookup = argStruct.specialRunParameters.Depth
    
    indexZero = np.where(logicalIndexZero[greenMapping])
    ampLookup[indexZero] = 0
    
    toLocalPhase = long * tideModel.DoodsonNum[0] * np.pi / 180
    toLocal = np.cos(toLocalPhase) + 1j * np.sin(toLocalPhase)
    phase = 'Local'
    if modo == 'g':
        toLocal = np.ones(toLocalPhase.shape) + 0j
        phase = 'Greenwich'
    
    grav = np.conj(toLocal * 1e8 * grav)
    pothi = np.conj(1e3 * toLocal * pothi)
    dispVar = np.conj(-1e3 * np.tile(toLocal, (1, 3)) * np.array([elt, ell, ett]).T)
    dispVar[:, 2] *= -1
    tilt = np.conj(1e9 * np.tile(toLocal, (1, 2)) * np.array([tilt[:, 1], tilt[:, 0]]).T)
    
    strain = np.zeros((strain.shape[0], 3))
    strain[:, 0] = ampLookupCos2 * ett - 2 * ampLookupSinCos * elt + ampLookupSin2 * ell
    strain[:, 1] = ampLookupSin2 * ett + ampLookupCos2 * ell + 2 * ampLookupSinCos * elt
    strain[:, 2] = ampLookupSinCos * (ett - ell)
    
    output = {
        'Station': station,
        'Type': 'LoadTide',
        'Phase': phase,
        'Info': None,
        'Darwin': tideModel.Darwin,
        'DoodsonNum': tideModel.DoodsonNum,
        'gravLoadTide': grav,
        'potentialHeight': pothi,
        'displacement': dispVar,
        'tilt': tilt,
        'strain': strain
    }
    
    yn = ['No', 'Yes']
    output['Info'] = {
        'GreenFile': green[0].filename,
        'TideModel': tideModel.file,
        'PolyFile': argStruct.polygon.file,
        'PInclude': argStruct.polygon.default,
        'Phase': phase,
        'OverrideLndSea': yn[argStruct.overrideLndSea],
        'OverrideCourse': yn[argStruct.overrideCourse],
        'PhiResolution': argStruct.phiResolutionFactor,
        'DelResolution': argStruct.delResolutionFactor,
        'SpecialRunMode': argStruct.specialRunMode,
        'PointSourceData': argStruct.specialRunParameters
    }

    return output
    '''
    return


def integrated_greens_function(dell, stp, g, index=None):
    
    """
    
    Calculate the integrated green value from [del-stp/2, del+stp/2] in concentric rings.

    Parameters:
    - del: Current radius.
    - stp: The step size.
    - g: An object with green function parameters, including G2, em, Gn, eps, eps1, eps2, eps3, and plc.
    - index: Optional array of indices for green function parameters.

    Returns:
    - grav: The integrated gravitational effect.
    - tilt: The integrated tilt effect.
    - pot: The integrated potential height.

    If del is greater than or equal to 0.05, the function uses one approximation,
    and if del is less than 0.05, it uses another.

    For more details on the green function parameters, refer to the 'g' parameter.

    If 'index' is not provided, it defaults to a sequence from 1 to the length of g.eps.

    Note: The function assumes numpy arrays or compatible data structures for input.

    Example usage:
    grav, tilt, pot = integratedGreensFunction(0.1, 0.02, green_parameters)
    
    """

    if index is None:
        index = np.arange(1, len(g.eps) + 1)

    c1 = np.sin(0.25 * stp)
    sthl = 0.5 * stp

    if dell >= 0.05:
        grav = -g.G2[index - 1] * np.cos(dell / 2) * c1 * (1 + g.eps[index - 1] / (np.cos(sthl) - np.cos(dell)))
        tilt = g.em[index - 1] * (-2.0 * np.sin(dell / 2) * c1 + np.log(np.tan((dell + sthl) / 4) / np.tan((dell - sthl) / 4)))
    else:  # del < 0.05
        d1 = dell + stp / 2
        d2 = dell - stp / 2
        grav = -g.Gn*(g.eps1[index - 1]*d1*d1 - 2*g.eps[index - 1])/(g.eps3[index - 1]*np.sqrt(g.eps1[index - 1]*d1*d1 + g.eps2[index - 1])) + g.Gn*(g.eps1[index - 1]*d2*d2-2*g.eps[index - 1])/(g.eps3[index - 1]*np.sqrt(g.eps1[index - 1]*d2*d2+g.eps2[index - 1]))
 #       tilt = g.em[index - 1] * (d2 / np.sqrt(g.eps2[index - 1] + d2 * d2) - d1 / np.sqrt(g.eps2[index - 1] + d1 * d1) + np.log((d1 + np.sqrt(g.eps2[index - 1] + d1 * d1)) / (d2 + np.sqrt(g.eps2[index - 1] + d2 * d2)))
    pot = g.plc[index - 1] * np.cos(dell / 2) * c1

    return grav, tilt, pot

def constructGFParams(height, ct):
    """
    Set up the variables for the Greens function integration lookup.

    Parameters:
    - height: An array of sensor heights.
    - ct: An array of cosines of the sensor latitudes.

    Returns:
    - g: A dictionary containing various green function parameters.

    This function calculates green function parameters based on sensor heights
    and cosines of latitudes. The parameters include:
    - N: The number of sensors.
    - eps: Height-to-radius ratio.
    - eps1: 1 + eps.
    - eps2: eps^2.
    - eps3: 2 * eps1^2.
    - G2: The gravitational constant factor.
    - localg: Local gravitational acceleration.
    - em: Gravitational constant divided by local gravity.
    - a: Earth's radius.
    - plc: The partial height correction parameter.
    - Gn: The gravitational constant.

    The formulas used to compute these parameters are based on relevant physical constants.

    Example usage:
    green_params = constructGFParams(sensor_heights, cos_sensor_latitudes)
    """
    N = len(ct)

    # Partial Height Corrections
    a = 6.371e6  # Radius of the Earth
    Gn = 6.67e-11  # Gravitational Constant

    EpsN = np.zeros((N, 4))
    EpsN[:, 0] = height / a
    EpsN[:, 1] = 1 + EpsN[:, 0]
    EpsN[:, 2] = EpsN[:, 0] ** 2
    EpsN[:, 3] = 2 * EpsN[:, 1] ** 2
    G2 = 2 * Gn / (1 + 1.5 * EpsN[:, 0])
    localg = 9.7803327 * (1 + 0.005279 * ct ** 2) - 3.08e-6 * height
    em = Gn / localg
    plc = 4 * a * em

    g = {
        'N': N,
        'eps': EpsN[:, 0],
        'eps1': EpsN[:, 1],
        'eps2': EpsN[:, 2],
        'eps3': EpsN[:, 3],
        'G2': G2,
        'localg': localg,
        'em': em,
        'a': a,
        'plc': plc,
        'Gn': Gn
    }

    return g

def make_station(*args):
    
    """
    Function makeStation: creates, checks, or indexes a station structure.
    
    A station variable is a structure of arrays with the fields: Name, Lat, Long, and Height.
    Name can either be a string or a list of strings of the same size as Lat. 
    Lat, Long, and Height are Nx1 NumPy arrays that must be the same size.

    Args:
        *args: Variable arguments to create, check, or index a station structure.

    Returns:
        dict: A dictionary representing the station structure with fields Name, Lat, Long, and Height.
    """

    nargPass = 0

    if len(args) == 1 and isinstance(args[0], dict):
        station = args[0]
        fieldsExp = ['Name', 'Lat', 'Long', 'Height']
        fieldsFound = sorted(station.keys())

        if len(fieldsExp) == len(fieldsFound) and all(field in fieldsFound for field in fieldsExp):
            if len(station) == 1:
                return station  # The input was a valid station, just return it.
            else:  # It's valid, but we want to concatenate it.
                Name = station['Name']
                Lat = np.concatenate(station['Lat'])
                Long = np.concatenate(station['Long'])
                Height = np.concatenate(station['Height'])
                station = make_station(Name, Lat, Long, Height)
                return station
        else:
            print(station)
            raise ValueError('The only input into makeStation was a structure, but it was not a station structure')
    
    if len(args) == 2 and isinstance(args[0], dict) and isinstance(args[1], int) and len(args[0]) == 1:
        station = make_station(args[0])
        index = args[1]
        if isinstance(station['Name'], list) and len(station['Name']) == len(station['Lat']):
            station['Name'] = [station['Name'][index]]
        station['Lat'] = station['Lat'][index]
        station['Long'] = station['Long'][index]
        station['Height'] = station['Height'][index]
        return station
    
    if len(args) == 2 and isinstance(args[0], dict) and len(args[0]) > 1:
        stn = args[0]
        make_station(stn)  # Error Check
        arg = args[1]
        
        if isinstance(arg, int):
            if max(arg) <= len(stn):
                ind = arg
            else:
                raise ValueError('Attempting to index station for concatenation, but the index is out of bounds')
        elif isinstance(arg, str):
            if arg.lower() == 'all':
                ind = list(range(1, len(stn) + 1))
            else:
                N = len(arg)
                ind = [i for i in range(1, len(stn) + 1) if str(stn[i - 1]['Name']).startswith(arg, 0, N)]
        elif isinstance(arg, list):
            ind = []
            for item in arg:
                for s in stn:
                    N = len(item)
                    if str(s['Name']).startswith(item, 0, N):
                        ind.append(stn.index(s))
            ind = list(set(ind))
        else:
            raise ValueError('Incorrect Call.')
        
        station = {
            'Name': [str(stn[i - 1]['Name']) for i in ind],
            'Lat': np.concatenate([stn[i - 1]['Lat'] for i in ind]),
            'Long': np.concatenate([stn[i - 1]['Long'] for i in ind]),
            'Height': np.concatenate([stn[i - 1]['Height'] for i in ind])
        }
        return station

    if not nargPass:
        if len(args) != 4:
            raise TypeError("makeStation() takes 4 arguments but received {}.".format(len(args)))

    Name, Lat, Long, Height = args

    # Size Check
    if not np.shape(Lat)==np.shape(Long)==np.shape(Height):
        print(np.shape(Lat), np.shape(Long), np.shape(Height))
        raise ValueError("Attempting to construct a sensor structure, but the sizes aren't consistent.")

    station = {
        'Name': [],
        'Lat': Lat,
        'Long': Long,
        'Height': Height
    }
    station['Name'] = Name
    if isinstance(Name, list) and len(Name) == 1:
        station['Name'] = Name[0]

    return station

def read_polygon_file(file):
    """
    Reads a polygon file and returns the information in a dictionary.

    The polygon file is expected to have the following structure:
    - Line 1: A text descriptor of the polygon file.
    - Line 2: The number of polygons in the file (only the first is read).
    - Line 3: The name of the polygon.
    - Line 4: The number of points to define the polygon.
    - Line 5: Either '+' or '-' to set the default include/exclude.
    - Line 5+N: The N long/lat pairs to define the vertices of the polygon, one vertex per line.

    The output dictionary has the following keys:
    - 'file': The name of the polygon file read.
    - 'name': The description of the polygon found on line 3 of the polygon file.
    - 'default': The include/exclude symbol. Either '+' or '-'. If a blank polygon is returned, this value will be '!'.
    - 'points': A list of (lon, lat) pairs for the vertices.

    Args:
        file (str): The name of the polygon file to read.

    Returns:
        dict: A dictionary containing the polygon information.

    Raises:
        ValueError: If the file is not found or has invalid format.

    Example:
        p = read_polygon_file('polygon.txt')
    """

    # Default empty polygon.
    if not file:
        p = {
            'file': 'None',
            'name': None,
            'default': '!',
            'points': []
        }
        return p

    if isinstance(file, dict):
        p = file
        expected_keys = ['file', 'name', 'default', 'points']
        if set(p.keys()) == set(expected_keys):
            return p
        else:
            raise ValueError('The input structure is not a valid polygon structure.')

    if not isinstance(file, str):
        raise ValueError('The polygon file must be a string.')

    try:
        with open(file, 'r') as fptr:
            # Start reading the file.
            p = {'file': file}
            fptr.readline()  # Remove the descriptor line
            Npoly = int(fptr.readline())
            if Npoly != 1:
                print('WARNING: There were multiple polygons found in the file. Only the first is used.')
            p['name'] = fptr.readline().strip()
            N = int(fptr.readline())
            p['default'] = fptr.readline().strip()
            points = []
            for _ in range(N):
                pLine = fptr.readline().strip()
                value = list(map(float, pLine.split()))
                if len(value) != 2:
                    raise ValueError(f'Invalid long/lat pair line in polygon file {file}: {pLine}')
                points.append(value)
            p['points'] = points
    except FileNotFoundError:
        raise ValueError(f'The polygon file {file} was not found.')

    return p
    
def is_empty_dict(d):
        return all(len(v) == 0 for v in d.values())

def read_tide_model_file(tide_file):
    """
    Reads one or multiple tide files. The extension on tide_file can either be txt, mat, or none.
    The extension will be ignored in any case.

    This program first looks for the existence of a .mat file and will attempt to load that file before reading
    in a .txt file. If there is no .mat file, it will read in the .txt file and save it as a .mat file before
    returning the tide model.

    If you give the full filename (e.g., o2.csr4tr.txt), it will load only that one file. If you remove the
    Darwin label (e.g., csr4tr), then all models with that basename will be loaded.

    If tide_file is a list of filenames, each filename will be recursively fed into tide_file.

    If you input a tide dictionary, this function will check if it's a valid tide structure.
    If it passes, it will return that dictionary. If it fails, an error is raised.

    Args:
        tide_file (str, list, dict): The tide model filename(s) or a dictionary containing tide model data.

    Returns:
        list of dict: A list of dictionaries, where each dictionary represents a tide model section.
        The dictionary contains the same keys as in the original data.

    Raises:
        ValueError: If the file is not found or has an invalid format.

    Example:
        tide_models = read_tide_model_file('o2.csr4tr.txt')
    """

    if isinstance(tide_file, dict):  # Check if it's valid
        tidmod = tide_file
        expected_keys = [
            'Name', 'file', 'Darwin', 'DoodsonNum', 'NSEWBoundry', 'dims', 'tideval'
        ]
        if set(tidmod.keys()) == set(expected_keys):
            return [tidmod]
        else:
            raise ValueError('Invalid Tide Model.')
    
    if isinstance(tide_file, list):
        tide_models = []
        for filename in tide_file:
            tide_models += read_tide_model_file(filename)
        return tide_models

    if not isinstance(tide_file, str):
        raise ValueError('The input in readTideModelFile was not a string.')

    # Get the base name without extension
    base, ext = os.path.splitext(tide_file)
    base = base.strip()

    exist_txt = os.path.exists(base + '.txt')
    exist_pkl = os.path.exists(base + '.pkl')

    # Search to see if a mat file exists. If so, load the mat file instead of reading the text file.
    # If not, read the text file and then save it as a mat file.
    if not (exist_txt or exist_pkl):
        start_directory = os.getcwd()
        os.chdir('../tidmod/')
        file_list = [file_name for file_name in os.listdir() if base in file_name]
        os.chdir(start_directory)
        if len(file_list) == 0:
            raise ValueError(f'The tide model file: {tide_file} could not be found.')
        
        for file_name in file_list:
            tide_models = read_tide_model_file(file_name)
        return tide_models

    if exist_pkl:
        try:
            with open(base+'.pkl', 'rb') as f:
                tidmod = pickle.load(f)
            return tidmod
        except Exception as e:
            print(e)
            raise ValueError(f'The mat file {base}.mat failed to load.')

    txt_file = base + '.txt'
    
    try:
        with open(txt_file, 'r') as fptr:
            tidmod = {
                'Name': '',
                'file': base,
                'Darwin': '',
                'DoodsonNum': [0] * 6,
                'NSEWBoundry': [0.0, 0.0, 0.0, 0.0],
                'dims': [0, 0],
                'tideval': []
            }
            tidmod['Darwin'] = fptr.readline().strip()
            if len(tidmod['Darwin']) != 2:
                print(f'Nonstandard Darwin Symbol ::{tidmod["Darwin"]}:: in file {txt_file}')
            
            doodson_line = fptr.readline()
            doodson_line = doodson_line.split()
            if len(doodson_line) != 6:
                raise ValueError(f'The Doodson number is not the expected 6 integers in file {txt_file}')
            
            tidmod['DoodsonNum'] = [int(val) for val in doodson_line]

            north = fptr.readline()
            north_vals = [float(val) for val in north.split()]
            num_north = north_vals[0] + north_vals[1] / 1000

            south = fptr.readline()
            south_vals = [float(val) for val in south.split()]
            num_south = south_vals[0] + south_vals[1] / 1000

            east = fptr.readline()
            east_vals = [float(val) for val in east.split()]
            num_east = east_vals[0] + east_vals[1] / 1000

            west = fptr.readline()
            west_vals = [float(val) for val in west.split()]
            num_west = west_vals[0] + west_vals[1] / 1000

            if None in [num_north, num_south, num_east, num_west]:
                raise ValueError('Bad/No boundary values in tide file: {txt_file}')

            tidmod['NSEWBoundry'] = [num_north, num_south, num_east, num_west]

            dimensions = fptr.readline().split()
            if len(dimensions) != 2:
                raise ValueError(f'Bad dimensions in tide file {txt_file}')

            tidmod['dims'] = [int(dim) for dim in dimensions]
            tidmod['Name'] = fptr.readline().strip()

            # Read the tidal data. This allows for a variable amount of data per line.
            i_data = 0
            n = tidmod['dims'][0] * tidmod['dims'][1]
            data = np.zeros(2 * n)
            while True:
                line = fptr.readline()
                if not line:
                    break
                line = line.split()
                n_line = len(line)
                line = [float(val) for val in line]
                data[i_data:i_data + n_line] = line
                i_data += n_line

            if 2 * n != i_data:
                raise ValueError(f'The tide file {txt_file} didn''t have the expected number of data points.')

            tide_values = 0.001 * (data[:n] + 1j * data[n:])
            tidmod['tideval'] = tide_values

            # Save the tidmod to a mat file for faster call next time.
            current_dir = os.getcwd()
            os.chdir('../tidmod/')
            with open(base+'.pkl', 'wb') as f:
                pickle.dump(tidmod, f)
            os.chdir(current_dir)

        return tidmod

    except Exception as e:
        raise ValueError(f'Failed to read file: {txt_file}. {e}')

def read_greens_function_file(file):

    """
    Read Green's Function Data from a Text File.

    This function reads Green's function data from a text file, performs
    some consistency checks, and returns the data as a structured dictionary.

    If the input is a dictionary, this function checks if it's a valid Green's
    function structure and returns it. If not, it raises an error.

    The file must be formatted as follows:

    - Line 1: A string description of the Green's function.
    - Line 2: Eight space-separated values: NGreen, Num, Ntot, Ngr, Rin, Rout, Spc, Fingrd
      - NGreen: Number of columns in the data line (typically 7).
      - Num: Current Green's function section.
      - Ntot: Total number of Green's function sections in the file.
      - Ngr: Number of lines of data this section contains.
      - Rin: Inner radius for this section in degrees.
      - Rout: Outer radius for this section in degrees.
      - Spc: Radius step size, (Ngr-1)*Spc == Rout-Rin.
      - Fingrd: Interpolation type ('F' for file scale or 'C' for coarse).
    - Line 3+: The integrated Green's function data with NGreen whitespace-separated entries
      for each line and Ngr lines.

    Inputs:
    - file: The Green's function filename (str).

    Output:
    - green: A structured dictionary containing Green's function data with fields:
      - filename: The name of the file.
      - description: Description of the Green's function.
      - ngreen: Number of columns in the data.
      - num: Current Green's function section.
      - ntot: Total number of sections.
      - ngr: Number of data lines in this section.
      - rin: Inner radius (in degrees).
      - rout: Outer radius (in degrees).
      - spc: Radius step size.
      - fingrd: Interpolation type ('F' for file scale or 'C' for coarse).
      - data: A 2D NumPy array containing the integrated Green's function data.

    Raises:
    - ValueError: If the input structure is invalid or if there are any errors
      in reading and formatting the Green's function data.

    Written by Charlie Sievers @ UNAVCO on 12/11.
    """
    # Function code goes here.
    
    if isinstance(file, dict):  # Check if it's valid
        green = file
        expected_keys = [
        'filename', 'description', 'ngreen', 'num', 'ntot', 'ngr', 'rin', 'rout', 'spc', 'fingrd', 'data'
        ]
        if set(green[0].keys()) == set(expected_keys):
            return [green]
        else:
            raise ValueError('Invalid Green function structure. Wrong fields.')

    if not isinstance(file, str):
        raise ValueError('The Green function file must be a string.')

    try:
        fptr = open(file, "r+")
        green = {}
        ii = -1
        description = fptr.readline()
        while True: 
            line = fptr.readline()
            if not line:  # This is effectively feof is true
                break
            params = line.split()
            if len(params) != 8:
                break
            ii += 1
            green[ii] = {}
            green[ii]['filename'] = file
            green[ii]['description'] = description
            green[ii]['ngreen'] = int(params[0])
            green[ii]['num'] = int(params[1])
            green[ii]['ntot'] = int(params[2])
            green[ii]['ngr'] = int(params[3])
            green[ii]['rin'] = float(params[4])
            green[ii]['rout'] = float(params[5])
            green[ii]['spc'] = float(params[6])
            green[ii]['fingrd'] = params[7].upper()
            data = np.zeros((green[ii]['ngr'], green[ii]['ngreen']))
            try:
                for jj in range(green[ii]['ngr']):
                    line = fptr.readline()
                    data[jj, :] = [float(value) for value in line.split()]
            except:
                raise ValueError(f'The green function file {file} was improperly formatted. '
                                 f'The {jj} data line of the {ii} section was expected to be a '
                                 f'data line, but instead it read {line}.')
            green[ii]['data'] = data
    
    
    except FileNotFoundError:
        raise ValueError(f'Unable to open Green\'s function file: {file}.')
    
    for kk in range(1,ii+1):
        # Consistency check for the Green's function
        expected_keys = ['ngreen', 'num', 'ntot', 'ngr', 'rin', 'rout', 'spc', 'fingrd', 'data']
        for key in expected_keys:
            if key not in green[kk]:
                raise ValueError(f'The field {key} in the data is missing.')
        
        for section in green[kk]['data']:
            M, N = green[kk]['data'].shape  # All the data's there
            if M != green[kk]['ngr'] or N != green[kk]['ngreen']:
                raise ValueError('Consistency check failure. Sizes not compatible.')
        
            # The radii add up correctly
            if abs(green[kk]['rin'] + green[kk]['spc'] * (green[kk]['ngr'] - 1) - green[kk]['rout']) > 1e-6:
                raise ValueError('Consistency check failure. rout - rin != spc * (ngr - 1).')
    
    # All the sections are there
    if len(green)!=green[1]['ntot']:
        raise ValueError("The Green's function file {file} doesn't have the expected number of sections.")

    return green