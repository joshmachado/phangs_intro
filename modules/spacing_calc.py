import numpy as np
from astropy.io import fits
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.table import Table
import pandas as pd
import math
import random as r
import scipy.stats
import os.path
import sys

###PATH WHERE DEPROJECT MODULE IS STORED
sys.path.append('/home/machado.35/projects/modules')
from deprojectGalaxy import deproject




def arraySort(variable, distance):
    # sorts variable and distance list by shortest distance
    # variable, distance: lists holding measurement of variable and distance to that measurement
    pattern = distance.argsort()
    dist = distance[pattern]
    var = variable[pattern]
    return (var, dist)

def findNearest(varArray, value, distArray):
    # sorts variable by distance, returning the closest that has the given value
    # varArray, distArray: arrays holding measurement of variable and distance to that measurement
    # value: float.
    var, dist = arraySort(varArray, distArray)
    ind = np.where(var >= value)

    if len(dist[ind]) > 0:
        nearestVal  = np.argmin(dist[ind])
        nearestDist = dist[ind][nearestVal]
        varVal     = var[ind][nearestVal]

    else:
        varVal     = float('nan')
        nearestDist = float('nan')

    return(varVal, nearestDist)

def printNearest(inten, SNR, dist_kpc, value, SNRcutoff = 0.0):
    # returns distance to the nearest molecular cloud and the intensity value found

    #apply SNR cutoff
    intenCut, dist_kpcCut = [],[]
    for i in range(len(inten)):
        if (SNR[i] >= SNRcutoff):
            intenCut.append(inten[i])
            dist_kpcCut.append(dist_kpc[i])

    valFound, nearestMC = findNearest(np.array(intenCut), value, np.array(dist_kpcCut))

    return(nearestMC, valFound)

def distanceCalculator(x1, x2, y1, y2, galDist):
    #calculate distance between two points (in kpc)
    #x1, y1 = xprime and yprime, x2, y2 = SN coords, dist = distance to galaxy (kpc)
    d = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    x = galDist * np.tan(d*np.pi/180.0)
    return(x)

def normalize(weightsArray):
    # takes a list of weights and normalizes them
    prob_factor = 1 / sum(weightsArray)
    return [prob_factor * p for p in weightsArray]

def place_nclouds_in_map(image, n_clouds, outfile_name, galRA, galDEC, pa, incl, galDist):

    hdulist   = pyfits.open(image)
    intMap    = hdulist[0].data

    wcs      = WCS(hdulist[0].header, naxis=2)
    naxis    = wcs._naxis # size of image naxis[0] = x and [1] = y
    grid     = np.indices((naxis[1],naxis[0]))
    ra, dec  = wcs.wcs_pix2world(grid[1],grid[0],0)

    centerCoord = [galRA, galDEC]

    #deproject ra and dec to dx and dy
    radius, projang, dx, dy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=ra, dec=dec, return_offset=True)


    #flatten data structures
    f_int  = intMap.flatten()
    f_ra   = ra.flatten()
    f_dec  = dec.flatten()
    f_dx   = dx.flatten()
    f_dy   = dy.flatten()

    #remove nans
    keep  = np.where(np.isfinite(f_int))
    inten = f_int[keep]
    ra    = f_ra[keep]
    dec   = f_dec[keep]
    dx    = f_dx[keep]
    dy    = f_dy[keep]

    if np.any(inten < 0):
        inten[np.where(inten < 0)] = 0

    #Convert image into probabilities
    total = sum(inten)
    prob  = inten/total

    #Generates index for each coordinate pair
    nX = len(dx)
    indicies = np.arange(nX, dtype=int)


    #Probability weighted random placement
    randInt = np.random.choice(indicies, p=prob, size=int(n_clouds))

    rX   = dx[randInt]
    rY   = dy[randInt]
    rRA  = ra[randInt]
    rDEC = dec[randInt]

    df = pd.DataFrame()
    df['rX'] = rX
    df['rY'] = rY
    df['rRA'] = rRA
    df['rDEC'] = rDEC
    df.to_csv(outfile_name)

    return()

def pair_dist_calc(x,y, neighbor_num):
    n = len(x)
    x1 = np.tile(x, (n,1)) #constant x along column
    y1 = np.tile(y, (n,1)) #constant y along column
    x2 = np.transpose(x1) #constant x along rows
    y2 = np.transpose(y1) #constant y along rows
    dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)

    np.ndarray.sort(dist, axis=1)
    if n < neighbor_num:
        mindist = np.zeros(n)*np.nan
    else:
        mindist = dist[:,neighbor_num]

    return(mindist)

def rad_density(x,y,rad,galdist):
    n = len(x)
    x1 = np.tile(x, (n,1)) #constant x along column
    y1 = np.tile(y, (n,1)) #constant y along column
    x2 = np.transpose(x1) #constant x along rows
    y2 = np.transpose(y1) #constant y along rows
    dist = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    dist = np.deg2rad(dist)*galdist
    inrad = dist < rad
    counts = np.sum(inrad, axis=0) - 1 #subtract one because each cloud has 0 distance to itself and will be counted
    return(counts)
