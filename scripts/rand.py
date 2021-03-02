# import packages and initialize settings

import numpy as np
from astropy.io import fits
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.table import Table

# This block of code contains the function to read fits file
# and deproject galaxy using Adam's deproject.py file


import sys
sys.path.append('/Users/josh/projects/scripts')
from deprojectGalaxy import deproject

def deprojectMap(image, errImage, galRA, galDEC, pa, incl, galDist):
    # convert pixel map to x' & y' and ra & dec
    # image and errImage are the image and noise file pathways
    # galRA, galDEC, pa, incl, galDist are all found within the PHANGS datatabels

    hdulist   = pyfits.open(image)
    intMap    = hdulist[0].data
    hdulist2  = pyfits.open(errImage)
    errMap    = hdulist2[0].data

    wcs      = WCS(hdulist[0].header, naxis=2)
    naxis    = wcs._naxis # size of image naxis[0] = x and [1] = y
    grid     = np.indices((naxis[1],naxis[0]))
    ra, dec  = wcs.wcs_pix2world(grid[1],grid[0],0)

    centerCoord = [galRA, galDEC]

    #deproject ra and dec to dx and dy
    radius, projang, dx, dy = deproject(center_coord=centerCoord, incl=incl, pa=pa, ra=ra, dec=dec, return_offset=True)

    # show galaxy image
    # plt.imshow(map, origin='lower', interpolation='nearest', zorder=1)

    #flatten data structures
    f_int  = intMap.flatten()
    f_err  = errMap.flatten()
    f_ra   = ra.flatten()
    f_dec  = dec.flatten()
    f_dx   = dx.flatten()
    f_dy   = dy.flatten()

    #remove nans
    keep  = np.where(np.isfinite(f_int))
    inten = f_int[keep]
    err   = f_err[keep]
    ra    = f_ra[keep]
    dec   = f_dec[keep]
    dx    = f_dx[keep]
    dy    = f_dy[keep]

    # calculate SNR
    SNR = []
    for i in range(len(inten)):
        if err[i] == 0.0:
            SNR.append(0.0)
        elif inten[i] < 0.0:
            SNR.append(0.0)
        else:
            SNR.append(inten[i]/err[i])

    return(inten, err, SNR, ra, dec, dx, dy)


# This block of code contains the supporting functions to generate
# a random distribution, you should not need to call any of these
# funcitons specifically.

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

def printNearest(inten, SNR, dist_kpc, value, SNRcutoff = 3.0):
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

# This block of code contains the random function,
# which performs your random placement experiment using the
# support functions above.

def random(numTests, weightsArray, dx, dy, ra, dec, galDist):

    # takes number of tests, array of weights, deprojected and
    #       official coordinates, and distance to galaxy...and returns
    #       list of numTests # of random values, their nearest distance,
    #       and their randomly generated coordinates
    # numTests: int , number of randomly generated values
    # weightsArray: list of weights
    # dx&dy:  deprojected coordinates
    # ra&dec: official coordinates
    # galDist: distance to galaxy in kpc
    val, dist, randX, randY, randRA, randDEC = [],[],[],[],[],[]

    for i in range(numTests):

        total = sum(weightsArray)
        prob  = weightsArray/total
        prob  = normalize(prob)

        nX = len(dx)
        indicies = np.arange(nX, dtype=int)
        randInt = np.random.choice(indicies, p=prob)

        rX   = dx[randInt]
        rY   = dy[randInt]
        rRA  = ra[randInt]
        rDEC = dec[randInt]

        distRand = []
        distRand = distanceCalculator(dx, rX, dy, rY, galDist)

        # set sorting cutoff
        intenCutoff = 0

        distance, value = printNearest(inten, SNR, distRand, intenCutoff, SNRcutoff = 3.0)
        val.append(value)
        dist.append(distance)
        randX.append(rX)
        randY.append(rY)
        randRA.append(rRA)
        randDEC.append(rDEC)
    return(val, dist, randX, randY, randRA, randDEC)
