"""

This module contains tools for extracting info from MeerKAT images

"""

import os
import sys
import numpy as np
import astropy.io.fits as pyfits
import astropy.stats as apyStats
from astLib import *
from . import catalogs

#------------------------------------------------------------------------------------------------------------
def getImagesStats(imgFileName, radiusArcmin = 12):
    """Read the given MeerKAT image and return stats such as the image centre coords,
       effective frequency (GHz), RMS in uJy/beam, etc.

    Args:
        imgFileName (:obj:`str`): Path to the FITS images.
        radiusArcmin (:obj:`float`, optional): Radius in arcmin within which stats will
            be calculated.

    Returns:
        Dictionary of image statistics.

    """

    with pyfits.open(imgFileName) as img:
        d=img[0].data
        if d.ndim == 4:
            d=d[0, 0]
        assert(d.ndim == 2)
        wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
    RADeg, decDeg=wcs.getCentreWCSCoords()
    RAMin, RAMax, decMin, decMax=astCoords.calcRADecSearchBox(RADeg, decDeg, radiusArcmin/60)
    clip=astImages.clipUsingRADecCoords(d, wcs, RAMin, RAMax, decMin, decMax)
    d=clip['data']
    wcs=clip['wcs']
    sigma=1e6
    for i in range(10):
        mask=np.logical_and(np.greater(d, d.mean()-3*sigma), np.less(d, d.mean()+3*sigma))
        sigma=np.std(d[mask])
    # print(">>> Image: %s - radiusArcmin = %.2f" % (sys.argv[1], radiusArcmin))
    # print("    clipped stdev image RMS = %.3f uJy/beam" % (sigma*1e6))
    # sbi=apyStats.biweight_scale(d, c = 9.0, modify_sample_size = True)
    # print("    biweight scale image RMS = %.3f uJy/beam" % (sbi*1e6))
    statsDict={'path': imgFileName,
               'object': wcs.header['OBJECT'],
               'centre_RADeg': RADeg,
               'centre_decDeg': decDeg,
               'RMS_uJy/beam': sigma*1e6,
               'freqGHz': wcs.header['CRVAL3']/1e9}

    return statsDict

