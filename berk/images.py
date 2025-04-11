"""

This module contains tools for extracting info from MeerKAT images

"""

import os
import sys
import numpy as np
import astropy.io.fits as pyfits
import astropy.stats as apyStats
from astLib import *
from . import plotSettings
import matplotlib.pyplot as plt
from astropy.wcs import WCS

#------------------------------------------------------------------------------------
def getImagesStats(imgFileName, radiusArcmin = 12):
    """Read the given MeerKAT image and return stats such as the image centre coords,
       effective frequency (GHz), RMS in uJy/beam, sky area in sq. deg. etc.

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

    # calculating area
    radiusRA = abs(wcs.header['NAXIS1']*wcs.header['CDELT1']*0.5)
    radiusDec = abs(wcs.header['NAXIS2']*wcs.header['CDELT2']*0.5)
    skyAreaSqDeg = np.pi*radiusRA*radiusDec
    
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
               'skyArea_sqDeg': skyAreaSqDeg,
               'RMS_uJy/beam': sigma*1e6,
               'dynamicRange': d.max()/sigma,
               'freqGHz': wcs.header['CRVAL3']/1e9}

    return statsDict

#------------------------------------------------------------------------------------------
def plotImages(imgFilePath, outDirName, colorMap = 'viridis', vmin = -2.e-5, vmax = 2.e-4,
               axLabelDeg = False, showGrid=True, statsDict = None, overwrite=False):
    """Read the given MeerKAT image and write an output plot of it in PNG format.

    Args:
        imgFilePath (:obj:`str`): Path to the FITS image.
        outDirName (:obj:`str`): Path to the output directory where png files are to be saved.
        colorMap (:obj:`str`, optional): The colormap to use for the image. Default is 'viridis'.
        vmin (:obj:`float`, optional): Minimum data value to anchor the colormap. Default is -2.e-5.
        vmax (:obj:`float`, optional): Maximum data value to anchor the colormap. Default is 2.e-4.
        axLabelDeg (:obj: `bool`, optional): Whether to label the axis coordinates in the units of degrees.
            Default is False.
        showGrid (:obj: `bool`, optional): Whether to show grids. Default is True.
        statsDict (:obj:`dict`, optional): Dictionary containing image statistics to overlay on the plot.
            If provided, a small textbox with key statistics will be displayed on the image. 
            Expected keys:
                - 'freqGHz' (:obj:`float`): Central frequency of the observation in GHz.
                - 'skyArea_sqDeg' (:obj:`float`): Sky area covered by the image in square degrees.
                - 'RMS_uJy/beam' (:obj:`float`): Image RMS noise in microJy/beam.
                - 'dynamicRange' (:obj:`float`): Dynamic range of the image, typically peak/RMS.
        overwrite (:obj: `bool`, optional): Whether to replace the plot, if it exists. Default is False.
    
    Returns:
        None

    """
    
    imgFileName = imgFilePath.split(os.path.sep)[-1].replace(".fits", "")
    imgOutName = outDirName+os.path.sep+imgFileName+".png"
    
    # Skip plotting if file exists and overwrite is not allowed
    if not overwrite and os.path.exists(imgOutName):
        return
    
    with pyfits.open(imgFilePath) as img:
        imageData=img[0].data
        imageHeader=img[0].header
        if imageData.ndim == 4:
            imageData=imageData[0, 0]
        assert(imageData.ndim == 2)
        
    imageData = imageData * 1e6 # converting from Jy/beam to microJy/beam
    
    # finding vmin and vmax
    imageDataClean = np.nan_to_num(imageData, nan=-99., posinf=-99., neginf=-99.)
    vmin = 0.0 #np.percentile(imageDataClean, 5)
    vmax = np.percentile(imageDataClean, 95)
    
    wcs = WCS(imageHeader, naxis=2)
    
    plt.figure(figsize=(8, 6))
    ax = plt.subplot(projection=wcs)
    
    im = ax.imshow(imageData, cmap=colorMap, vmin=vmin, vmax=vmax, origin='lower')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(label=r'$\mu$Jy/beam')
    
    plt.title(imgFileName, fontsize=9)
    plt.xlabel("RA (J2000)")
    plt.ylabel("Dec (J2000)")
    
    if axLabelDeg:
        lon = ax.coords[0]
        lat = ax.coords[1]
        lon.set_major_formatter('d.dd')
        lat.set_major_formatter('d.dd')
        
    if statsDict:
        text = (
        f"Freq: {statsDict['freqGHz']:.2f} GHz\n"
        f"Area: {statsDict['skyArea_sqDeg']:.2f} sq. deg.\n"
        f"RMS: {statsDict['RMS_uJy/beam']:.2f} $\mu$Jy/beam\n"
        f"Dyn. Ran.: {statsDict['dynamicRange']:.2f}"
        )
        
        plt.gca().text(0.02, 0.98, text, fontsize=8, transform=plt.gca().transAxes, ha='left', va='top',  bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.5'))
        
    if showGrid:
        plt.grid(color='white', linestyle='--', linewidth=0.5)

    plt.savefig(imgOutName, dpi=300, bbox_inches = 'tight')
    plt.close()
    
    


