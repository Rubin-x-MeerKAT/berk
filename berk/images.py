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
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm
from astropy.wcs import WCS


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
               'dynamicRange': d.max()/sigma,
               'freqGHz': wcs.header['CRVAL3']/1e9}

    return statsDict
    
def plotImages(imgFilePath, outDirName, colorMap = 'viridis', vmin=-2.e-5, vmax=2.e-4, ax_label_deg=False, show_grid=True):
    """Read the given MeerKAT image and plots it in png format.

    Args:
        imgFilePath (:obj:`str`): Path to the FITS image.
        outDirName (:obj:'str'): Path to the output directory where png files are to be saved.
        colorMap (:obj:'str', optional): The colormap to use for the image. Default is 'viridis'.
        vmin (:obj:'float', optional): Minimum data value to anchor the colormap. Default is -2.e-5.
        vmax (:obj:'float', optional): Maximum data value to anchor the colormap. Default is 2.e-4.
        ax_label_deg (:obj: 'bool', optional): Whether to label the axis coordinates in the units of degrees. Default is False.
        show_grid (:obj: 'bool', optional): Whether to show grids. Default is True.
    
    Returns:
        None
    """
    
    imgFileName = imgFilePath.split(os.path.sep)[-1].replace(".fits", "")
    imgOutName = outDirName+os.path.sep+imgFileName+".png"
    
    if os.path.exists(imgOutName) == True:
        return
    
    captID = imgFileName.split('_')[3]
    target = imgFileName.split('_')[7][:-3]
 
    with pyfits.open(imgFilePath) as img:
        image_data=img[0].data
        image_header=img[0].header
        if image_data.ndim == 4:
            image_data=image_data[0, 0]
        assert(image_data.ndim == 2)

    image_data = image_data * 1e6 # converting from Jy/beam to microJy/beam
    
    # finding vmin and vmax
    image_data_clean = np.nan_to_num(image_data, nan=-99., posinf=-99., neginf=-99.)
    vmin = 0.0 #np.percentile(image_data_clean, 5)
    vmax = np.percentile(image_data_clean, 95)
    
    wcs = WCS(image_header, naxis=2)
    
    plt.figure(figsize=(8, 6))
    ax = plt.subplot(projection=wcs)
    
    norm = simple_norm(image_data, 'log')
    im = ax.imshow(image_data, cmap='viridis', vmin=vmin, vmax=vmax, origin='lower')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(label=r'$\mu$Jy/beam')
    
    plt.title("Capture Block: %s  Target: %s" %(captID, target))
    plt.xlabel("RA (J2000)")
    plt.ylabel("Dec (J2000)")
    
    if(ax_label_deg == True):
        lon = ax.coords[0]
        lat = ax.coords[1]
	
        lon.set_major_formatter('d.dd')
        lat.set_major_formatter('d.dd')

    if(show_grid):
        plt.grid(color='white', linestyle='--', linewidth=0.5)

    pngFileName = imgFileName.split('_')[3]
	
    plt.savefig(imgOutName , dpi=300, bbox_inches = 'tight')
    plt.close()

    
    return None
    

