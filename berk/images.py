"""

This module contains tools for extracting info from MeerKAT images

"""

import os
import sys
import numpy as np
import astropy.io.fits as pyfits
import astropy.stats as apyStats
from astLib import *
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd

#------------------------------------------------------------------------------------
def getImageAreaSqDeg(imgFileName):
    """Calculate the total sky area of a given FITS image
    file(image or rms) in square degrees.

    Args:
        imgFileName (:obj:`str`): Path to the FITS image file.

    Returns:
        float: Total sky area covered by the image in square degrees.
    """

    with pyfits.open(imgFileName) as img:
        imgData=img[0].data
        if imgData.ndim == 4:
            imgData=imgData[0, 0]
        assert(imgData.ndim == 2)
        wcs=astWCS.WCS(img[0].header, mode = 'pyfits')

    imgDataFlat = imgData.flatten()
    imgDataFlatNonNan = imgDataFlat[~np.isnan(imgDataFlat)] # ignoring pixels with NaN
    totalNPixels = len(imgDataFlatNonNan) # Number of non Nan pixels

    pixelAreaSqDeg = abs(wcs.header['CDELT1']) * abs(wcs.header['CDELT2']) # area of a pixel in sq. deg.

    skyAreaSqDeg = pixelAreaSqDeg * totalNPixels

    return skyAreaSqDeg

#------------------------------------------------------------------------------------
def findOverlappingPointings(imagesTab, raCol='centre_RADeg', decCol='centre_decDeg', areaCol='skyArea_sqDeg', bandCol='band'):
    """Identify groups of overlapping pointings based on their sky positions
    and approximate sky coverage.

    Args:
        imagesTab (:obj:`astropy.table.Table`): Table containing pointing
            information, including right ascension, declination, and sky area.
        raCol (:obj:`str`, optional): Name of the column containing the
            pointing right ascension in degrees. Defaults to 'centre_RADeg'.
        decCol (:obj:`str`, optional): Name of the column containing the
            pointing declination in degrees. Defaults to 'centre_decDeg'.
        areaCol (:obj:`str`, optional): Name of the column containing the
            sky area covered by each pointing in square degrees. Defaults to
            'skyArea_sqDeg'.
        bandCol (:obj:`str`, optional): Name of the column containing the
            band information for each pointing. Defaults to 'band'.

    Returns:
        tuple: A tuple containing:
            - overlapGroups (:obj:`list`): List of sets, where each set
              contains the indices of pointings belonging to the same
              overlapping group.
            - isolatedIndices (:obj:`list`): List of indices corresponding
              to pointings that do not overlap with any other pointing.
    """
    coords = SkyCoord(ra=imagesTab[raCol]*u.deg, dec=imagesTab[decCol]*u.deg)
    n = len(imagesTab)
    radii = np.array(np.sqrt(imagesTab[areaCol] / np.pi)) # not accurate, but good enough for our purpose

    # Build adjacency: two pointings overlap if their separation < sum of their radii
    overlaps = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i+1, n):
            # Only compare pointings within the same band
            if imagesTab[bandCol][i] != imagesTab[bandCol][j]:
                continue
            sep = coords[i].separation(coords[j]).deg
            if sep < (radii[i] + radii[j]):
                overlaps[i, j] = True
                overlaps[j, i] = True

    # Group overlapping pointings using connected components
    visited = np.zeros(n, dtype=bool)
    overlapGroups = []
    isolatedIndices = []

    for i in range(n):
        if visited[i]:
            continue
        # find all pointings connected to i
        group = set()
        queue = [i]
        while queue:
            current = queue.pop()
            if visited[current]:
                continue
            visited[current] = True
            group.add(current)
            neighbours = np.where(overlaps[current])[0]
            for nb in neighbours:
                if not visited[nb]:
                    queue.append(nb)
        if len(group) == 1:
            isolatedIndices.append(i)
        else:
            overlapGroups.append(group)

    print("\nFound %d overlapping groups and %d isolated pointings\n" % (len(overlapGroups), len(isolatedIndices)))

    overlapGroupsByBand = {}

    for group in overlapGroups:
        band = imagesTab[bandCol][list(group)[0]]

        if band not in overlapGroupsByBand:
            overlapGroupsByBand[band] = []

        overlapGroupsByBand[band].append(group)

    isolatedIndicesByBand = {}

    for idx in isolatedIndices:
        band = imagesTab[bandCol][idx]

        if band not in isolatedIndicesByBand:
            isolatedIndicesByBand[band] = []

        isolatedIndicesByBand[band].append(idx)

    return overlapGroupsByBand, isolatedIndicesByBand

#------------------------------------------------------------------------------------
def getImagesStats(imgFileName, radiusArcmin = 12):
    """Read the given MeerKAT image and return stats such as the image centre coords,
       effective frequency (GHz), RMS in uJy/beam, sky area in sq. deg. etc.

    Args:
        imgFileName (:obj:`str`): Path to the FITS images.
        radiusArcmin (:obj:`float`, optional): Radius in arcmin within which
            RMS and dynamic range will be calculated.

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
    # radiusRA = abs(wcs.header['NAXIS1']*wcs.header['CDELT1']*0.5)
    # radiusDec = abs(wcs.header['NAXIS2']*wcs.header['CDELT2']*0.5)
    # skyAreaSqDeg = np.pi*radiusRA*radiusDec

    skyAreaSqDeg = getImageAreaSqDeg(imgFileName)

    targetObject = wcs.header.get('OBJECT', 'Unspecified')

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
               'object': targetObject,
               'centre_RADeg': RADeg,
               'centre_decDeg': decDeg,
               'skyArea_sqDeg': skyAreaSqDeg,
               'RMS_uJy/beam': sigma*1e6,
               'dynamicRange': d.max()/sigma,
               'freqGHz': wcs.header['CRVAL3']/1e9} #TODO: frequncy header changes during DDFacet

    return statsDict

#------------------------------------------------------------------------------------------
def plotImages(imgFilePath, outDirName=os.getcwd(), colorMap = 'viridis', vmin = None, vmax = None,
               axLabelDeg = False, showGrid=True, statsDict = None, plotTitle=None, overwrite=False):
    """Read the given MeerKAT image and write an output plot of it in PNG format.

    Args:
        imgFilePath (:obj:`str`): Path to the FITS image.
        outDirName (:obj:`str`): Path to the output directory where png files are to be saved.
                                Default is current working directory
        colorMap (:obj:`str`, optional): The colormap to use for the image. Default is 'viridis'.
        vmin (:obj:`float`, optional): Minimum data value to anchor the colormap. Default is 0.
        vmax (:obj:`float`, optional): Maximum data value to anchor the colormap. Default is 95th percentile.
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
        plotTitle (:obj: `str`, optional): Title of the plot. Default is the image file name.
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

    if 'BUNIT' in imageHeader:
        fluxUnit = imageHeader.get('BUNIT').strip()
    else:
        fluxUnit = None

    if not fluxUnit:
        print("Unit of flux not found in header, assuming it to be Jy/beam")
        fluxUnit = 'Jy/beam'

    fluxUnit = fluxUnit.lower()
    if fluxUnit == 'jy/beam':
        imageData = imageData * 1e6 # converting from Jy/beam to microJy/beam
    elif fluxUnit == 'mjy/beam':
        imageData = imageData * 1e3 # converting from mJy/beam to microJy/beam

    fluxUnit = 'microjy/beam'
    fluxUnitLab = r'$\mu$Jy/beam'

    # finding vmin and vmax
    imageDataClean = np.nan_to_num(imageData, nan=-99., posinf=-99., neginf=-99.)

    if not vmin:
        vmin = 0.0 #np.percentile(imageDataClean, 5)
    if not vmax:
        vmax = np.percentile(imageDataClean, 95)

    wcs = WCS(imageHeader, naxis=2)

    plt.figure(figsize=(8, 6))
    ax = plt.subplot(projection=wcs)

    im = ax.imshow(imageData, cmap=colorMap, vmin=vmin, vmax=vmax, origin='lower')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(label=fluxUnitLab)

    if plotTitle:
        plt.title(plotTitle, fontsize=9)
    else:
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

#------------------------------------------------------------------------------------------
def mosaicRMSMaps(rmsFileList, outputFile, resolutionArcmin=1.0):
    """
    Mosaics a list of RMS maps by reprojecting onto a common WCS grid
    and taking the minimum RMS value per pixel across all overlapping pointings.

    Args:
        rmsFileList (list): List of paths to RMS FITS files.
        outputFile (str): Path to save the output mosaicked RMS FITS file.
        resolutionArcmin (float): Output pixel scale in arcminutes. Default is 1.0
            arcmin, which is sufficient for area calculations and avoids memory issues.
            The native pixel scale of MeerKAT images (~2 arcsec) is much finer but
            unnecessary for footprint estimation.

    Returns:
        str: Path to the output mosaicked RMS FITS file.
    """

    hdus = []
    for rmsFile in rmsFileList:
        with pyfits.open(rmsFile) as hdul:
            data = np.squeeze(hdul[0].data).astype(np.float32)
            header = hdul[0].header

        wcs4d = WCS(header)
        wcs2d = wcs4d.celestial

        # Validate the celestial WCS before adding
        # by checking if the corner pixels project to finite sky coordinates
        ny, nx = data.shape
        corners = wcs2d.pixel_to_world_values(
            [0, nx-1, 0, nx-1],
            [0, 0, ny-1, ny-1]
        )
        if not np.all(np.isfinite(corners)):
            print("WARNING: Skipping %s — WCS projects to NaN sky coordinates" % rmsFile)
            continue

        header2d = wcs2d.to_header()
        header2d['NAXIS']  = 2
        header2d['NAXIS1'] = data.shape[1]
        header2d['NAXIS2'] = data.shape[0]

        hdu2d = pyfits.PrimaryHDU(data=data, header=header2d)
        hdus.append(hdu2d)

    if len(hdus) == 0:
        print("ERROR: No valid HDUs found. Cannot mosaic.")
        return None

    wcsOut, shapeOut = find_optimal_celestial_wcs(
        hdus,
        resolution=resolutionArcmin * u.arcmin  # coarser resolution = much smaller array
    )

    mosaicData, footprint = reproject_and_coadd(
        hdus,
        wcsOut,
        shape_out=shapeOut,
        reproject_function=reproject_interp,
        combine_function='min'
    )

    hduOut = pyfits.PrimaryHDU(data=mosaicData.astype(np.float32),
                             header=wcsOut.to_header())
    hduOut.writeto(outputFile, overwrite=True)

    return outputFile

