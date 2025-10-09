"""

This module contains tools for creating the summary plots

"""

import os
import numpy as np
from astropy.io import fits
import astropy.table as atpy
from . import startup, tasks
import matplotlib.pyplot as plt
import glob

def plotSkyCoverage(fullImagesTab, bandColorDict, plotOutPath, plotProjection='aitoff'):
    """Plot the sky coverage of MeerKAT pointings for each band on an Aitoff projection.

    Args:
        fullImagesTab (:obj:`astropy.table.Table`): Table containing image metadata including RA/Dec and band info.
        bandColorDict (:obj:`dict`): Dictionary mapping band names to color codes.
        plotOutPath (:obj:`str`): Path to save the output plot.
        plotProjection (:obj:`str`, optional): Projection to use for plotting (default is 'aitoff').

    Returns:
        None. Saves the sky coverage plot to the specified path.
    """

    plt.figure(figsize=(10, 5))
    ax = plt.subplot(111, projection=plotProjection)

    orderedBands = list(bandColorDict.keys())

    for band in orderedBands:
        bandMask = fullImagesTab['band'] == band
        bandData = fullImagesTab[bandMask]

        # Plotting sky coverage

        raBandDataDeg = np.array(bandData['centre_RADeg'])
        decBandDataDeg = np.array(bandData['centre_decDeg'])

        # The aitoff projection in Matplotlib expects RA values to be in the range of [-pi, pi] (or [-180, 180] degrees) where 0 is at the center.
        raBandDataRad = np.radians(np.remainder(raBandDataDeg + 360 - 180, 360) - 180)
        decBandDataRad = np.radians(decBandDataDeg)

        ax.scatter(raBandDataRad, decBandDataRad, color="black", fc=bandColorDict[band], s=70, alpha=1, label="%s" %band, zorder=5)

    ax.grid(True, color='gray', linestyle='--', linewidth=0.7, zorder=0)

    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")
    plt.legend(loc="upper right")
    plt.savefig(plotOutPath, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nMeerKAT processed pointings plotted!\n")

def _getS2p5dNdS(SMean, NSources, SBinWidth, skyAreaSqDeg):
    """Compute Euclidean-normalized differential source counts.

    Args:
        SMean (:obj:`np.ndarray`): Mean flux per bin.
        NSources (:obj:`np.ndarray`): Number of sources per bin.
        SBinWidth (:obj:`np.ndarray`): Width of each flux bin.
        skyAreaSqDeg (:obj:`float`): Total survey area in square degrees.

    Returns:
        :obj:`np.ndarray`: Euclidean-normalized source count per bin.
    """

    skyAreaSterdian = ((np.pi/180.)**2)*skyAreaSqDeg # sq. deg to steredian
    sourceCount = SMean**2.5 * (NSources/(SBinWidth*skyAreaSterdian))
    return sourceCount

def _computeSourceCount(fluxVals, skyAreaSqDeg, nFluxBins=20):
    """Compute source counts normalized by S^2.5 for log-spaced flux bins.

    Args:
        fluxVals (:obj:`np.ndarray`): Array of flux values.
        skyAreaSqDeg (:obj:`float`): Sky area in square degrees.
        nFluxBins (:obj:`int`, optional): Number of flux bins (default is 20).

    Returns:
        tuple: (bin centers, source count values, source count errors).
    """

    fluxMin = np.min(fluxVals)
    fluxMax = np.max(fluxVals)
    fluxBins = np.logspace(np.log10(fluxMin), np.log10(fluxMax), nFluxBins)

    SCounts, binEdges = np.histogram(fluxVals, bins=fluxBins)

    SMean = (binEdges[:-1] + binEdges[1:]) / 2.
    SBinWidths = np.diff(binEdges)

    sourceCountValues = _getS2p5dNdS(SMean, SCounts, SBinWidths, skyAreaSqDeg)
    sourceCountErr = _getS2p5dNdS(SMean, np.sqrt(SCounts), SBinWidths, skyAreaSqDeg)

    return SMean, sourceCountValues, sourceCountErr

def plotSourceCounts(fullImagesTab, fluxCol, nFluxBins, bandColorDict, plotOutPath, plotMALS=False, plotLOFAR=False):
    """Plot Euclidean-normalized source counts for each band using survey catalogs.

    Args:
        fullImagesTab (:obj:`astropy.table.Table`): Table containing metadata for all images.
        fluxCol (:obj:`str`): Name of the flux column in the catalog table.
        nFluxBins (:obj:`int`): Number of flux bins.
        bandColorDict (:obj:`dict`): Dictionary mapping band names to color codes.
        plotOutPath (:obj:`str`): Path to save the output plot.

    Returns:
        None. Saves the source count plot to the specified path.
    """

    plt.figure(figsize=(8, 5))
    ax = plt.subplot(111)

    orderedBands = list(bandColorDict.keys())

    for band in orderedBands:
        bandMask = fullImagesTab['band'] == band
        bandData = fullImagesTab[bandMask]
        bandTotalArea = bandData['skyArea_sqDeg'].sum()

        catFileName = startup.config['productsDir']+os.path.sep+"survey_catalog_%s.fits" %band
        catalogTab = atpy.Table().read(catFileName)

        fluxVals = catalogTab[fluxCol].value
        fluxUnitLabel = catalogTab[fluxCol].unit

        sourceCountFlux, sourceCountVal, sourceCountErr =  _computeSourceCount(fluxVals=fluxVals, skyAreaSqDeg=bandTotalArea, nFluxBins=nFluxBins)

        ax.errorbar(sourceCountFlux, sourceCountVal, sourceCountErr, mec='k', mfc=bandColorDict[band], ecolor=bandColorDict[band], marker='o', ms=7, alpha=1, ls='None', label="MeerKAT %s-band" %band)

    if plotMALS is True:

        malsData = np.loadtxt(startup.config['productsDir']+os.path.sep+"source_counts_MALS.txt", skiprows=1)

        SmJyMALS = malsData[:, 0]
        SJyMALS = SmJyMALS * 1e-3
        S5dNdSMALS = malsData[:, 1]
        S5dNdSErrMALS = malsData[:, 2]

        ax.errorbar(SJyMALS, S5dNdSMALS, S5dNdSErrMALS, mec='k', mfc='#F0B13B', ecolor='#F0B13B', marker='s', ms=7, alpha=1, ls='None', label="MALS L-band (Wagenveld+23)")

    if plotLOFAR is True:
        lofarData = np.loadtxt(startup.config['productsDir']+os.path.sep+"source_counts_LOFAR.txt", skiprows=1)

        SmJyLOFAR = lofarData[:, 0]
        SJyLOFAR = SmJyLOFAR * 1e-3
        S5dNdSLOFAR = lofarData[:, 1]
        S5dNdSLowerErrLOFAR = lofarData[:, 2]
        S5dNdSUpperErrLOFAR = lofarData[:, 3]

        ax.errorbar(SJyLOFAR, S5dNdSLOFAR, [S5dNdSLowerErrLOFAR, S5dNdSUpperErrLOFAR], mec='k', mfc='#FDA5D5', ecolor='#FDA5D5', marker='s', ms=7, alpha=1, ls='None', label="LOFAR 150 MHz (Williams+16)")



    ax.set_xlabel("Total Flux (%s)" %fluxUnitLabel)
    ax.set_ylabel(r"$S^{5/2} \mathrm{d}N/\mathrm{d}S$ $(\mathrm{Jy}^{3/2} \mathrm{sr}^{-1})$")

    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.legend(loc="upper left")
    plt.savefig(plotOutPath, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nSourcecounts plotted!\n")

def _getRMSCoverage(rmsMapFile):
    """Compute the CDF of RMS noise from a given RMS map.

    Args:
        rmsMapFile (:obj:`str`): Path to the FITS file containing RMS map.

    Returns:
        tuple: (bin centers, CDF values) of the RMS distribution.
    """

    with fits.open(rmsMapFile, memmap=True) as hdu:
        rmsData = hdu[0].data

    rmsData = rmsData.flatten()
    rmsData = rmsData[~np.isnan(rmsData)]

    rmsCount, fluxBinsEdges = np.histogram(rmsData, bins=100)
    pdfRMS = rmsCount / sum(rmsCount)
    cdfRMS = np.cumsum(pdfRMS)
    fluxBinCenter = (fluxBinsEdges[:-1] + fluxBinsEdges[1:]) / 2

    return fluxBinCenter, cdfRMS

def plotRMSCoverage(plotOutPath):
    """Plot RMS noise coverage from all RMS map files in the configured directory.

    Args:
        plotOutPath (:obj:`str`): Path to save the output RMS coverage plot.

    Returns:
        None. Saves the RMS coverage plot to the specified path.
    """

    rmsDirPath = startup.config['productsDir']+os.path.sep+'rms'+os.path.sep+"*rms.fits"

    rmsFiles=glob.glob(rmsDirPath)

    fluxUnitLabel = 'Jy'

    plt.figure(figsize=(8, 5))
    ax = plt.subplot(111)

    for rmsFile in rmsFiles:
        fluxBinCenter, cdfRMS = _getRMSCoverage(rmsFile)
        ax.plot(fluxBinCenter, cdfRMS)

    ax.set_xlabel("RMS Noise (%s/beam)" %fluxUnitLabel)
    ax.set_ylabel("Coverage")

    ax.set_xscale('log')

    plt.savefig(plotOutPath, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nRMS Coverage plotted!\n")

def _getRMSHistogram(rmsFile, rmsBins):
    """Compute histogram and corresponding area for a given RMS map.

    Args:
        rmsFile (:obj:`str`): Path to the FITS file containing RMS map.
        rmsBins (:obj:`np.ndarray`): Array of RMS bin edges.

    Returns:
        tuple: (counts per bin, area per bin in square degrees).
    """

    with fits.open(rmsFile, memmap=True) as hdu:
        data = np.squeeze(hdu[0].data).flatten()
        hdr = hdu[0].header
        pixelAreaSqDeg = abs(hdr['CDELT1']) * abs(hdr['CDELT2'])
        data = data[~np.isnan(data)]
        data = data[data > 0]

    countsInBins, _ = np.histogram(data, bins=rmsBins)
    areaInBins = countsInBins * pixelAreaSqDeg
    return countsInBins, areaInBins

def plotRMSAreaCoverageDifferential(plotOutPath, bandColorDict, nRMSBins=30):
    """Plot RMS noise vs sky area coverage for all bands, using cached histograms if available.

    Args:
        plotOutPath (:obj:`str`): Path to save the output plot.
        bandColorDict (:obj:`dict`): Dictionary mapping band names to color codes.
        nRMSBins (:obj:`int`, optional): Number of RMS bins to use (default is 30).

    Returns:
        None. Saves the RMS area coverage plot to the specified path.
    """

    rmsDirPath = startup.config['productsDir']+os.path.sep+'rms'

    rmsFiles = sorted(glob.glob(rmsDirPath+os.path.sep+"*rms.fits"))

    rmsBins = np.logspace(-6, -2, nRMSBins)

    binCentres = 0.5 * (rmsBins[1:] + rmsBins[:-1])
    binCentresMuJy = binCentres * 1e6
    fluxUnitLabel = r'$\mu$Jy'

    totalCounts = np.zeros(len(binCentres))
    totalAreaSqDeg = np.zeros(len(binCentres))

    orderedBands = list(bandColorDict.keys())
    globalRMSCountDict={'L': np.zeros(len(binCentres)), 'UHF': np.zeros(len(binCentres)), 'S': np.zeros(len(binCentres))}
    globalRMSAreaDict={'L': np.zeros(len(binCentres)), 'UHF': np.zeros(len(binCentres)), 'S': np.zeros(len(binCentres))}

    for rmsFile in rmsFiles:

        rmsFileHeader = fits.getheader(rmsFile)
        bandFreqHz = rmsFileHeader['CRVAL4']
        bandKey = tasks.getBandKey(bandFreqHz*1e-9)

        if bandKey not in globalRMSCountDict.keys():
            print("Skipping unrecognised band:", rmsFile)
            continue

        baseName = os.path.basename(rmsFile).split('.fits')[0]
        rmsHistFile = os.path.join(rmsDirPath, baseName + "_rmshist.txt")

        if os.path.exists(rmsHistFile):
            rmsHistFromFile = np.loadtxt(rmsHistFile)
            binCentresFile, countsFile, areaSqDegFile = rmsHistFromFile [:, 0], rmsHistFromFile [:, 1], rmsHistFromFile [:, 2]

            if len(binCentres) != len(binCentresFile) or not np.allclose(binCentres, binCentresFile):
                countsInBins, areaSqDegInBins = _getRMSHistogram(rmsFile, rmsBins)
                np.savetxt(rmsHistFile, np.column_stack((binCentres, countsInBins, areaSqDegInBins)))
            else:
                countsInBins = countsFile
                areaSqDegInBins = areaSqDegFile

        else:
            countsInBins, areaSqDegInBins = _getRMSHistogram(rmsFile, rmsBins)
            np.savetxt(rmsHistFile, np.column_stack((binCentres, countsInBins, areaSqDegInBins)))

        # Add to total histogram
        totalCounts += countsInBins
        totalAreaSqDeg += areaSqDegInBins

        globalRMSCountDict[bandKey] += countsInBins
        globalRMSAreaDict[bandKey] += areaSqDegInBins

    #plt.figure(figsize=(8, 5))
    #ax = plt.subplot(111)

    fig,ax=plt.subplots(nrows=1,ncols=3,sharex=True, sharey=False)
    fig.set_size_inches(16,5)

    for bandi, band in enumerate(orderedBands):
        positiveMask = globalRMSCountDict[band] > 0

        ax[bandi].plot(binCentresMuJy[positiveMask], globalRMSAreaDict[band][positiveMask], drawstyle='steps-mid', color=bandColorDict[band],  label="%s" %band)

        ax[bandi].set_xlabel("RMS Noise (%s/beam)" %fluxUnitLabel)
        ax[bandi].set_xscale('log')
        ax[bandi].set_yscale('log')

    ax[0].set_ylabel("Area (sq. deg.)")

    plt.legend(loc="upper right")
    plt.savefig(plotOutPath, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nRMS Area coverage plotted!\n")

def plotRMSAreaCoverageCumulative(coveragePlotOutPath, fractionPlotOutPath, bandColorDict, bandTotalAreaDict, nRMSBins=30):
    """Plot cumulative sky area as a function of RMS noise for all bands.

    Args:
        plotOutPath (str): Path to save the output plot.
        bandColorDict (dict): Dictionary mapping band names to color codes.
        nRMSBins (int, optional): Number of RMS bins to use (default is 30).

    Returns:
        None. Saves the cumulative RMS area plot to the specified path.
    """

    rmsDirPath = startup.config['productsDir'] + os.path.sep + 'rms'
    rmsFiles = sorted(glob.glob(rmsDirPath + os.path.sep + "*rms.fits"))

    rmsBins = np.linspace(0, 2e-4, nRMSBins)
    binCentres = 0.5 * (rmsBins[1:] + rmsBins[:-1])
    binCentresMuJy = binCentres * 1e6
    fluxUnitLabel = r'$\mu$Jy'

    orderedBands = list(bandColorDict.keys())
    globalRMSAreaDict = {band: np.zeros(len(binCentres)) for band in orderedBands}

    # Loop over RMS files
    for rmsFile in rmsFiles:
        rmsFileHeader = fits.getheader(rmsFile)
        bandFreqHz = rmsFileHeader['CRVAL4']
        bandKey = tasks.getBandKey(bandFreqHz * 1e-9)

        if bandKey not in globalRMSAreaDict:
            print("Skipping unrecognised band:", rmsFile)
            continue

        baseName = os.path.basename(rmsFile).split('.fits')[0]
        rmsHistFile = os.path.join(rmsDirPath, baseName + "_rmshist.txt")

        if os.path.exists(rmsHistFile):
            rmsHistFromFile = np.loadtxt(rmsHistFile)
            binCentresFile, _, areaSqDegFile = rmsHistFromFile[:, 0], rmsHistFromFile[:, 1], rmsHistFromFile[:, 2]

            if len(binCentres) != len(binCentresFile) or not np.allclose(binCentres, binCentresFile):
                countsInBins, areaSqDegInBins = _getRMSHistogram(rmsFile, rmsBins)
                np.savetxt(rmsHistFile, np.column_stack((binCentres, countsInBins, areaSqDegInBins)))
            else:
                areaSqDegInBins = areaSqDegFile
        else:
            countsInBins, areaSqDegInBins = _getRMSHistogram(rmsFile, rmsBins)
            np.savetxt(rmsHistFile, np.column_stack((binCentres, countsInBins, areaSqDegInBins)))

        globalRMSAreaDict[bandKey] += areaSqDegInBins

    # Plot cumulative area

    fig,ax=plt.subplots(nrows=1,ncols=3,sharex=True, sharey=False)
    fig.set_size_inches(12,3)

    for bandi, band in enumerate(orderedBands):
        cumulativeArea = np.cumsum(globalRMSAreaDict[band])

        positiveMask = cumulativeArea > 0

        ax[bandi].plot(binCentresMuJy[positiveMask], cumulativeArea[positiveMask],
                drawstyle='steps-mid', color=bandColorDict[band])

        ax[bandi].set_xlabel("RMS Noise (%s/beam)" % fluxUnitLabel)
        ax[bandi].set_xlim(0, 200)
        ax[bandi].text(0.5,0.10, "%s band" %(band),transform=ax[bandi].transAxes,ha='center')


    ax[0].set_ylabel("Cumulative Area (sq. deg.)")
    #ax.set_xscale('log')
    #ax.set_yscale('log')

    plt.savefig(coveragePlotOutPath, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nCumulative RMS area plotted!\n")

    # Plot cumulative area fraction

    fig,ax=plt.subplots(nrows=1,ncols=3,sharex=True, sharey=False)
    fig.set_size_inches(12,3)

    for bandi, band in enumerate(orderedBands):
        cumulativeArea = np.cumsum(globalRMSAreaDict[band])

        positiveMask = cumulativeArea > 0

        ax[bandi].plot(binCentresMuJy[positiveMask], cumulativeArea[positiveMask]/bandTotalAreaDict[band],
                drawstyle='steps-mid', color=bandColorDict[band])

        ax[bandi].set_xlabel("RMS Noise (%s/beam)" % fluxUnitLabel)
        #ax[bandi].set_xlim(0, 200)
        ax[bandi].text(0.5,0.10, "%s band" %(band),transform=ax[bandi].transAxes,ha='center')


    ax[0].set_ylabel("Fraction of cumulative Area")
    #ax.set_xscale('log')
    #ax.set_yscale('log')

    plt.savefig(fractionPlotOutPath, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nCumulative RMS area fraction plotted!\n")
