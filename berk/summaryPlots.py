"""

This module contains tools for creating the summary plots

"""

import os
import numpy as np
from astropy.io import fits
import astropy.table as atpy
from . import startup, tasks, images
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

def getS2p5dNdS(SMean, NSources, SBinWidth, skyAreaSqDeg):
    """Compute Euclidean-normalized differential source counts.

    Args:
        SMean (:obj:`np.ndarray`): Mean flux per bin.
        NSources (:obj:`np.ndarray`): Number of sources per bin.
        SBinWidth (:obj:`np.ndarray`): Width of each flux bin.
        skyAreaSqDeg (:obj:`float`): Array of area (sq. deg) covered in each flux bin.

    Returns:
        :obj:`np.ndarray`: Euclidean-normalized source count per bin.
    """
    skyAreaSqDeg = np.array(skyAreaSqDeg)
    skyAreaSterdian = ((np.pi/180.)**2)*skyAreaSqDeg # sq. deg to steredian

    sourceCount = SMean**2.5 * (NSources/(SBinWidth*skyAreaSterdian))
    return sourceCount

def computeSourceCount(fluxVals, skyAreaSqDeg, fluxMin=None, fluxMax=None, nFluxBins=20):
    """Compute source counts normalized by S^2.5 for log-spaced flux bins.

    Args:
        fluxVals (:obj:`np.ndarray`): Array of flux values.
        skyAreaSqDeg (:obj:`float`): Array of area (sq. deg) covered in each flux bin.
        nFluxBins (:obj:`int`, optional): Number of flux bins (default is 20).

    Returns:
        tuple: (bin centers, source count values, source count errors).
    """

    fluxMin = np.min(fluxVals) if fluxMin is None else fluxMin
    fluxMax = np.max(fluxVals) if fluxMax is None else fluxMax
    fluxBins = np.logspace(np.log10(fluxMin), np.log10(fluxMax), nFluxBins+1)

    SCounts, binEdges = np.histogram(fluxVals, bins=fluxBins)

    SMean = (binEdges[:-1] + binEdges[1:]) / 2.
    SBinWidths = np.diff(binEdges)

    sourceCountValues = getS2p5dNdS(SMean, SCounts, SBinWidths, skyAreaSqDeg)
    sourceCountErr = getS2p5dNdS(SMean, np.sqrt(SCounts), SBinWidths, skyAreaSqDeg)

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

        sourceCountFlux, sourceCountVal, sourceCountErr =  computeSourceCount(fluxVals=fluxVals, skyAreaSqDeg=bandTotalArea, nFluxBins=nFluxBins)

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

def getRMSAreaCoverage(rmsFile, rmsBins):
    """Compute histogram and corresponding area for a given RMS map.

    Args:
        rmsFile (:obj:`str`): Path to the FITS file containing RMS map.
        rmsBins (:obj:`np.ndarray`): Array of RMS bin edges.

    Returns:
        tuple: (counts per bin, area per bin in square degrees).
    """

    with fits.open(rmsFile, memmap=True) as hdu:
        dataFull = np.squeeze(hdu[0].data).flatten()
        data = dataFull[~np.isnan(dataFull)]
        hdr = hdu[0].header
        pixelAreaSqDeg = abs(hdr['CDELT1']) * abs(hdr['CDELT2'])
        data = data[~np.isnan(data)]
        data = data[data > 0]

    countsInBins, _ = np.histogram(data, bins=rmsBins)
    areaInBins = countsInBins * pixelAreaSqDeg
    return countsInBins, areaInBins

def plotRMSAreaCoverageCumulative(areaCoveragePlotOutName, bandColorDict, nRMSBins=30):
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

    rmsBins = np.logspace(-10, 0, nRMSBins)
    binCentres = 0.5 * (rmsBins[1:] + rmsBins[:-1])

    orderedBands = list(bandColorDict.keys())
    globalRMSNPixelsInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSCumulativeNPixelsInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSAreaInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSCumulativeAreaInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    totalRMSAreaDict = {band: 0 for band in orderedBands}

    # Loop over RMS files
    for rmsFile in rmsFiles:
        rmsFileHeader = fits.getheader(rmsFile)
        bandFreqHz = rmsFileHeader['CRVAL4']
        bandKey = tasks.getBandKey(bandFreqHz * 1e-9)

        if bandKey not in globalRMSAreaInBinsDict:
            print("Skipping unrecognised band:", rmsFile)
            continue

        rmsFileSkyArea = images.getImageAreaSqDeg(rmsFile)
        totalRMSAreaDict[bandKey] += rmsFileSkyArea

        baseName = os.path.basename(rmsFile).split('.fits')[0]

        rmsPlotFile = os.path.join(rmsDirPath, baseName + "_histplot.png")

        rmsHistFile = os.path.join(rmsDirPath, baseName + "_rmshist.txt")

        countsInBins = np.zeros(len(binCentres))

        if os.path.exists(rmsHistFile):
            rmsHistFromFile = np.loadtxt(rmsHistFile)
            binCentresFile, _, areaSqDegFile, cumulativeAreaSqDegFile = rmsHistFromFile[:, 0], rmsHistFromFile[:, 1], rmsHistFromFile[:, 2], rmsHistFromFile[:, 3]

            if len(binCentres) != len(binCentresFile) or not np.allclose(binCentres, binCentresFile):
                countsInBins, areaSqDegInBins = getRMSAreaCoverage(rmsFile, rmsBins)
            else:
                areaSqDegInBins = areaSqDegFile
                cumulativeAreaSqDegInBins = cumulativeAreaSqDegFile
        else:
            countsInBins, areaSqDegInBins = getRMSAreaCoverage(rmsFile, rmsBins)

        cumulativeNPixInBins = np.cumsum(countsInBins)
        cumulativeAreaSqDegInBins = np.cumsum(areaSqDegInBins)
        np.savetxt(rmsHistFile,
                   np.column_stack((binCentres, countsInBins, areaSqDegInBins, cumulativeNPixInBins, cumulativeAreaSqDegInBins)),
                   fmt='%.6f\t%d\t%.6f\t%d\t%.6f',
                   header='RMS(Jy/beam)\tNPixel\tArea(sq.deg.)\tCumulativeNPixel\tCumulativeArea(sq.deg)')

        globalRMSNPixelsInBinsDict[bandKey] += countsInBins
        globalRMSCumulativeNPixelsInBinsDict[bandKey] += cumulativeNPixInBins
        globalRMSAreaInBinsDict[bandKey] += areaSqDegInBins
        globalRMSCumulativeAreaInBinsDict[bandKey] += cumulativeAreaSqDegInBins

        if not os.path.exists(rmsPlotFile):

            fig,ax=plt.subplots(nrows=1,ncols=1)
            fig.set_size_inches(5,4)

            ax.plot(binCentres, cumulativeAreaSqDegInBins/rmsFileSkyArea)
            ax.axhline(y=1.0, linestyle='dashed')
            ax.set_xscale('log')
            ax.set_xlabel("RMS Noise (Jy/beam)")
            ax.set_ylabel("Fraction of area")
            ax.set_xlim(1E-7, 1E-1)
            plt.savefig(rmsPlotFile, dpi=700, bbox_inches='tight')
            plt.close()

    # Plot cumulative area fraction

    fig,ax=plt.subplots(nrows=1,ncols=3,sharex=True, sharey=False)
    fig.set_size_inches(12,3)

    for bandi, band in enumerate(orderedBands):

        areaCoverageTxtOutName = areaCoveragePlotOutName.replace(".png", "")

        np.savetxt(areaCoverageTxtOutName+'_%s.txt' %band,
                   np.column_stack((binCentres,
                                    globalRMSNPixelsInBinsDict[band],
                                    globalRMSAreaInBinsDict[band],
                                    globalRMSCumulativeNPixelsInBinsDict[band],
                                    globalRMSCumulativeAreaInBinsDict[band])),
                   fmt='%.6f\t%d\t%.6f\t%d\t%.6f',
                   header='RMS(Jy/beam)\tNPixel\tArea(sq.deg.)\tCumulativeNPixel\tCumulativeArea(sq.deg)')

        cumulativeArea = globalRMSCumulativeAreaInBinsDict[band]

        ax[bandi].plot(binCentres, cumulativeArea/totalRMSAreaDict[band], color=bandColorDict[band])

        ax[bandi].axhline(y=1.0, linestyle='dashed', color='k')

        ax[bandi].set_xlabel("RMS Noise (Jy/beam)")
        ax[bandi].set_xlim(1E-7, 1E-1)
        ax[bandi].text(0.5,0.10, "%s band" %(band),transform=ax[bandi].transAxes,ha='center')
        ax[bandi].set_xscale('log')


    ax[0].set_ylabel("Fraction of Cumulative Area")


    plt.savefig(areaCoveragePlotOutName, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nCumulative RMS area fraction plotted!\n")
