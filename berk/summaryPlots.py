"""

This module contains tools for creating the summary plots

"""

import os
import numpy as np
from astropy.io import fits
import astropy.table as atpy
from . import startup, tasks, images, catalogs
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
    if len(orderedBands) > 1:
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

    mask = (SBinWidth > 0) & (skyAreaSterdian > 0)
    sourceCount = np.full_like(NSources, np.nan, dtype=float)
    sourceCount[mask] = SMean[mask]**2.5 * (NSources[mask] / (SBinWidth[mask]*skyAreaSterdian[mask]))

    return sourceCount

def getEffectiveAreaInFluxBinsfromRMS(fluxBinCentres, rmsBinsCentres, cumArea, sigmaDetection=5.0):
    """Compute the effective survey area accessible for each flux bin using the RMS histogram.

    Args:
        fluxBinCentres (:obj:`np.ndarray`): Array of flux-bin centres (Jy).
        rmsBinsCentres (:obj:`np.ndarray`): Array of RMS-bin centres (Jy).
        cumArea (:obj:`np.ndarray`): Cumulative area (e.g., in deg²) for each RMS bin.
        sigmaDetection (:obj:`float`, optional): Detection threshold in units of RMS.
            Default is 5.0.

    Returns:
        :obj:`np.ndarray`: Effective survey area corresponding to each flux bin.
    """

    # Ensure rmsBinsCentres is ascending
    if not np.all(np.diff(rmsBinsCentres) > 0):
        sortIdx = np.argsort(rmsBinsCentres)
        rmsBinsCentres = rmsBinsCentres[sortIdx]
        cumArea = cumArea[sortIdx]

    effArea = np.interp(fluxBinCentres / sigmaDetection, rmsBinsCentres, cumArea)

    return effArea

def computeSourceCount(fluxVals, fluxBins, fluxCompleteness, rmsBinCentreJy=None, cumAreaSqDeg=None, corrRMSCoverage=True, commonAreaSqDeg=None, corrCompleteness=True):
    """Compute source counts normalized by S^2.5 for log-spaced flux bins.

    Args:
        fluxVals (:obj:`np.ndarray`): Array of source flux densities (Jy).
        fluxBins (:obj:`np.ndarray`): Array of flux-bin edges (Jy). Must be monotonic.
        fluxCompleteness (:obj:`np.ndarray`): Array of completeness values corresponding to each source in fluxVals.
        rmsBinCentreJy (:obj:`np.ndarray`, optional): Array of RMS-bin centres (Jy). Required only if corrRMSCoverage=True.
        cumAreaSqDeg (:obj:`np.ndarray`, optional): Cumulative survey area (sq.deg.) corresponding to each RMS bin.
            Required if corrRMSCoverage=True OR if corrRMSCoverage=False and
            no commonAreaSqDeg is supplied.
        corrRMSCoverage (:obj:`bool`, optional): Whether to apply RMS-dependent area correction. Default: True.
        commonAreaSqDeg (:obj:`float`, optional): A constant survey area (sq.deg.) used when corrRMSCoverage=False and
            cumAreaSqDeg is not supplied.

    Returns:
        tuple: (bin centers, row counts, row count errors, source count values, source count errors, effective areas).
    """

    fluxBinCentre = np.sqrt(fluxBins[:-1] * fluxBins[1:])

    if corrCompleteness is True:
        fluxVals = fluxVals[fluxCompleteness > 0.0] # to avoid division by zero
        fluxCompleteness = fluxCompleteness[fluxCompleteness > 0.0]
        weights = 1.0/fluxCompleteness

        fluxCounts, binEdges = np.histogram(fluxVals, bins=fluxBins, weights=weights)
        fluxCountsErr = np.sqrt(np.histogram(fluxVals, bins=fluxBins, weights=weights**2)[0])

    else:
        fluxCounts, binEdges = np.histogram(fluxVals, bins=fluxBins)
        fluxCountsErr = np.sqrt(fluxCounts) # Poisson error

    fluxBinWidths = np.diff(binEdges)

    if corrRMSCoverage is True:
        effAreaInFluxBins = getEffectiveAreaInFluxBinsfromRMS(fluxBinCentre, rmsBinCentreJy, cumAreaSqDeg, sigmaDetection=5.0)
    else:
        if cumAreaSqDeg is not None:
            effAreaInFluxBins = [np.max(cumAreaSqDeg) for i in range(len(fluxBinCentre))]
        elif commonAreaSqDeg is not None:
            effAreaInFluxBins = [commonAreaSqDeg for i in range(len(fluxBinCentre))]
        else:
            raise ValueError("Either cumAreaSqDeg or commonAreaSqDeg must be provided when corrRMSCoverage=False.")

    sourceCountValues = getS2p5dNdS(fluxBinCentre, fluxCounts, fluxBinWidths, effAreaInFluxBins)
    sourceCountErr = getS2p5dNdS(fluxBinCentre, fluxCountsErr, fluxBinWidths, effAreaInFluxBins)

    return fluxBinCentre, fluxCounts, fluxCountsErr, sourceCountValues, sourceCountErr, effAreaInFluxBins

def plotSourceCounts(surveyCatDir, catSubScript, fluxCol='Total_flux', fluxMin=None, fluxMax=None, nFluxBins=50, bandColorDict=None, plotOutPath=None):
    """Plot Euclidean-normalized source counts for each band using survey catalogs.

    Args:
        surveyCatDir (:obj:`str`): Directory containing the survey catalog FITS files for each band.
        catSubScript (:obj:`str`): Subscript for the catalog file names.
        fluxCol (:obj:`str`): Name of the flux column in the catalog table.
        fluxMin (:obj:`float`): Minimum flux to be considered.
        fluxMax (:obj:`float`): Maximum flux to be considered.
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

        catFileName = surveyCatDir+os.path.sep+"survey_catalog_%s%s.fits" %(band, catSubScript)
        catalogTab = atpy.Table().read(catFileName)

        # Converting fluxes to 1.4 GHz
        catalogTab['%s_1p4GHz' %fluxCol] = catalogs.convertFluxFreq(catalogTab[fluxCol].value, catalogTab['freqGHz'].value, 1.4, alpha=0.7)
        fluxVals = catalogTab['%s_1p4GHz' %fluxCol].value
        fluxUnitLabel = catalogTab[fluxCol].unit
        completeness = catalogTab['completeness'] if 'completeness' in catalogTab.colnames else np.ones_like(fluxVals)

        fluxMin = np.min(fluxVals) if fluxMin is None else fluxMin
        fluxMax = np.max(fluxVals) if fluxMax is None else fluxMax
        fluxBins = np.logspace(np.log10(fluxMin), np.log10(fluxMax), nFluxBins+1)

        RMSAreaCoverageName = surveyCatDir+os.path.sep+"MeerKAT_RMS_area_coverage_unique%s_%s.txt" %(catSubScript, band)

        RMSAreaCoverage = np.loadtxt(RMSAreaCoverageName)
        rmsBinCentreJy = RMSAreaCoverage[:, 0]
        cumAreaSqDeg = RMSAreaCoverage[:, 4]

        fluxBinCentre, rowCount, rawCountErr, sourceCountUncorr, sourceCountErrUncorr, effAreaInFluxBins = computeSourceCount(fluxVals, fluxBins, completeness, rmsBinCentreJy, cumAreaSqDeg, corrRMSCoverage=False, corrCompleteness=True)
        fluxBinCentre, rowCount, rawCountErr, sourceCountCorr, sourceCountErrCorr, effAreaInFluxBins = computeSourceCount(fluxVals, fluxBins, completeness, rmsBinCentreJy, cumAreaSqDeg, corrRMSCoverage=True, corrCompleteness=True)

        ax.errorbar(fluxBinCentre, sourceCountCorr, sourceCountErrCorr, mec='k', mew=0.5, mfc=bandColorDict[band], ecolor=bandColorDict[band], marker='o', ms=6, alpha=1, ls='None', label="%s-band" %band)

        textFileName = plotOutPath.split('.png')[0]+'_%s.txt' %band
        np.savetxt(textFileName,
                   np.column_stack((fluxBinCentre,
                                    rowCount,
                                    rawCountErr,
                                    sourceCountUncorr, sourceCountErrUncorr,
                                    sourceCountCorr, sourceCountErrCorr)),
                   fmt='%.6f\t%d\t%d\t%.6f\t%.6f\t%.6f\t%.6f',
                   header='Flux(Jy)\tRawCount\tRawCountErr\tSourceCountUncorr\tSourceCountUncorrErr\tSourceCountCorr\tSourceCountCorrErr')

    ax.set_xlabel("Total Flux (%s)" %fluxUnitLabel)
    ax.set_ylabel(r"$S^{2.5} \ \mathrm{d}N/\mathrm{d}S$ $(\mathrm{Jy}^{1.5} \mathrm{sr}^{-1})$")

    ax.set_xscale('log')
    ax.set_yscale('log')
    if len(orderedBands) > 1:
        plt.legend(loc="lower right")
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

def plotRMSAreaCoverageCumulative(rmsDirPath, areaCoveragePlotOutName, bandColorDict, nRMSBins=100):
    """Plot cumulative sky area as a function of RMS noise for all bands.

    Args:
        areaCoveragePlotOutName (str): Path to save the output plot.
        bandColorDict (dict): Dictionary mapping band names to color codes.
        nRMSBins (int, optional): Number of RMS bins to use (default is 100).

    Returns:
        None. Saves the cumulative RMS area plot to the specified path.
    """

    rmsFiles = sorted(glob.glob(rmsDirPath + os.path.sep + "*rms.fits"))

    rmsBins = np.logspace(-10, 0, nRMSBins)
    binCentres = 0.5 * (rmsBins[1:] + rmsBins[:-1])

    orderedBands = list(bandColorDict.keys())
    numBands = len(orderedBands)
    globalRMSNPixelsInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSCumulativeNPixelsInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSAreaInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSCumulativeAreaInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    totalRMSAreaDict = {band: 0 for band in orderedBands}

    print("\nCalculating RMS area coverage ... \n")

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
            binCentresFile, _, areaSqDegFile, _, cumulativeAreaSqDegFile = rmsHistFromFile[:, 0], rmsHistFromFile[:, 1], rmsHistFromFile[:, 2], rmsHistFromFile[:, 3], rmsHistFromFile[:, 4]

            if len(binCentres) != len(binCentresFile) or not np.allclose(np.round(binCentres,6), np.round(binCentresFile,6)):
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

            ax.plot(binCentres, cumulativeAreaSqDegInBins)
            ax.axhline(y=rmsFileSkyArea, linestyle='dashed')
            ax.set_xscale('log')
            ax.set_xlabel("RMS Noise (Jy/beam)")
            ax.set_ylabel("Sky Area (sq. deg.)")
            ax.set_xlim(1E-7, 1E-1)
            plt.savefig(rmsPlotFile, dpi=700, bbox_inches='tight')
            plt.close()

    # Plot cumulative area

    figWidthPerPanel = 4
    figHeight = 3
    fig,ax=plt.subplots(nrows=1,ncols=numBands,sharex=True, sharey=False, squeeze=False)
    fig.set_size_inches(figWidthPerPanel * numBands, figHeight)
    ax = ax.flatten()

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
        totalArea = cumulativeArea[-1]

        ax[bandi].plot(binCentres, cumulativeArea, color=bandColorDict[band])

        #ax[bandi].axhline(y=totalRMSAreaDict[band], linestyle='dashed', color='k')
        ax[bandi].axhline(y=totalArea, linestyle='dashed', color='k')

        ax[bandi].set_xlabel("RMS Noise (Jy/beam)")
        ax[bandi].set_xlim(1E-7, 1E-1)
        ax[bandi].text(0.75,0.15, "%s band\n(%0.2f sq. deg.)" %(band, totalRMSAreaDict[band]),transform=ax[bandi].transAxes,ha='center')
        ax[bandi].set_xscale('log')

    ax[0].set_ylabel("Sky Area (sq. deg.)")

    plt.savefig(areaCoveragePlotOutName, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nCumulative RMS area plotted!\n")

def plotRMSAreaCoverageCumulativeUnique(rmsDirPath, areaCoveragePlotOutName, bandColorDict, imagesTab, overlapGroupsByBand, isolatedIndicesByBand, nRMSBins=100):
    """Compute and plot the cumulative unique sky area as a function of RMS
    noise for different observing bands, accounting for overlapping pointings.

    Overlapping pointings within each band are combined into RMS mosaics before
    calculating the area coverage, while isolated pointings are treated
    independently. The resulting RMS-dependent area coverage is computed
    separately for each band and plotted.

    Args:
        rmsDirPath (:obj:`str`): Directory containing RMS FITS images for individual pointings.
        areaCoveragePlotOutName (:obj:`str`): Path to save the output area coverage plot.
        bandColorDict (:obj:`dict`): Dictionary mapping band names to colours used for plotting.
        imagesTab (:obj:`astropy.table.Table`): Table containing pointing information, including RMS image/catalogue paths.
        overlapGroupsByBand (:obj:`dict`): Dictionary containing lists of overlapping pointing groups for each band. Each group contains
            indices corresponding to rows in ``imagesTab``.
        isolatedIndicesByBand (:obj:`dict`): Dictionary containing lists of isolated pointing indices for each band. Indices correspond to rows
            in ``imagesTab``.
        nRMSBins (:obj:`int`, optional): Number of RMS noise bins used for calculating the area coverage. Defaults to 100.

    Returns:
        None: Saves the cumulative RMS area coverage plot to ``areaCoveragePlotOutName``.
    """
    rmsBins = np.logspace(-10, 0, nRMSBins)
    binCentres = 0.5 * (rmsBins[1:] + rmsBins[:-1])

    orderedBands = list(bandColorDict.keys())
    numBands = len(orderedBands)
    globalRMSNPixelsInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSCumulativeNPixelsInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSAreaInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}
    globalRMSCumulativeAreaInBinsDict = {band: np.zeros(len(binCentres)) for band in orderedBands}

    # Process isolated pointings band-by-band
    for band, isolatedIndices in isolatedIndicesByBand.items():

        print("Processing %d isolated %s-band pointings..." %
              (len(isolatedIndices), band))

        for i in isolatedIndices:

            rmsFile = os.path.join(
                rmsDirPath,
                os.path.basename(imagesTab['radioCatPath'][i])
                .replace('_srl_bdsfcat.fits', '_rms.fits')
            )

            countsInBins, areaSqDegInBins = getRMSAreaCoverage(rmsFile=rmsFile, rmsBins=rmsBins)

            cumulativeNPixInBins = np.cumsum(countsInBins)
            cumulativeAreaSqDegInBins = np.cumsum(areaSqDegInBins)

            globalRMSNPixelsInBinsDict[band] += countsInBins
            globalRMSCumulativeNPixelsInBinsDict[band] += cumulativeNPixInBins
            globalRMSAreaInBinsDict[band] += areaSqDegInBins
            globalRMSCumulativeAreaInBinsDict[band] += cumulativeAreaSqDegInBins

    # Process overlapping groups band-by-band
    for band, overlapGroups in overlapGroupsByBand.items():
        for gi, group in enumerate(overlapGroups):

            rmsFileList = [os.path.join(rmsDirPath, os.path.basename(imagesTab['radioCatPath'][i]).replace('_srl_bdsfcat.fits', '_rms.fits')) for i in group]

            mosaicFilesDir = os.path.join(rmsDirPath, "mosaics")
            os.makedirs(mosaicFilesDir, exist_ok=True)
            mosaicFile = os.path.join(mosaicFilesDir, "%s_mosaic_group%d.fits" %(band, gi))

            images.mosaicRMSMaps(rmsFileList, outputFile=mosaicFile, resolutionArcmin=1.0)

            countsInBins, areaSqDegInBins = getRMSAreaCoverage(rmsFile=mosaicFile,rmsBins=rmsBins)

            cumulativeNPixInBins = np.cumsum(countsInBins)
            cumulativeAreaSqDegInBins = np.cumsum(areaSqDegInBins)

            globalRMSNPixelsInBinsDict[band] += countsInBins
            globalRMSCumulativeNPixelsInBinsDict[band] += cumulativeNPixInBins
            globalRMSAreaInBinsDict[band] += areaSqDegInBins
            globalRMSCumulativeAreaInBinsDict[band] += cumulativeAreaSqDegInBins

    # Plot cumulative area

    figWidthPerPanel = 4
    figHeight = 3
    fig,ax=plt.subplots(nrows=1,ncols=numBands,sharex=True, sharey=False, squeeze=False)
    fig.set_size_inches(figWidthPerPanel * numBands, figHeight)
    ax = ax.flatten()

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
        totalArea = cumulativeArea[-1]

        ax[bandi].plot(binCentres, cumulativeArea, color=bandColorDict[band])

        ax[bandi].axhline(y=totalArea, linestyle='dashed', color='k')

        ax[bandi].set_xlabel("RMS Noise (Jy/beam)")
        ax[bandi].set_xlim(1E-7, 1E-1)
        ax[bandi].text(0.75,0.15, "%s band\n(%0.2f sq. deg.)" %(band, totalArea),transform=ax[bandi].transAxes,ha='center')
        ax[bandi].set_xscale('log')

    ax[0].set_ylabel("Sky Area (sq. deg.)")

    plt.savefig(areaCoveragePlotOutName, dpi=700, bbox_inches='tight')
    plt.close()
    print("\nCumulative RMS area plotted for unique coverage!\n")