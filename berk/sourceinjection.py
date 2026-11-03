"""

Routines for estimating the point-source completenss.

"""

import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import copy
import warnings
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.modeling.models import Gaussian2D
from multiprocessing import Pool
from berk import catalogs, crossmatch, startup
import bdsf

# ---------------------------------------------------------------------------
# PyBDSF settings - we might need to take it from the cache dir later on
# ---------------------------------------------------------------------------

pybdsfArgs = {
    "process_image": {
        "thresh_isl":        3.0,
        "thresh_pix":        5.0,
        "group_by_isl":      True,
        "rms_box":           (150, 15),
        "adaptive_rms_box":  True,
        "adaptive_thresh":   100.0,
        "rms_box_bright":    (40, 15),
        "atrous_do":         True,
        "atrous_orig_isl":   True,
        "atrous_jmax":       3,
    },
    "export_image": {
        "ch0":           False,
        "rms":           False,
        "mean":          False,
        "pi":            False,
        "gaus_resid":    False,
        "shap_resid":    False,
        "psf_major":     False,
        "psf_minor":     False,
        "psf_pa":        False,
        "psf_ratio":     False,
        "psf_ratio_aper":False,
        "island_mask":   False,
    },
}

def sourceExtract(imageFile, pybdsfBaseName, saveRMSMeanResidualMaps=False,
                  outDir=None, rmsFileName=None, meanFileName=None):
    """
    Extracts sources from the given image using PyBDSF and saves the catalogue and optionally the RMS, mean and residual maps.
     Args:
        imageFile: Path to the input image file (FITS format).
        pybdsfBaseName: Base name for the output files (catalogue and maps).
        saveRMSMeanResidualMaps: If True, saves the RMS, mean and residual maps as well.
        outDir: Directory to save the output files. Default is current working directory.
        rmsFileName: Optional path to the RMS map file to be used by PyBDSF. If None, it will not be used.
        meanFileName: Optional path to the mean map file to be used by PyBDSF. If None, it will not be used.
    """

    if outDir is None:
        outDir = os.getcwd()

    outCatalog = os.path.join(outDir, pybdsfBaseName+'_srl_bdsfcat.fits')

    if os.path.exists(outCatalog):
        print("\nPyBDSF catalogue already exists. Skipping %s image\n" %imageFile)
        return

    print("\nExtracting sources from %s\n" %imageFile)

    pybdsfArgsCopy = copy.deepcopy(pybdsfArgs)

    if (rmsFileName is not None) and (meanFileName is not None):
        pybdsfArgsCopy['process_image']['advanced_opts'] = True
        pybdsfArgsCopy['process_image']['rmsmean_map_filename'] = [
            meanFileName,
            rmsFileName
        ]   # have to check this - but looks like this advance settings
            # do not take the absolute path to rms and mean files.
            # We need to have a copy in the working directory.

    if saveRMSMeanResidualMaps is True:
        for mapType in ['rms', 'mean','gaus_resid']:
            pybdsfArgsCopy['export_image'][mapType] = True

    try:
        img = bdsf.process_image(imageFile, **pybdsfArgsCopy['process_image'])
        for imgType in pybdsfArgsCopy['export_image']:
            if pybdsfArgsCopy['export_image'][imgType] is True:
                img.export_image(outfile=os.path.join(outDir, pybdsfBaseName+'_'+imgType+'.fits'), clobber=True, img_type=imgType)

        img.write_catalog(outfile=outCatalog, format='fits',
                          catalog_type='srl', clobber=True)

    except Exception as e:
        print("PyBDSF FAILED on %s" % imageFile)
        print("Error: %s" % e)
        print("Skipping this file and moving on...")
        return None

def makeExtendedGaussian(localX, localY, x0, y0, sourceMajDeg, sourceMinDeg, sourcePADeg, beamMajDeg, beamMinDeg, beamPADeg, wcsObj, fluxJy):
    """
    Creates a beam-convolved elliptical Gaussian source in Jy/beam.

    Args:
        - localX: A 2D numpy array containing the x pixel coordinates of the source stamp.
        - localY: A 2D numpy array containing the y pixel coordinates of the source stamp.
        - x0: The x pixel coordinate of the source centre.
        - y0: The y pixel coordinate of the source centre.
        - sourceMajDeg: The intrinsic major-axis FWHM of the source in degrees.
        - sourceMinDeg: The intrinsic minor-axis FWHM of the source in degrees.
        - sourcePADeg: The intrinsic source position angle in degrees east of north.
        - beamMajDeg: The restoring beam major-axis FWHM in degrees.
        - beamMinDeg: The restoring beam minor-axis FWHM in degrees.
        - beamPADeg: The restoring beam position angle in degrees east of north.
        - wcsObj: A two-dimensional Astropy WCS object describing the celestial coordinate system.
        - fluxJy: The integrated flux density of the source in Jy.

    Returns:
        - A 2D numpy array containing the beam-convolved source profile in Jy/beam.
    """

    fwhmToSigma = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    def skyCovariance(majorDeg, minorDeg, paDeg):
        major = majorDeg * fwhmToSigma
        minor = minorDeg * fwhmToSigma
        pa = np.deg2rad(paDeg)

        # Coordinates are east and north on the sky.
        majorVector = np.array([np.sin(pa), np.cos(pa)])
        minorVector = np.array([np.cos(pa), -np.sin(pa)])

        return (
            major**2 * np.outer(majorVector, majorVector)
            + minor**2 * np.outer(minorVector, minorVector)
        )

    sourceCovSky = skyCovariance(
        sourceMajDeg, sourceMinDeg, sourcePADeg
    )
    beamCovSky = skyCovariance(
        beamMajDeg, beamMinDeg, beamPADeg
    )

    # Local Jacobian mapping pixel offsets (dx, dy) to
    # tangent-plane sky offsets (east, north), in degrees.
    eps = 1e-3  # pixel

    ra0, dec0 = wcsObj.all_pix2world(x0, y0, 0)

    raX, decX = wcsObj.all_pix2world(x0 + eps, y0, 0)
    raY, decY = wcsObj.all_pix2world(x0, y0 + eps, 0)

    # Wrap RA differences and project them onto the local east direction.
    draX = ((raX - ra0 + 180.0) % 360.0 - 180.0) * np.cos(np.deg2rad(dec0))
    draY = ((raY - ra0 + 180.0) % 360.0 - 180.0) * np.cos(np.deg2rad(dec0))

    jacobian = np.array([
        [draX / eps, draY / eps],
        [(decX - dec0) / eps, (decY - dec0) / eps],
    ])

    if not np.all(np.isfinite(jacobian)) or abs(np.linalg.det(jacobian)) < 1e-20:
        raise ValueError("Invalid local WCS Jacobian")

    invJacobian = np.linalg.inv(jacobian)

    sourceCovPix = invJacobian @ sourceCovSky @ invJacobian.T
    beamCovPix = invJacobian @ beamCovSky @ invJacobian.T

    observedCovPix = sourceCovPix + beamCovPix

    dx = localX - x0
    dy = localY - y0

    invCov = np.linalg.inv(observedCovPix)

    exponent = (
        invCov[0, 0] * dx**2
        + 2.0 * invCov[0, 1] * dx * dy
        + invCov[1, 1] * dy**2
    )

    # Normalise so the integrated Gaussian flux is fluxJy.
    peakJyBeam = fluxJy * np.sqrt(
        np.linalg.det(beamCovPix) / np.linalg.det(observedCovPix)
    )

    return peakJyBeam * np.exp(-0.5 * exponent)

ARCSEC_TO_DEG = 1.0 / 3600.0

def drawRandomSourceSize(rng, sizeCfg):
    """
    Draws the intrinsic (deconvolved) size of one injected source.

    Args:
        rng: numpy random Generator.
        sizeCfg: None for a point source, otherwise a dict with the keys
            - minSizeArcsec, maxSizeArcsec: range of the deconvolved major-axis FWHM (arcsec).
            - spacing: 'log' (default, uniform in log10) or 'linear'.
            - minAxisRatio: the minor/major axis ratio is drawn uniformly in [minAxisRatio, 1].
            The position angle is drawn uniformly in [0, 180) degrees.

    Returns:
        (major axis FWHM [deg], minor axis FWHM [deg], position angle [deg])
    """
    if sizeCfg is None:
        return 0.0, 0.0, 0.0

    lo, hi = sizeCfg['minSizeArcsec'], sizeCfg['maxSizeArcsec']
    if sizeCfg.get('spacing', 'log') == 'log':
        majArcsec = 10 ** rng.uniform(np.log10(lo), np.log10(hi))
    else:
        majArcsec = rng.uniform(lo, hi)
    axisRatio = rng.uniform(sizeCfg.get('minAxisRatio', 0.3), 1.0)
    paDeg = rng.uniform(0.0, 180.0)

    return majArcsec * ARCSEC_TO_DEG, majArcsec * axisRatio * ARCSEC_TO_DEG, paDeg

def _sizeBinEdges(sizeCfg, nSizeBins):
    """Bin edges (arcsec) for the size axis of the 2D completeness grid."""
    lo, hi = sizeCfg['minSizeArcsec'], sizeCfg['maxSizeArcsec']
    if sizeCfg.get('spacing', 'log') == 'log':
        return np.logspace(np.log10(lo), np.log10(hi), nSizeBins + 1)
    return np.linspace(lo, hi, nSizeBins + 1)

def _runSingleInjectionToResidual(args):
    """
    Runs a single injection of fake sources to the residual image and saves the injected image and the corresponding fake source catalogue.

    Args:
     args: A tuple containing the following elements:
        - rep: The repetition number (used for random seed).
        - img: The residual image data as a 2D numpy array.
        - hdr: The FITS header of the residual image.
        - beam: A Gaussian2D model representing the beam of the image.
        - nInjectionSources: The number of fake sources to inject.
        - fieldRACentre: The right ascension of the field centre in degrees.
        - fieldDecCentre: The declination of the field centre in degrees.
        - injectRadius: The radius within which to inject sources in degrees.
        - minFluxJyInj: The minimum flux density of injected sources in Jy.
        - maxFluxJyInj: The maximum flux density of injected sources in Jy.
        - sinjectDir: The directory to save the injected images and catalogues.
        - imageBaseName: The base name for the injected image and catalogue files.
        - rCutOff: The cutoff radius in pixels for the PSF stamp around each injected source.

    """

    (rep, img, hdr, beam, nInjectionSources,
     fieldRACentre, fieldDecCentre,
     injectRadius, minFluxJyInj, maxFluxJyInj,
     sinjectDir, imageBaseName, rCutOff, sizeCfg, crowdingFactor) = args

    rng = np.random.default_rng(seed=rep)

    try:

        sInjectOutFile = os.path.join(
            sinjectDir, "%s_sinjected_%dsources_%d_image.fits" %(imageBaseName,
                                                                nInjectionSources,
                                                                rep)
        )

        if os.path.exists(sInjectOutFile):
            print("\n%s file already exists. Skipping this injection\n" %sInjectOutFile)
            return

        wcsObj = WCS(hdr, naxis=2)
        updatedHeader = hdr.copy()

        print("\nInjecting %d sources on repetition %d ...\n" %(nInjectionSources, rep+1))
        copyImg = img.copy()

        ny, nx = copyImg.shape

        logMin = np.log10(minFluxJyInj)
        logMax = np.log10(maxFluxJyInj)
        beamMajDeg = max(hdr["BMAJ"], hdr["BMIN"])

        # to keep track of previously injected sources to avoid crowding
        injectedCoords = []
        injectedObsFwhmDeg = []   # observed (beam-convolved) FWHM of each injected source

        nSourcesInjected = 0

        randomRAInjected = []
        randomDecInjected = []
        randomFluxInjected= []
        randomSourceMajDegInjected = []
        randomSourceMinDegInjected = []
        randomSourcePADegInjected = []

        maxAttempts = nInjectionSources * 10  # maximum attempts to inject sources
        nAttempts = 0

        while nSourcesInjected < nInjectionSources:

            nAttempts += 1
            if nAttempts > maxAttempts:
                print("Warning: could only inject %d / %d sources before hitting max attempts" %(nSourcesInjected, nInjectionSources))
                break

            randomRA, randomDec = crossmatch.randomPointsInCircleExactNPoints(fieldRACentre, fieldDecCentre, injectRadius, 1, rng=rng)

            # Random intrinsic (deconvolved) size of this source, drawn from the configured range
            sourceMajDeg, sourceMinDeg, sourcePADeg = drawRandomSourceSize(rng, sizeCfg)
            obsFwhmDeg = float(np.hypot(beamMajDeg, sourceMajDeg))   # observed major-axis FWHM

            # Minimum separation to the sources injected so far. It scales with the observed
            # size of both sources: crowdingFactor * (FWHM_i + FWHM_j) / 2, so that two point
            # sources are separated by crowdingFactor beams. (Flat-sky approximation.)
            if len(injectedCoords) > 0:
                prev = np.asarray(injectedCoords)
                dRA = (prev[:, 0] - randomRA[0] + 180.0) % 360.0 - 180.0
                dRA = dRA * np.cos(np.deg2rad(randomDec[0]))
                dDec = prev[:, 1] - randomDec[0]
                seps = np.hypot(dRA, dDec)
                minSeps = crowdingFactor * 0.5 * (obsFwhmDeg + np.asarray(injectedObsFwhmDeg))
                if np.any(seps < minSeps):
                    continue

            xPix, yPix = wcsObj.wcs_world2pix(randomRA, randomDec, 0)

            x0 = xPix[0]
            y0 = yPix[0]

            xMin = max(0, int(np.floor(x0 - rCutOff)))
            xMax = min(nx, int(np.ceil(x0 + rCutOff)))
            yMin = max(0, int(np.floor(y0 - rCutOff)))
            yMax = min(ny, int(np.ceil(y0 + rCutOff)))

            if xMax <= xMin or yMax <= yMin:
                continue # source fell outside the image

            injectedCoords.append([randomRA[0], randomDec[0]])
            injectedObsFwhmDeg.append(obsFwhmDeg)

            randomFlux = 10 ** rng.uniform(logMin, logMax)

            localY, localX = np.mgrid[yMin:yMax, xMin:xMax]

            if sourceMajDeg < 1e-10:
                psfStamp = beam(localX - x0, localY - y0)
                copyImg[yMin:yMax, xMin:xMax] += randomFlux * psfStamp
            else:
                extendedStamp = makeExtendedGaussian(localX, localY, x0, y0,sourceMajDeg, sourceMinDeg, sourcePADeg,hdr["BMAJ"], hdr["BMIN"], hdr["BPA"], wcsObj, randomFlux)
                copyImg[yMin:yMax, xMin:xMax] += extendedStamp

            nSourcesInjected += 1
            randomRAInjected.append(float(randomRA[0]))
            randomDecInjected.append(float(randomDec[0]))
            randomFluxInjected.append(float(randomFlux))
            randomSourceMajDegInjected.append(float(sourceMajDeg))
            randomSourceMinDegInjected.append(float(sourceMinDeg))
            randomSourcePADegInjected.append(float(sourcePADeg))

        # Save the injected image

        updatedHeader.update(wcsObj.to_header())
        updatedHeader['BUNIT']='JY/BEAM'
        hdu = fits.PrimaryHDU(data=copyImg, header=updatedHeader)
        hdu.writeto(sInjectOutFile, overwrite=True)

        # Save the fake source catalogue

        fakeSourcesTab = Table({
            "RADeg_injected":  randomRAInjected,
            "decDeg_injected": randomDecInjected,
            "fluxJy_injected": randomFluxInjected,
            "sourceMajDeg_injected": randomSourceMajDegInjected,
            "sourceMinDeg_injected": randomSourceMinDegInjected,
            "sourcePADeg_injected": randomSourcePADegInjected
        })
        fakeSourcesOutFile = os.path.join(sinjectDir, "%s_sinjected_%dsources_%d_fakesources.fits" %(imageBaseName, nInjectionSources, rep))
        fakeSourcesTab.write(fakeSourcesOutFile, overwrite=True)

        print("%d sources injected to %s\n" %(nSourcesInjected, sInjectOutFile))

    except Exception as e:
        print("Error during injection of %s: %s" % (sInjectOutFile, e))
        print("Skipping this injection and moving on...")
        return None

def injectImage(imageToInjectFileName, fieldRACentre, fieldDecCentre, injectRadius,
                sinjectDir, imageBaseName, nInjectionSources=100,
                minFluxJyInj=1E-5, maxFluxJyInj=1.0, nRepetitions=100, sizeCfg=None, crowdingFactor=5.0):
    """
    Injects fake sources into the residual image and saves the injected images and the corresponding fake source catalogues.

    Args:
        imageToInjectFileName: Name of the residual image file to inject sources into.
        fieldRACentre: Right ascension of the field centre (deg).
        fieldDecCentre: Declination of the field centre (deg).
        injectRadius: Radius within which to inject sources (deg).
        sinjectDir: Directory to save the injected images and catalogues.
        imageBaseName: Base name for the injected image and catalogue files.
        nInjectionSources: Number of fake sources to inject.
        minFluxJyInj: Minimum flux of injected sources (Jy).
        maxFluxJyInj: Maximum flux of injected sources (Jy).
        nRepetitions: Number of repetitions for each injection.

    Returns:
        None
    """

    print("\nInjecting to residual image %s\n" %imageToInjectFileName)

    try:
        hdul = fits.open(os.path.join(sinjectDir, imageToInjectFileName))
    except Exception as e:
        print("Error: %s" % e)
        print("Skipping this file and moving on...")
        return None

    img = hdul[0].data.squeeze()
    hdr = hdul[0].header

    bmajDeg = hdr["BMAJ"]     # Beam major axis (deg)
    bminDeg = hdr["BMIN"]     # Beam minor axis (deg)
    bpaDeg  = hdr["BPA"]      # Beam position angle (deg)
    pixelScale = min(abs(hdr["CDELT1"]), abs(hdr["CDELT2"]))
    if sizeCfg is not None:
        maxSourceMajDeg = sizeCfg['maxSizeArcsec'] * ARCSEC_TO_DEG
    else:
        maxSourceMajDeg = 0.0
    maxObservedFwhmDeg = np.sqrt(bmajDeg**2 + maxSourceMajDeg**2)
    sigmaMaxPix = (maxObservedFwhmDeg / pixelScale / (2.0 * np.sqrt(2.0 * np.log(2.0))))
    rCutOff = int(np.ceil(7.0 * sigmaMaxPix))

    cDelt = np.abs(hdr["CDELT1"]) # pixel scale (deg/pix)

    # Converting beam FWHM from deg to pixels
    bmajPix = bmajDeg / cDelt
    bminPix = bminDeg / cDelt

    # Converting beam FWHM to Gaussian sigma
    sigmaX = bmajPix / (2.0 * np.sqrt(2 * np.log(2.0)))
    sigmaY = bminPix / (2.0 * np.sqrt(2 * np.log(2.0)))

    # Gaussian beam definition
    beam = Gaussian2D(1.0, 0, 0, sigmaX, sigmaY, theta=np.deg2rad(bpaDeg))

    # Stamp cutoff: 5 * FWHM in pixels
    #rCutOff = int(5 * max(sigmaX, sigmaY)) + 1

    argsList = [
        (rep, img, hdr, beam, nInjectionSources, fieldRACentre, fieldDecCentre, injectRadius,
         minFluxJyInj, maxFluxJyInj, sinjectDir, imageBaseName, rCutOff, sizeCfg, crowdingFactor)
        for rep in range(nRepetitions)
    ]

    nProcesses = int(os.cpu_count())
    with Pool(processes=nProcesses) as pool:
        pool.map(_runSingleInjectionToResidual, argsList)

    print("\nAll injections done. Results in %s\n" %sinjectDir)

def matchFakeToRecovered(fakeTab, recovTab, matchRadDeg, fluxTolerance=0.5):
    """
    Matches the injected fake sources to the recovered sources based on their sky positions.
    A fake source is considered recovered if there is a recovered source within matchRadDeg degrees.

     Args:
        fakeTab: Table containing the injected fake sources with columns "RADeg_injected" and "decDeg_injected".
        recovTab: Table containing the recovered sources with columns "RA" and "DEC".
        matchRadDeg: Matching radius in degrees to consider a fake source as recovered.
        fluxTolerance: Fractional tolerance for flux matching. Default is 0.5

     Returns:
        injectedRecoveredMask: Boolean array indicating which injected fake sources were recovered.
        matchedFakeTab: Table of the injected fake sources that were matched to recovered sources.
    """

    #TODO: fluxTolerance criteria is something to play around.

    catFake  = SkyCoord(ra=fakeTab["RADeg_injected"].value * u.deg,
                      dec=fakeTab["decDeg_injected"].value * u.deg)
    catRecov = SkyCoord(ra=recovTab["RA"].value * u.deg,
                      dec=recovTab["DEC"].value * u.deg)

    # For each fake source, find nearest recovered source
    idx, sep2d, _ = catFake.match_to_catalog_sky(catRecov)

    # Position match
    matchRadiusAngle = np.asarray(matchRadDeg, dtype=float) * u.deg   # scalar or one value per fake source
    withinRadius = sep2d < matchRadiusAngle

    # Flux match
    if fluxTolerance is not None:
        recoveredFlux = np.array(recovTab["Total_flux"][idx])
        injectedFlux = np.array(fakeTab["fluxJy_injected"])
        fractionalDiff = np.abs(recoveredFlux - injectedFlux) / injectedFlux
        fluxMatch = fractionalDiff < fluxTolerance
    else:
        fluxMatch = np.ones(len(fakeTab), dtype=bool)

    injectedRecoveredMask = withinRadius & fluxMatch

    matchedFakeTab  = fakeTab[injectedRecoveredMask]

    return injectedRecoveredMask, matchedFakeTab

def getBeamParamsFromFits(fitsFile):
    """
    Extracts the beam parameters (BMAJ, BMIN) from the FITS file header. It first looks for these parameters in the second HDU's COMMENT fields, and if not found, it looks in the primary HDU's header.
     Args:
        fitsFile: Path to the FITS file.

     Returns:
        tuple: A tuple containing the beam major axis (BMAJ) and minor axis (BMIN) in degrees.
    """

    with fits.open(fitsFile) as hdul:
        if len(hdul) > 1:
            hdr = hdul[1].header
            comments = list(hdr.get("COMMENT", []))
            bmaj = bmin = None
            for line in comments:
                m_bmaj = re.search(r"BMAJ\s*=\s*([0-9.eE+\-]+)", str(line))
                m_bmin = re.search(r"BMIN\s*=\s*([0-9.eE+\-]+)", str(line))
                if m_bmaj:
                    bmaj = float(m_bmaj.group(1))
                if m_bmin:
                    bmin = float(m_bmin.group(1))
            if bmaj is not None and bmin is not None:
                return bmaj, bmin

        hdr0 = hdul[0].header
        bmaj = hdr0.get("BMAJ", None)
        bmin = hdr0.get("BMIN", None)

    if bmaj is None or bmin is None:
        raise ValueError("Cannot find BMAJ / BMIN in %s" % fitsFile)

    return bmaj, bmin

def calculateCompleteness(sInjectedCatList, pybdsfCat, imageName, sinjectDir,
                          minFluxJyInj, maxFluxJyInj, nJyBins=30, sizeCfg=None,
                          nFluxBins2D=40, nSizeBins=15, minCountPerCell=10,
                          nSizeGroups=4, sizePlotFluxRangesJy=None,
                          matchRadiusScalesWithSize=True):
    """
    Calculates the completeness (fraction of injected sources recovered) as a function of flux density
    and, if sources with random sizes were injected, also as a function of source size.

    Note that when sizes are random the 1D completeness vs flux density is averaged over the injected
    size distribution, so it depends on the chosen size range.

    Args:
        sInjectedCatList: List of file paths to the PyBDSF catalogues of the injected images.
        pybdsfCat: The PyBDSF catalogue of the original image (without injections).
        imageName: Name of the image being processed (used for plot titles).
        sinjectDir: Directory where the injection results are stored (used for saving plots).
        minFluxJyInj: Minimum flux density of injected sources in Jy.
        maxFluxJyInj: Maximum flux density of injected sources in Jy.
        nJyBins: Number of flux bins to use for the 1D completeness curve.
        sizeCfg: Size configuration used for the injection (None for point sources). If given, the 2D
            (flux density x size) completeness is also computed and plotted.
        nFluxBins2D, nSizeBins: Number of flux and size bins of the 2D grid.
        minCountPerCell: Cells of the 2D grid with fewer injected sources (summed over all repetitions) are masked.
        nSizeGroups: Number of size groups for the "completeness vs flux density in size bins" plot.
        sizePlotFluxRangesJy: List of (min, max) flux ranges in Jy for the "completeness vs size" plot.
            Default: four octaves starting at half the flux density where the size-averaged completeness reaches 50%.
        matchRadiusScalesWithSize: If True, the radius for matching an injected to a recovered source is
            max(beam major axis, 0.5 * observed FWHM of the injected source), so that large sources whose fitted
            centroid is offset from the injected position are not counted as missed. Unchanged for point sources.

    Returns:
        meanFraction: Array of mean completeness values for each flux bin.
        stdFraction: Array of standard deviation of completeness for each flux bin.
        allBinCentres: Array of flux bin centre values.
    """

    allFractionRecovered = []
    allBinCentres        = None

    fluxBins    = np.logspace(np.log10(minFluxJyInj), np.log10(maxFluxJyInj), nJyBins + 1)
    binCentres  = np.sqrt(fluxBins[:-1] * fluxBins[1:])
    allBinCentres = binCentres

    # 2D (flux density x size) bookkeeping, one array per injection run
    allInj2D, allRec2D = [], []
    fluxEdges2D = sizeEdges = None
    if sizeCfg is not None:
        fluxEdges2D = np.logspace(np.log10(minFluxJyInj), np.log10(maxFluxJyInj), nFluxBins2D + 1)
        sizeEdges = _sizeBinEdges(sizeCfg, nSizeBins)
    beamMajArcsec = None

    def accumulate2D(fakeTab, recoveredMask):
        if sizeCfg is None:
            return
        fluxInj = np.asarray(fakeTab["fluxJy_injected"], dtype=float)
        sizeInj = np.asarray(fakeTab["sourceMajDeg_injected"], dtype=float) * 3600.0   # arcsec
        nInj, _, _ = np.histogram2d(fluxInj, sizeInj, bins=[fluxEdges2D, sizeEdges])
        nRec, _, _ = np.histogram2d(fluxInj[recoveredMask], sizeInj[recoveredMask],
                                    bins=[fluxEdges2D, sizeEdges])
        allInj2D.append(nInj)
        allRec2D.append(nRec)

    for catIdx, sInjectedCatFile in enumerate(sInjectedCatList):

        print("\nProcessing %s (%d/%d)" % (sInjectedCatFile, catIdx+1, len(sInjectedCatList)))

        fakeSourceFile = sInjectedCatFile.replace("srl_bdsfcat", "fakesources")
        if not os.path.exists(fakeSourceFile):
            print("Fake-source file not found - skipping.")
            continue

        recovTab    = Table.read(sInjectedCatFile)
        fakeTab     = Table.read(fakeSourceFile)

        tabWrapAngle = catalogs.detectWrapAngle(recovTab['RA'])
        recovTab = catalogs.fixRA(recovTab, raCol='RA', wrapAngle=tabWrapAngle)

        if len(recovTab) == 0:
            print("No sources recovered by PyBDSF - skipping.")
            allFractionRecovered.append(np.zeros(nJyBins))
            accumulate2D(fakeTab, np.zeros(len(fakeTab), dtype=bool))
            continue

        try:
            beamMajor, beamMinor = getBeamParamsFromFits(sInjectedCatFile)
        except ValueError as e:
            print("Error: %s - skipping." % e)
            continue
        beamMajArcsec = beamMajor * 3600.0

        # in case we want to inject sources to an inner area of the image #TODO
        matchRadDeg = max(beamMajor, beamMinor)
        if matchRadiusScalesWithSize:
            obsFwhmDeg = np.hypot(beamMajor, np.asarray(fakeTab["sourceMajDeg_injected"], dtype=float))
            matchRadDeg = np.maximum(matchRadDeg, 0.5 * obsFwhmDeg)   # one radius per fake source

        # injectedRecoveredMask : boolean mask for the injected source table that were recovered by PyBDSF
        # matchedFakeTab : subset of the injected source table that were recovered by PyBDSF
        injectedRecoveredMask, matchedFakeTab = matchFakeToRecovered(fakeTab, recovTab, matchRadDeg, fluxTolerance=None)

        nRecovered = injectedRecoveredMask.sum()
        print("%d / %d injected fake sources recovered" % (nRecovered, len(fakeTab)))

        accumulate2D(fakeTab, injectedRecoveredMask)

        # Sky plot: injected but NOT recovered
        brightMask = fakeTab["fluxJy_injected"] > 1e-2   # Jy
        brightNotRecovered = fakeTab[brightMask & ~injectedRecoveredMask]

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(pybdsfCat["RA"], pybdsfCat["DEC"],
                   s=1, color="gray", alpha=0.3, label="Real sources")
        ax.scatter(fakeTab["RADeg_injected"], fakeTab["decDeg_injected"],
                   s=2, color="blue", alpha=0.5, label="Injected")
        ax.scatter(brightNotRecovered["RADeg_injected"],
                   brightNotRecovered["decDeg_injected"],
                   s=15, color="red", label="Not recovered (>10 mJy)")
        ax.invert_xaxis()
        ax.set_xlabel("RA (deg)")
        ax.set_ylabel("Dec (deg)")
        ax.set_title("%s\nNot recovered (flux > 10 mJy) — run %d" % (os.path.basename(imageName), catIdx))
        ax.legend(markerscale=3)
        plotFile = os.path.join(sinjectDir, sinjectDir, "notRecovered_bright_%03d.png" % catIdx)
        fig.savefig(plotFile, dpi=150, bbox_inches="tight")
        plt.close(fig)

        # Completeness per flux bin
        recoveredCounts, _ = np.histogram(
            matchedFakeTab["fluxJy_injected"], bins=fluxBins
        )

        injectedCounts, _ = np.histogram(
            fakeTab["fluxJy_injected"], bins=fluxBins
        )

        with np.errstate(invalid="ignore", divide="ignore"):
            fractionRecovered = np.where(
                injectedCounts > 0,
                recoveredCounts / injectedCounts,
                np.nan
            )

        allFractionRecovered.append(fractionRecovered)

    if len(allFractionRecovered) == 0:
        print("No valid injection runs found. Nothing to plot.")
        return

    allFractionRecovered = np.array(allFractionRecovered)
    meanFraction = np.nanmean(allFractionRecovered, axis=0)
    stdFraction  = np.nanstd(allFractionRecovered, axis=0)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(allBinCentres, meanFraction, color="steelblue", lw=2)
    ax.fill_between(allBinCentres,
                    meanFraction - stdFraction,
                    meanFraction + stdFraction,
                    alpha=0.3, color="steelblue")
    ax.axhline(1.0, linestyle="--", color="k", lw=1)
    ax.set_xscale("log")
    ax.set_ylim(-0.05, 1.15)
    ax.set_xlabel("Flux density S / Jy", fontsize=13)
    ax.set_ylabel("Completeness  (recovered / injected)", fontsize=13)
    title = os.path.basename(imageName)
    if sizeCfg is not None:
        title += "\n(averaged over sizes %s-uniform in %.1f-%.1f arcsec)" % (
            sizeCfg.get('spacing', 'log'), sizeCfg['minSizeArcsec'], sizeCfg['maxSizeArcsec'])
    ax.set_title(title, fontsize=10)
    ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
    plotName = os.path.join(sinjectDir, "complPlot_%s.png" %os.path.basename(imageName))
    print(plotName)
    fig.savefig(plotName, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("\nCompleteness plot saved to %s" % plotName)

    # Saving the completeness data to a text file
    completenessDataFile = os.path.join(sinjectDir, "completeness_%s.txt" %os.path.basename(imageName))
    with open(completenessDataFile, "w") as f:
        f.write("# FluxBinCentre_Jy  MeanCompleteness  StdDevCompleteness\n")
        for centre, mean, std in zip(allBinCentres, meanFraction, stdFraction):
            f.write(f"{centre:.6e} {mean:.6f} {std:.6f}\n")
    print("\nCompleteness data saved to %s" % completenessDataFile)

    if sizeCfg is not None and len(allInj2D) > 0:
        analyse2DCompleteness(allInj2D, allRec2D, fluxEdges2D, sizeEdges, sizeCfg, imageName, sinjectDir,
                              beamMajArcsec=beamMajArcsec, minCountPerCell=minCountPerCell,
                              nSizeGroups=nSizeGroups, sizePlotFluxRangesJy=sizePlotFluxRangesJy)

    return meanFraction, stdFraction, allBinCentres

def analyse2DCompleteness(allInj2D, allRec2D, fluxEdges, sizeEdges, sizeCfg, imageName, sinjectDir,
                          beamMajArcsec=None, minCountPerCell=10, nSizeGroups=4, sizePlotFluxRangesJy=None):
    """
    Completeness as a function of flux density and source size.

    Saves (i) the 2D grid (.npz and a plain text table), (ii) a heatmap of completeness in the flux density -
    size plane, (iii) completeness vs flux density in groups of source size and (iv) completeness vs size in
    ranges of flux density.

    Args:
        allInj2D, allRec2D: Lists with one (nFluxBins, nSizeBins) array per injection run, with the number of
            injected and recovered sources.
        fluxEdges, sizeEdges: Bin edges of the grid (flux density in Jy, deconvolved major-axis FWHM in arcsec).
        sizeCfg: Size configuration used for the injection.
        beamMajArcsec: Beam major axis in arcsec. If given, a second y axis shows the size in units of the beam.
        minCountPerCell: Cells with fewer injected sources (summed over runs) are masked in the heatmap.
        nSizeGroups: Number of size groups for plot (iii).
        sizePlotFluxRangesJy: List of (min, max) flux ranges (Jy) for plot (iv). Default: four octaves starting
            at half the flux density where the size-averaged completeness reaches 50%.
    """

    inj = np.asarray(allInj2D)      # (nRuns, nFlux, nSize)
    rec = np.asarray(allRec2D)
    logSize = sizeCfg.get('spacing', 'log') == 'log'

    fluxCen = np.sqrt(fluxEdges[:-1] * fluxEdges[1:])
    sizeCen = np.sqrt(sizeEdges[:-1] * sizeEdges[1:]) if logSize else 0.5 * (sizeEdges[:-1] + sizeEdges[1:])

    # Pooled completeness (sum over runs): more stable than averaging per-run fractions in sparse cells
    injTot, recTot = inj.sum(axis=0), rec.sum(axis=0)
    comp = np.full(injTot.shape, np.nan)
    good = injTot >= minCountPerCell
    comp[good] = recTot[good] / injTot[good]
    compMasked = np.ma.masked_invalid(comp)

    baseName = os.path.basename(imageName)

    # ---- save the grid
    npzFile = os.path.join(sinjectDir, "completeness2D_%s.npz" % baseName)
    np.savez(npzFile, fluxEdgesJy=fluxEdges, sizeEdgesArcsec=sizeEdges, nInjected=injTot,
             nRecovered=recTot, completeness=comp)
    txtFile = os.path.join(sinjectDir, "completeness2D_%s.txt" % baseName)
    with open(txtFile, "w") as f:
        f.write("# FluxCentre_Jy  SizeCentre_arcsec(deconvolved major FWHM)  N_injected  N_recovered  Completeness\n")
        for i in range(len(fluxCen)):
            for j in range(len(sizeCen)):
                f.write("%.6e %.4f %d %d %.6f\n" % (fluxCen[i], sizeCen[j], injTot[i, j], recTot[i, j], comp[i, j]))
    print("\n2D completeness saved to %s and %s" % (npzFile, txtFile))

    def meanStd(numer, denom):
        """Mean and scatter over runs of the per-run completeness numer/denom; arrays (nRuns, nBins)."""
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.where(denom > 0, numer / denom, np.nan)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            return np.nanmean(frac, axis=0), np.nanstd(frac, axis=0)

    # ---- (ii) heatmap
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    pcm = ax.pcolormesh(fluxEdges, sizeEdges, compMasked.T, cmap="viridis", vmin=0.0, vmax=1.0,
                        shading="flat")
    cbar = fig.colorbar(pcm, ax=ax, pad=0.14)
    cbar.set_label("Completeness  (recovered / injected)", fontsize=12)
    try:
        cs = ax.contour(fluxCen, sizeCen, compMasked.T, levels=[0.5, 0.9], colors=["w", "k"], linewidths=1.3)
        ax.clabel(cs, fmt="%.1f", fontsize=9)
    except Exception as e:
        print("Could not draw completeness contours: %s" % e)
    ax.set_xscale("log")
    if logSize:
        ax.set_yscale("log")
    ax.set_xlabel("Injected flux density S / Jy", fontsize=13)
    ax.set_ylabel("Deconvolved major axis (FWHM) / arcsec", fontsize=13)
    ax.set_title(baseName, fontsize=10)
    if beamMajArcsec is not None and beamMajArcsec > 0:
        secAx = ax.secondary_yaxis("right", functions=(lambda x: x / beamMajArcsec,
                                                       lambda x: x * beamMajArcsec))
        secAx.set_ylabel(r"$\theta_{\rm maj}$ / beam major axis", fontsize=12)
    heatName = os.path.join(sinjectDir, "complHeatmap_%s.png" % baseName)
    fig.savefig(heatName, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Completeness heatmap saved to %s" % heatName)

    # ---- (iii) completeness vs flux density in groups of size
    groups = np.array_split(np.arange(len(sizeCen)), nSizeGroups)
    colours = plt.cm.viridis(np.linspace(0.0, 0.9, len(groups)))
    fig, ax = plt.subplots(figsize=(8, 6))
    for idx, col in zip(groups, colours):
        if len(idx) == 0:
            continue
        mean, std = meanStd(rec[:, :, idx].sum(axis=2), inj[:, :, idx].sum(axis=2))
        lab = "%.1f - %.1f arcsec" % (sizeEdges[idx[0]], sizeEdges[idx[-1] + 1])
        ax.plot(fluxCen, mean, color=col, lw=2, label=lab)
        ax.fill_between(fluxCen, mean - std, mean + std, color=col, alpha=0.2)
    mean, std = meanStd(rec.sum(axis=2), inj.sum(axis=2))
    ax.plot(fluxCen, mean, color="k", lw=2, ls="--", label="all sizes")
    ax.axhline(1.0, linestyle=":", color="k", lw=1)
    ax.set_xscale("log")
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("Injected flux density S / Jy", fontsize=13)
    ax.set_ylabel("Completeness  (recovered / injected)", fontsize=13)
    ax.legend(title="Deconvolved major axis", frameon=False)
    ax.set_title(baseName, fontsize=10)
    ax.xaxis.set_major_formatter(ticker.LogFormatterSciNotation())
    byName = os.path.join(sinjectDir, "complVsFluxBySize_%s.png" % baseName)
    fig.savefig(byName, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Completeness vs flux density in size groups saved to %s" % byName)

    # ---- (iv) completeness vs size in ranges of flux density
    if sizePlotFluxRangesJy is None:
        allFrac = rec.sum(axis=(0, 2)) / np.maximum(inj.sum(axis=(0, 2)), 1.0)   # size-averaged, pooled
        above = np.where(allFrac >= 0.5)[0]
        if len(above) == 0:
            print("Size-averaged completeness never reaches 50%: skipping completeness vs size plot.")
            return
        i = above[0]
        if i == 0 or allFrac[i] == allFrac[i - 1]:
            f50 = fluxCen[i]
        else:
            lf = np.log10(fluxCen[i - 1:i + 1])
            c0, c1 = allFrac[i - 1], allFrac[i]
            f50 = 10 ** (lf[0] + (0.5 - c0) * (lf[1] - lf[0]) / (c1 - c0))
        sizePlotFluxRangesJy = [(f50 * 2.0 ** k, f50 * 2.0 ** (k + 1)) for k in range(-1, 3)]
        print("Size-averaged completeness reaches 50%% at S = %.3g Jy" % f50)

    colours = plt.cm.plasma(np.linspace(0.0, 0.85, len(sizePlotFluxRangesJy)))
    fig, ax = plt.subplots(figsize=(8, 6))
    for (lo, hi), col in zip(sizePlotFluxRangesJy, colours):
        sel = (fluxCen >= lo) & (fluxCen < hi)
        if sel.sum() == 0:
            continue
        mean, std = meanStd(rec[:, sel, :].sum(axis=1), inj[:, sel, :].sum(axis=1))
        ax.plot(sizeCen, mean, "-o", ms=3, color=col, lw=1.8,
                label="%.3g - %.3g mJy" % (lo * 1e3, hi * 1e3))
        ax.fill_between(sizeCen, mean - std, mean + std, color=col, alpha=0.2)
    if logSize:
        ax.set_xscale("log")
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("Deconvolved major axis (FWHM) / arcsec", fontsize=13)
    ax.set_ylabel("Completeness  (recovered / injected)", fontsize=13)
    ax.legend(title="Injected flux density", frameon=False)
    ax.set_title(baseName, fontsize=10)
    sizeName = os.path.join(sinjectDir, "complVsSize_%s.png" % baseName)
    fig.savefig(sizeName, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("Completeness vs size saved to %s" % sizeName)

def executeSingle(imageName, nInjectionSources=5000, nRepetitions=1, minFluxJyInj=1e-5, maxFluxJyInj=1.0,
                  radiusFactorToInject=0.95, outDir=None, injectUnresolvedOnly=False,
                  minSizeArcsec=1.0, maxSizeArcsec=60.0, minAxisRatio=0.3, crowdingFactor=5.0):
    """
    Executes the source injection and completeness analysis for a single image.

    Args:
    - imageName: Name of the image file to process.
    - nInjectionSources: Number of sources to inject (default: 10000).
    - nRepetitions: Number of repetitions for the injection (default: 100).
    - minFluxJyInj: Minimum flux density of injected sources in Jy (default: 1e-5).
    - maxFluxJyInj: Maximum flux density of injected sources in Jy (default: 1.0).
    - radiusFactorToInject: Factor to multiply the band radius for defining the injection area (default: 0.95).
    - outDir: Path to the directory to save source injection files. Default: current working directory
    - injectUnresolvedOnly: If True, only unresolved sources will be injected. If False, sources with random sizes are injected (default: False).
    - minSizeArcsec, maxSizeArcsec: range of the injected deconvolved major-axis FWHM in arcsec; sizes are drawn log-uniformly (default: 1-60).
    - minAxisRatio: minor/major axis ratio is drawn uniformly in [minAxisRatio, 1] (default: 0.3). The position angle is uniform in [0, 180) deg.
    - crowdingFactor: minimum separation between injected sources in units of their mean observed FWHM (default: 5).
    """

    currentDir = os.getcwd()
    imageBaseName = os.path.basename(imageName).split(".fits")[0]
    catBaseName = os.path.basename(imageName).split(".")[0]

    localProdDir  = startup.config['productsDir']

    needExtraction = False

    pybdsfCatFileName = "%s_srl_bdsfcat.fits" %catBaseName
    residualFileName = "%s_gaus_resid.fits" %catBaseName
    rmsFileName  = "%s_rms.fits" %catBaseName
    meanFileName = "%s_mean.fits" %catBaseName

    pybdsfCatFilePath = None

    if os.path.exists(os.path.join(currentDir, pybdsfCatFileName)):
        pybdsfCatFilePath = os.path.join(currentDir, pybdsfCatFileName)
        print("\nUsing existing PyBDSF catalogue from current directory: %s\n" % pybdsfCatFilePath)
    elif os.path.exists(os.path.join(localProdDir, "catalogs", pybdsfCatFileName)):
        pybdsfCatFilePath = os.path.join(localProdDir, "catalogs", pybdsfCatFileName)
        print("\nUsing existing PyBDSF catalogue from berk products directory: %s\n" % pybdsfCatFilePath)
    else:
        needExtraction = True

    if os.path.exists(os.path.join(currentDir, residualFileName)):
        residualFilePath = os.path.join(currentDir, residualFileName)
        print("\nUsing existing residual image from current directory: %s\n" % residualFilePath)
    elif os.path.exists(os.path.join(localProdDir, "residual", residualFileName)):
        residualFilePath = os.path.join(localProdDir, "residual", residualFileName)
        print("\nUsing existing residual image from berk products directory: %s\n" % residualFilePath)
    else:
        needExtraction = True

    if os.path.exists(os.path.join(currentDir, rmsFileName)):
        rmsFilePath = os.path.join(currentDir, rmsFileName)
        print("\nUsing existing RMS map from current directory: %s\n" % rmsFilePath)
    elif os.path.exists(os.path.join(localProdDir, "rms",  rmsFileName)):
        rmsFilePath  = os.path.join(localProdDir, "rms",  rmsFileName)
        print("\nUsing existing RMS map from berk products directory: %s\n" % rmsFilePath)
    else:
        needExtraction = True

    if os.path.exists(os.path.join(currentDir, meanFileName)):
        meanFilePath = os.path.join(currentDir, meanFileName)
        print("\nUsing existing mean map from current directory: %s\n" % meanFilePath)
    elif os.path.exists(os.path.join(localProdDir, "mean", meanFileName)):
        meanFilePath = os.path.join(localProdDir, "mean", meanFileName)
        print("\nUsing existing mean map from berk products directory: %s\n" % meanFilePath)
    else:
        needExtraction = True

    if needExtraction is True:
        try:
            sourceExtract(imageFile=imageName, pybdsfBaseName=imageBaseName, saveRMSMeanResidualMaps=True, outDir=currentDir, rmsFileName=None, meanFileName=None)
            pybdsfCatFilePath = os.path.join(currentDir, pybdsfCatFileName)
            residualFilePath = os.path.join(currentDir, residualFileName)
            rmsFilePath = os.path.join(currentDir, rmsFileName)
            meanFilePath = os.path.join(currentDir, meanFileName)
        except Exception as e:
            print("Error during source extraction of %s: %s" % (imageName, e))
            print("Skipping this file and moving on...")
            return

    else:
        print("\nAll necessary files found. Proceeding with injection and completeness analysis...\n")

    if pybdsfCatFilePath is not None:
        pybdsfCat = Table.read(pybdsfCatFilePath, format='fits')
    else:
        print("Error: PyBDSF catalogue file not found. Cannot proceed.")
        return

    tabWrapAngle = catalogs.detectWrapAngle(pybdsfCat['RA'])
    pybdsfCat = catalogs.fixRA(pybdsfCat, raCol='RA', wrapAngle=tabWrapAngle)

    if injectUnresolvedOnly is True:
        sourceType = 'unresolved'
    else:
        sourceType = 'extended'

    if outDir is None:
        outDir = os.getcwd()
    sinjectDir = os.path.join(outDir, "sinjected_%s_%s_%dsources_%dreps" %(sourceType, imageBaseName, nInjectionSources, nRepetitions))

    os.makedirs(sinjectDir, exist_ok=True)

    # Symlink the residual and rms/mean maps into sinjectDir
    for srcPath, dstName in [
        (residualFilePath, os.path.basename(residualFilePath)),
        (rmsFilePath,  os.path.basename(rmsFilePath)),
        (meanFilePath, os.path.basename(meanFilePath)),
    ]:
        dst = os.path.join(sinjectDir, dstName)
        if not os.path.exists(dst):
            try:
                os.symlink(srcPath, dst)
            except OSError as e:
                print("Symlink warning: %s" % e)

    # Inject fake sources into the residual image
    fieldRACentre, fieldDecCentre, bandRadius = crossmatch.getCentreRadiusFromCatalog(
        pybdsfCat, radRACol="RA", radDecCol="DEC"
    )
    #fieldRACentre, fieldDecCentre, bandRadius = crossmatch.getCentreRadiusFromImagesTab(pybdsfCatFilePath)

    # Inject within certain % of the band radius to avoid edge effects
    radiusToInject = bandRadius * radiusFactorToInject

    if injectUnresolvedOnly is True:
        sizeCfg = None   # point sources only
    else:
        sizeCfg = {
            'minSizeArcsec': minSizeArcsec,
            'maxSizeArcsec': maxSizeArcsec,
            'spacing': 'log',
            'minAxisRatio': minAxisRatio,
        }

    injectImage(
        residualFileName,
        fieldRACentre, fieldDecCentre, radiusToInject,
        sinjectDir, imageBaseName,
        nInjectionSources=nInjectionSources,
        minFluxJyInj=minFluxJyInj,
        maxFluxJyInj=maxFluxJyInj,
        nRepetitions=nRepetitions,
        sizeCfg=sizeCfg,
        crowdingFactor=crowdingFactor,
    )

    # Run PyBDSF on each injected image
    sInjectedFileList = sorted(glob.glob(
        os.path.join(sinjectDir, "%s_sinjected_*_image.fits" % imageBaseName)
    ))

    for sInjectedImageFile in sInjectedFileList:
        sInjectedPybdsfBaseName = os.path.basename(sInjectedImageFile).replace("_image.fits", "")
        sourceExtract(
            sInjectedImageFile,
            sInjectedPybdsfBaseName,
            saveRMSMeanResidualMaps=False,
            outDir=sinjectDir,
            rmsFileName=rmsFileName,
            meanFileName=meanFileName,
        )

    # Analyse completeness
    sInjectedCatList = sorted(glob.glob(
        os.path.join(sinjectDir, "*_srl_bdsfcat.fits")
    ))

    print("\nFound %d injection catalogues to analyse.\n" % len(sInjectedCatList))

    calculateCompleteness(
        sInjectedCatList,
        pybdsfCat,
        imageName,
        sinjectDir,
        minFluxJyInj,
        maxFluxJyInj,
        nJyBins=100,
        sizeCfg=sizeCfg,
    )


def execute(imageName=None, nInjectionSources=5000, nRepetitions=1, minFluxJyInj=1e-5, maxFluxJyInj=1.0, radiusFactorToInject=0.95, injectUnresolvedOnly=False,
            minSizeArcsec=1.0, maxSizeArcsec=60.0, minAxisRatio=0.3, crowdingFactor=5.0):
    """
    Main function to execute the source injection and completeness analysis.

    Args:
    - imageName: Name of the image file to process. If None, all '*-image.fits' files in current directory will be processed.
    - pybdsfCatFilePath: Path to the PyBDSF catalogue file. If None, it will be searched in the current directory and then in the local products directory.
    - rmsFilePath: Path to the RMS map file. If None, it will be searched in the current directory and then in the local products directory.
    - meanFilePath: Path to the mean map file. If None, it will be searched in the current directory and then in the local products directory.
    - residualFilePath: Path to the residual image file. If None, it will be searched in the current directory and then in the local products directory.
    - nInjectionSources: Number of sources to inject (default: 10000).
    - nRepetitions: Number of repetitions for the injection (default: 100).
    - minFluxJyInj: Minimum flux density of injected sources in Jy (default: 1e-5).
    - maxFluxJyInj: Maximum flux density of injected sources in Jy (default: 1.0).
    - radiusFactorToInject: Factor to multiply the band radius for defining the injection area (default: 0.95).
    - injectUnresolvedOnly: If True, only unresolved sources will be injected. If False, sources with random sizes are injected (default: False).
    - minSizeArcsec, maxSizeArcsec: range of the injected deconvolved major-axis FWHM in arcsec; sizes are drawn log-uniformly (default: 1-60).
    - minAxisRatio: minor/major axis ratio is drawn uniformly in [minAxisRatio, 1] (default: 0.3). The position angle is uniform in [0, 180) deg.
    - crowdingFactor: minimum separation between injected sources in units of their mean observed FWHM (default: 5).
    """

    if imageName is not None:
        executeSingle(imageName, nInjectionSources, nRepetitions, minFluxJyInj, maxFluxJyInj, radiusFactorToInject, None, injectUnresolvedOnly, minSizeArcsec=minSizeArcsec, maxSizeArcsec=maxSizeArcsec, minAxisRatio=minAxisRatio, crowdingFactor=crowdingFactor)
        return

    print("No imageName provided. Searching for FITS images in current directory...\n")

    currentDir = os.getcwd()
    imageFilesInCwd = glob.glob(os.path.join(currentDir, "*-image.fits"))

    if len(imageFilesInCwd) == 0:
        print("\nNo '*-image.fits' file found in current directory. Exiting.")
        return

    print("\nFound %d image files:\n" % len(imageFilesInCwd))

    for imageFile in imageFilesInCwd:
        print("=" * 60)
        print("Working on %s" % imageFile)

        try:
            executeSingle(imageFile, nInjectionSources, nRepetitions, minFluxJyInj, maxFluxJyInj, radiusFactorToInject, None, injectUnresolvedOnly, minSizeArcsec=minSizeArcsec, maxSizeArcsec=maxSizeArcsec, minAxisRatio=minAxisRatio, crowdingFactor=crowdingFactor)
        except Exception as e:
            print("Error processing %s: %s" % (imageFile, e))
            print("Skipping this file and moving on...\n")
            continue