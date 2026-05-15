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
     sinjectDir, imageBaseName, rCutOff) = args
    
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

        # to keep track of previously injected sources to avoid crowding
        injectedCoords = []

        nSourcesInjected = 0

        randomRAInjected = []
        randomDecInjected = []
        randomFluxInjected= []

        while nSourcesInjected < nInjectionSources:

            randomRA, randomDec = crossmatch.randomPointsInCircleExactNPoints(fieldRACentre, fieldDecCentre, injectRadius, 1, rng=rng)

            randPostCoord = SkyCoord(ra=randomRA[0]*u.deg, dec=randomDec[0]*u.deg)
            if len(injectedCoords) > 0:
                prevCoords = SkyCoord(ra=np.array(injectedCoords)[:,0]*u.deg, dec=np.array(injectedCoords)[:,1]*u.deg)
                seps = randPostCoord.separation(prevCoords)
                # #TODO: harcoded 5 times the the beam size
                if seps.min() < (2 * 6.0 * u.arcsec):
                    #print("Source at RA=%.4f, Dec=%.4f is too close to a previously injected source. Skipping." % (randomRA[0], randomDec[0]))
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

            randomFlux = 10 ** rng.uniform(logMin, logMax)

            localY, localX = np.mgrid[yMin:yMax, xMin:xMax]

            psfStamp = beam(localX - x0, localY - y0)

            copyImg[yMin:yMax, xMin:xMax] += randomFlux * psfStamp

            nSourcesInjected += 1
            randomRAInjected.append(float(randomRA[0]))
            randomDecInjected.append(float(randomDec[0]))
            randomFluxInjected.append(float(randomFlux))
            
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
                minFluxJyInj=1E-5, maxFluxJyInj=1.0, nRepetitions=100):
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
    rCutOff = int(5 * max(sigmaX, sigmaY)) + 1

    argsList = [
        (rep, img, hdr, beam, nInjectionSources, fieldRACentre, fieldDecCentre, injectRadius,
         minFluxJyInj, maxFluxJyInj, sinjectDir, imageBaseName, rCutOff)
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
    matchRadiusAngle = matchRadDeg * u.deg
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
                        minFluxJyInj, maxFluxJyInj, nJyBins=30):
    """
    Calculates the completeness (fraction of injected sources recovered) as a function of flux density, and also the flux recovery ratio.
    Args:
        sInjectedCatList: List of file paths to the PyBDSF catalogues of the injected images.
        pybdsfCat: The PyBDSF catalogue of the original image (without injections).
        imageName: Name of the image being processed (used for plot titles).
        sinjectDir: Directory where the injection results are stored (used for saving plots).
        minFluxJyInj: Minimum flux density of injected sources in Jy.
        maxFluxJyInj: Maximum flux density of injected sources in Jy.
        nJyBins: Number of flux bins to use for completeness calculation.
        
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

    for catIdx, sInjectedCatFile in enumerate(sInjectedCatList):

        print("\nProcessing %s (%d/%d)" % (sInjectedCatFile, catIdx+1, len(sInjectedCatList)))

        fakeSourceFile = sInjectedCatFile.replace("srl_bdsfcat", "fakesources")
        if not os.path.exists(fakeSourceFile):
            print("Fake-source file not found - skipping.")
            continue

        recovTab    = Table.read(sInjectedCatFile)
        fakeTab     = Table.read(fakeSourceFile)

        if any(recovTab['RA'] < 0.0):
            recovTab = catalogs.fixRA(recovTab, raCol='RA', wrapAngle=360)

        if len(recovTab) == 0:
            print("No sources recovered by PyBDSF - skipping.")
            allFractionRecovered.append(np.zeros(nJyBins))
            continue

        try:
            beamMajor, beamMinor = getBeamParamsFromFits(sInjectedCatFile)
        except ValueError as e:
            print("Error: %s - skipping." % e)
            continue

        # in case we want to inject sources to an inner area of the image #TODO
        matchRadDeg = max(beamMajor, beamMinor)

        # injectedRecoveredMask : boolean mask for the injected source table that were recovered by PyBDSF
        # matchedFakeTab : subset of the injected source table that were recovered by PyBDSF
        injectedRecoveredMask, matchedFakeTab = matchFakeToRecovered(fakeTab, recovTab, matchRadDeg, fluxTolerance=None)

        nRecovered = injectedRecoveredMask.sum()
        print("%d / %d injected fake sources recovered" % (nRecovered, len(fakeTab)))

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
    ax.set_title(os.path.basename(imageName), fontsize=10)
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
    
    return meanFraction, stdFraction, allBinCentres

def executeSingle(imageName, nInjectionSources=5000, nRepetitions=1, minFluxJyInj=1e-5, maxFluxJyInj=1.0, radiusFactorToInject=1.0, outDir=None):
    """
    Executes the source injection and completeness analysis for a single image.
    
    Args:       
    - imageName: Name of the image file to process.
    - nInjectionSources: Number of sources to inject (default: 10000).
    - nRepetitions: Number of repetitions for the injection (default: 100).
    - minFluxJyInj: Minimum flux density of injected sources in Jy (default: 1e-5).
    - maxFluxJyInj: Maximum flux density of injected sources in Jy (default: 1.0).
    - radiusFactorToInject: Factor to multiply the band radius for defining the injection area (default: 1.0).
    - outDir: Path to the directory to save source injection files. Default: current working directory
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

    if any(pybdsfCat['RA'] < 0.0):
        pybdsfCat = catalogs.fixRA(pybdsfCat, raCol='RA', wrapAngle=360)

    if outDir is None:
        outDir = os.getcwd()
    sinjectDir = os.path.join(outDir, "sinjected_%s_%dsources_%dreps" %(imageBaseName, nInjectionSources, nRepetitions))

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

    injectImage(
        residualFileName,
        fieldRACentre, fieldDecCentre, radiusToInject,
        sinjectDir, imageBaseName,
        nInjectionSources=nInjectionSources,
        minFluxJyInj=minFluxJyInj,
        maxFluxJyInj=maxFluxJyInj,
        nRepetitions=nRepetitions,
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
    )


def execute(imageName=None, nInjectionSources=5000, nRepetitions=1, minFluxJyInj=1e-5, maxFluxJyInj=1.0, radiusFactorToInject=1.0):
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
    - radiusFactorToInject: Factor to multiply the band radius for defining the injection area (default: 1.0).
    """
    
    if imageName is not None:
        executeSingle(imageName, nInjectionSources, nRepetitions, minFluxJyInj, maxFluxJyInj, radiusFactorToInject)
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
            executeSingle(imageFile, nInjectionSources, nRepetitions, minFluxJyInj, maxFluxJyInj, radiusFactorToInject)
        except Exception as e:
            print("Error processing %s: %s" % (imageFile, e))
            print("Skipping this file and moving on...\n")
            continue
