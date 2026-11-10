"""

Tasks that can be performed by Berk

"""

import os
import sys
import subprocess
import glob
import datetime
import astropy.table as atpy
from astropy.io import fits
from astropy.table import Column
from . import startup, jobs, catalogs, images,  __version__, crossmatch, summaryPlots, sourceinjection
import shlex
from astropy.coordinates import SkyCoord
import astropy.units as u
import random
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from multiprocessing import Pool
from scipy.interpolate import PchipInterpolator

#------------------------------------------------------------------------------------------------------------
def fetch(captureBlockId):
    """Fetch...

    Probably need to run ``berk fetch -o rdb_link`` in GNU screen for this task.

    """

    captureBlockIdLink=captureBlockId
    captureBlockId=captureBlockIdLink.split("https://archive-gw-1.kat.ac.za/")[-1].split("/")[0]
    msPath=os.environ['BERK_MSCACHE']+os.path.sep+"%s_sdp_l0.ms" % (captureBlockId)
    fetchLogPath=os.environ['BERK_MSCACHE']+os.path.sep+"%s_fetch.log" % (captureBlockId)

    if os.path.exists(msPath) is True:
        print("Already fetching %s - if it failed, remove %s and try again" % (captureBlockId, msPath))
        sys.exit()

    cmd="mvftoms.py %s --flags cam,data_lost,ingest_rfi -o %s" % (captureBlockIdLink, msPath)
    cmdSubprocess=shlex.split(cmd)
    with subprocess.Popen(cmdSubprocess, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        with open(fetchLogPath, 'w') as logFile:
            for line in process.stdout:
                logFile.write(line)
                print(line, end='')

            for line in process.stderr:
                logFile.write(line)
                print(line, end='')
        process.wait()

    # Automatic setting of KATSDPTELSTATE_ALLOW_PICKLE=1
    katsString="you can allow the pickles to be loaded by setting KATSDPTELSTATE_ALLOW_PICKLE=1 in the environment"
    with open(fetchLogPath, 'r') as file:
        fetchLogMsg = file.read()
        if katsString in fetchLogMsg:
            if os.path.exists(msPath) is True:
                os.rmdir(msPath)
            os.environ["KATSDPTELSTATE_ALLOW_PICKLE"] = "1"
            with subprocess.Popen(cmdSubprocess, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
                with open(fetchLogPath, 'w') as logFile:
                    for line in process.stdout:
                        logFile.write(line)
                        print(line, end='')

                    for line in process.stderr:
                        logFile.write(line)
                        print(line, end='')
                process.wait()

    sys.exit()

#------------------------------------------------------------------------------------------------------------
def store(captureBlockId):
    """Store...

    """

    print("Task 'store' is not implemented yet.")
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def getBandKey(freqGHz):
    """Return band name ('L' or 'UHF') based on freqGHz.

    """
    if freqGHz >= 0.9 and freqGHz < 1.7:
        bandKey='L'
    elif freqGHz >= 0.6 and freqGHz < 0.9:
        bandKey='UHF'
    elif freqGHz >= 1.7 and freqGHz < 3.5:
        bandKey='S'
    else:
        raise ValueError("%f - Not sure what band this is - need to add a band key for it" %freqGHz)

    return bandKey

#------------------------------------------------------------------------------------------------------------
def listObservations():
    """List observations available on this machine, and check their processing status with the central list.

    """

    if 'BERK_INFO_FILE' in os.environ.keys():
        tab=atpy.Table().read(os.environ['BERK_INFO_FILE'])
    else:
        print("Set BERK_INFO_FILE environment variable to check processing status of observations against central list.")
        tab=None

    # MJH: Don't hard code in things like this, it is broken for me
    globXmatchTab=None
    # globXmatchLink = 'https://dl.dropbox.com/scl/fi/39gma12peyymhadc2fss2/xmatchCat_DECaLS_r_4p0asec.fits?rlkey=v4bazyqfy294ukt9eqmfaeed7&st=7e7nfcde&dl=0'
    # globXmatchTab=atpy.Table().read(globXmatchLink)

    msList=glob.glob(os.environ['BERK_MSCACHE']+os.path.sep+"*_sdp_l0.ms")
    msList.sort()

    print("Downloaded observations available locally [by captureBlockId]:")
    for ms in msList:
        captureBlockId=os.path.split(ms)[-1].split("_")[0]
        status="cached"
        # Test for setup done
        globStr=os.environ['BERK_ROOT']+os.path.sep+'processing'+os.path.sep+captureBlockId+os.path.sep+'project_info.json'
        if len(glob.glob(globStr)) > 0:
            status="setup"
        # Test for 1GC done
        globStr=os.environ['BERK_ROOT']+os.path.sep+'processing'+os.path.sep+captureBlockId+os.path.sep+'%s_sdp_*_1024ch_*.ms' % (captureBlockId)
        if len(glob.glob(globStr)) > 0:
            status="process1"
        # Test for 2GC done
        globStr=os.environ['BERK_ROOT']+os.path.sep+'processing'+os.path.sep+captureBlockId+os.path.sep+'IMAGES'+os.path.sep+'img_%s_*_pcalmask-MFS-image.fits' % (captureBlockId)
        if len(glob.glob(globStr)) > 0:
            status="process2"
        # Below only shows up _after_ we run 'analyse', then 'collect' and 'builddb' which updates the BERK_INFO_FILE table
        if tab is not None:
            if captureBlockId in tab['captureBlockId']:
                status="processed_and_analysed"
        if globXmatchTab is not None:
            if captureBlockId in globXmatchTab['captureBlockId']:
                status="processed_analysed_xmatched-DECaLS"

        print("   %s    %s" % (captureBlockId, status))

def samplePhotoZ68p95(zmed, zl68, zu68, zl95, zu95, nReal, rng=None, clip_min=1e-6):
    """
    Draw Monte Carlo samples from an approximate photometric-redshift PDF
    defined by the median and the 68% and 95% confidence intervals.

    Args:
        zmed (float): Median photometric redshift.
        zl68 (float): Lower 68% confidence limit.
        zu68 (float): Upper 68% confidence limit.
        zl95 (float): Lower 95% confidence limit.
        zu95 (float): Upper 95% confidence limit.
        nReal (int): Number of Monte Carlo samples.
        rng (numpy.random.Generator, optional): Random number generator.
        clip_min (float or None, optional): Minimum allowed redshift.

    Returns:
        numpy.ndarray: Array of sampled redshifts.
    """

    if rng is None:
        rng = np.random.default_rng()

    probs = np.array([0.025, 0.16, 0.50, 0.84, 0.975])
    zvals = np.array([zl95, zl68, zmed, zu68, zu95])

    # Monotonic interpolation of the inverse CDF
    invCDF = PchipInterpolator(probs, zvals, extrapolate=True)

    u = rng.uniform(0.025, 0.975, nReal)
    samples = invCDF(u)

    if clip_min is not None:
        samples = np.clip(samples, clip_min, None)

    return samples
#------------------------------------------------------------------------------------------------------------
def samplePhotoZ68(zmed, zl68, zu68, nReal, rng=None, clip_min=1e-6):
    """
    Draws Monte Carlo samples from an asymmetric (split-normal) approximation
    to a photometric redshift probability distribution.

    Args:
        zmed (float): Median photometric redshift.
        zl68 (float): Lower 68% confidence limit of the photometric redshift.
        zu68 (float): Upper 68% confidence limit of the photometric redshift.
        nReal (int): Number of Monte Carlo redshift realisations to generate.
        rng (numpy.random.Generator, optional): NumPy random number generator.
            If None, a new default generator is created.
        clip_min (float or None, optional): Minimum allowed sampled redshift.
            Default 1e-6.
    """

    if rng is None:
        rng = np.random.default_rng()

    sigmaLow = np.clip(zmed - zl68, 1e-6, None)
    sigmaHigh = np.clip(zu68 - zmed, 1e-6, None)

    normalDeviates = rng.standard_normal(size=nReal)

    samples = np.where(normalDeviates < 0, zmed + normalDeviates * sigmaLow, zmed + normalDeviates * sigmaHigh)

    if clip_min is not None:
        samples = np.clip(samples, clip_min, None)

    return samples

#------------------------------------------------------------------------------------------------------------
def processVmaxForXmatchFile(args):
    """
    Computes Vmax for all sources in a single cross-match catalogue file.

    Args:
        args (tuple): A tuple containing the following elements:
            - xmatchFile (str): Full path to the cross-match FITS file.
            - DRImages (astropy.table.Table): Table of DR images with RMS and sky area info.
            - DRDir (str): Path to the data release directory.
            - cosmologyDR (astropy.cosmology instance): Cosmology to use for calculations.

    Returns:
        tuple: A pair (xmatchTab, zmaxTab) as astropy Tables, or (None, None) if the
            zmax file already exists or the image row is not found.
    """
    xmatchFile, DRImages, DRDir, cosmologyDR, zmaxDirName, zColName, radCatName = args

    xmatchTab = atpy.Table().read(xmatchFile)

    zMaxDirPath = os.path.join(DRDir, zmaxDirName)
    os.makedirs(zMaxDirPath, exist_ok = True)

    zmaxFilePath = os.path.join(zMaxDirPath, 'zmax_%s.fits' %_normaliseCatName(radCatName))

    if os.path.exists(zmaxFilePath):
        zmaxTab = atpy.Table().read(zmaxFilePath)
        return xmatchTab, zmaxTab

    print("\nProcessing Vmax for %s..." % xmatchFile)

    if 'radCatName' not in xmatchTab.colnames:
        xmatchTab['radCatName'] = 'catalogs'+os.path.sep+radCatName

    if zColName not in xmatchTab.colnames:
        xmatchTab[zColName] = np.where(np.isfinite(xmatchTab['zspec_opt']) &
                                                    (xmatchTab['zspec_opt'] != -99),
                                                    xmatchTab['zspec_opt'],
                                                    xmatchTab['zphoto_median_opt'])

    zmaxTab = xmatchTab[xmatchTab[zColName] != -99]

    maxRedshift = np.max(zmaxTab[zColName])
    minRedshift = np.min(zmaxTab[zColName])

    zmaxTab['LuminosityWHz_rad'] = catalogs.calculateRadioLum(zmaxTab['Total_flux_rad'].value, zmaxTab[zColName].value,
                                                              spectralIndex=0.7, cosmology=cosmologyDR)

    # Finding effective area based on RMS area coverage curve

    rmsAreaFileName = _normaliseCatName(radCatName) + '_rms_rmshist.txt'
    rmsAreaFile = os.path.join(DRDir, 'rms', rmsAreaFileName)
    rmsHistFromFile = np.loadtxt(rmsAreaFile)
    rmsBinsCentres, cumArea = rmsHistFromFile[:, 0], rmsHistFromFile[:, 4]
    zmaxTab['EffectiveAreaSqDeg'] = summaryPlots.getEffectiveAreaInFluxBinsfromRMS(zmaxTab['Total_flux_rad'].value,
                                                   rmsBinsCentres, cumArea, sigmaDetection=5.0)
    zmaxTab['WholeAreaSqDeg'] = np.max(cumArea)

    # radioCatPath may hold the catalogue name with or without 'srl' (older images tables): compare base names
    DRImageRow = DRImages[np.array([_normaliseCatName(p) == _normaliseCatName(radCatName) for p in DRImages['radioCatPath']])]

    if len(DRImageRow) == 0:
        return None, None

    # Radio completeness: as a function of flux density and size if the 2D completeness is available,
    # otherwise as a function of flux density only
    radCatBaseName = _normaliseCatName(radCatName)
    completenessDir = _findSourceInjectionDir(os.path.join(DRDir, 'sinjection'), radCatBaseName)
    sizeMajArcsec = _getSizeArcsec(zmaxTab, [('DC_Maj_rad', 3600.0), ('sizeMajAsec_rad', 1.0), ('sizeMajAsec', 1.0)])
    completenessSources, completenessType = assignRadioCompleteness(zmaxTab['Total_flux_rad'], sizeMajArcsec,
                                                                    completenessDir, radCatBaseName)
    zmaxTab['completeness_rad'] = completenessSources
    zmaxTab['completenessType_rad'] = np.int16(completenessType)   # 2: flux + size, 1: flux only, 0: not available
    if completenessType == 0:
        print("Radio completeness is not available for %s" % radCatName)

    # Calculating optical completeness
    # Given by the fraction of radio sources in a flux bin
    # that has an optical counterpart

    radioCatTab = atpy.Table().read(os.path.join(DRDir, 'catalogs',radCatName))

    fluxBinsForOptComp = np.logspace(np.log10(np.nanmin(zmaxTab['Total_flux_rad'])), np.log10(np.nanmax(zmaxTab['Total_flux_rad'])), 35)
    radioSourcesHist, _ = np.histogram(radioCatTab['Total_flux'], bins=fluxBinsForOptComp)
    xmatchedSourcesHist, _ = np.histogram(zmaxTab['Total_flux_rad'], bins=fluxBinsForOptComp)
    optCompletenessInFluxBins = np.divide(xmatchedSourcesHist, radioSourcesHist, out=np.zeros_like(xmatchedSourcesHist, dtype=float),where=radioSourcesHist>0)

    RMS_uJypbeam = DRImageRow['RMS_uJy/beam'][0]
    skyArea_sqDeg = DRImageRow['skyArea_sqDeg'][0]

    sigmaTimes = 5 # 5 sigma detection limit

    Slim_Jypbeam = sigmaTimes*RMS_uJypbeam*1e-6 # flux limit in Jy/beam

    zMaxList = []
    VMaxList = []
    VMaxEffAreaList = []
    VMaxWholeAreaList = []
    optCompList = []

    for galNow in zmaxTab:

        indexFlux = (np.digitize(galNow['Total_flux_rad'], fluxBinsForOptComp) - 1)
        indexFlux = np.clip(indexFlux, 0, len(optCompletenessInFluxBins)-1)
        optCompNow = optCompletenessInFluxBins[indexFlux]
        optCompList.append(optCompNow)

        galRedshift = galNow[zColName]
        galLum_WHz = galNow['LuminosityWHz_rad']
        galEffArea = galNow['EffectiveAreaSqDeg']
        galWholeArea = galNow['WholeAreaSqDeg']
        Slim_uJypbeam = Slim_Jypbeam*1E6

        zMax = catalogs.calculateZmax(galRedshift, galLum_WHz, Slim_uJypbeam, alpha=0.7, zmaxLimit=1E10,cosmology=cosmologyDR)
        VMax_h3Mpc3 = catalogs.calculateComovingVolBetweenZ_h3Mpc3(skyArea_sqDeg, zMin=0.0, zMax=zMax, cosmology=cosmologyDR)
        VMax_h3Mpc3_effarea = catalogs.calculateComovingVolBetweenZ_h3Mpc3(galEffArea, zMin=0.0, zMax=zMax, cosmology=cosmologyDR)
        VMax_h3Mpc3_wholearea = catalogs.calculateComovingVolBetweenZ_h3Mpc3(galWholeArea, zMin=0.0, zMax=zMax, cosmology=cosmologyDR)

        zMaxList.append(zMax)
        VMaxList.append(VMax_h3Mpc3)
        VMaxEffAreaList.append(VMax_h3Mpc3_effarea)
        VMaxWholeAreaList.append(VMax_h3Mpc3_wholearea)

    zmaxTab['zmax_rad'] = zMaxList
    zmaxTab['Vmax_rad_h3Mpc3'] = VMaxList
    zmaxTab['Vmax_rad_h3Mpc3_effarea'] = VMaxEffAreaList
    zmaxTab['Vmax_rad_h3Mpc3_wholearea'] = VMaxWholeAreaList
    zmaxTab['completeness_opt'] = optCompList
    zmaxTab.write(zmaxFilePath, overwrite=True)

    return xmatchTab, zmaxTab

#------------------------------------------------------------------------------------------------------------
def _xmatchDirHasBestMatches(xmatchDir):
    """True if the cross-match folder exists and contains a best-matches table (xmatchtable_bestmatches_*.fits)."""
    return os.path.isdir(xmatchDir) and len(glob.glob(os.path.join(glob.escape(xmatchDir), 'xmatchtable_bestmatches_*.fits'))) > 0

#------------------------------------------------------------------------------------------------------------
def builddr(dataRelease='EDR', includeBands=None, includeQuality=None, removeDuplicates=True):
    """Build Data Release...

    Args:
        dataRelease (str): The name of the data release to build, default is 'EDR'.
        includeBands (list of str): List of bands to include, default is ['L'].
        includeQuality (list of int): List of quality flags to include, default is [0, 1].
        removeDuplicates (bool): Whether to remove duplicate entries, default is True.

    """

    if includeBands is None:
        includeBands = ['L']
    if includeQuality is None:
        includeQuality = [0, 1]

    print("\n" + "═" * 40)
    print("║ Building %s ║" %dataRelease)
    print("═" * 40 + "\n")

    parentImagesFile = os.path.join(startup.config['productsDir'], 'images.fits')

    parentImagesTab = atpy.Table().read(parentImagesFile)

    DRDir = os.path.join(startup.config['rootDir'], dataRelease)
    DRImageDir = os.path.join(DRDir, 'images')
    DRCatDir = os.path.join(DRDir, 'catalogs')
    DRRmsDir = os.path.join(DRDir, 'rms')
    DRXmatchDir = os.path.join(DRDir, 'xmatches_DECaLSDR10')

    os.makedirs(DRDir, exist_ok = True)
    os.makedirs(DRImageDir, exist_ok = True)
    os.makedirs(DRCatDir, exist_ok = True)
    os.makedirs(DRRmsDir, exist_ok = True)
    os.makedirs(DRXmatchDir, exist_ok = True)

    filteredRows = []

    for band in includeBands:
        for q in includeQuality:
            rows = parentImagesTab[(parentImagesTab['band'] == band) &
                                (parentImagesTab['quality'] == q)]
            filteredRows.append(rows)

    if filteredRows:
        filteredTable = atpy.vstack(filteredRows)
    else:
        filteredTable = atpy.Table()  # empty table

    nFiltered = len(filteredTable)

    # Remove duplicate sky positions by retaining the best observation per location:
    # prioritise lower quality flag, then lower RMS, then higher source count;
    # if all equal, select one randomly.

    if removeDuplicates is True:
        uniqueIdx = []
        raFilt = filteredTable["centre_RADeg"]
        decFilt = filteredTable["centre_decDeg"]
        coordsFilt = SkyCoord(ra=raFilt*u.deg, dec=decFilt*u.deg)
        for i in range(nFiltered):
            candidateIdx = i
            coordi = coordsFilt[i]
            matchedIndices = []
            for idx in uniqueIdx:
                sep = coordi.separation(coordsFilt[idx])
                if sep.arcsec < 1e-3:
                    matchedIndices.append(idx)
            if not matchedIndices:
                uniqueIdx.append(candidateIdx)
            else:
                bestIdx = candidateIdx
                for idx in matchedIndices:
                    qi = filteredTable[bestIdx]['quality']
                    qj = filteredTable[idx]['quality']
                    if qi < qj:
                        pass
                    elif qi > qj:
                        bestIdx = idx
                    else:
                        rmsi = filteredTable[bestIdx]['RMS_uJy/beam']
                        rmsj = filteredTable[idx]['RMS_uJy/beam']
                        if rmsi < rmsj:
                            pass
                        elif rmsi > rmsj:
                            bestIdx = idx
                        else:
                            nRadioSourcesi = filteredTable[bestIdx]['nRadioSources']
                            nRadioSourcesj = filteredTable[idx]['nRadioSources']
                            if nRadioSourcesi > nRadioSourcesj:
                                pass
                            elif nRadioSourcesi < nRadioSourcesj:
                                bestIdx = idx
                            else:
                                bestIdx = random.choice([bestIdx, idx])
                uniqueIdx = [idx for idx in uniqueIdx if idx not in matchedIndices]
                if bestIdx not in uniqueIdx:
                    uniqueIdx.append(bestIdx)
        DRImages = filteredTable[uniqueIdx]
    else:
        DRImages = filteredTable

    nDR = len(DRImages)

    print("\n%s has %d images.\n" % (dataRelease, nDR))

    for row in DRImages:

        fitsFile = os.path.join(startup.config['productsDir'], row['path'])
        pngFile = fitsFile.replace('.fits','.png')


        pybdsfCommonName = _normaliseCatName(row['radioCatPath'])

        catFile = os.path.join(startup.config['productsDir'], 'catalogs', pybdsfCommonName+'_bdsfcat.fits')
        rmsFile = os.path.join(startup.config['productsDir'], 'rms', pybdsfCommonName+'_rms.fits')
        rmsHistFile = os.path.join(startup.config['productsDir'], 'rms', pybdsfCommonName+'_rms_rmshist.txt')

        xmatchDir = os.path.join(startup.config['productsDir'],
                                 'xmatches_DECaLSDR10',
                                 'xmatch_%s_bdsfcat_DECaLSDR10_rband' %pybdsfCommonName)

        if os.path.exists(fitsFile) and not os.path.exists(os.path.join(DRDir, row['path'])):
            fitsFileRelPath = os.path.relpath(fitsFile, start=DRImageDir)
            pngFileRelPath = os.path.relpath(pngFile, start=DRImageDir)
            os.system("ln -s %s %s" %(fitsFileRelPath, DRImageDir))
            os.system("ln -s %s %s" %(pngFileRelPath, DRImageDir))
        if os.path.exists(catFile) and not os.path.exists(os.path.join(DRDir, row['radioCatPath'].replace('_srl',''))):
            catFileRelPath = os.path.relpath(catFile, start=DRCatDir)
            os.system("ln -s %s %s" %(catFileRelPath, DRCatDir))
        if os.path.exists(rmsFile) and not os.path.exists(os.path.join(DRRmsDir, os.path.basename(rmsFile))):
            rmsFileRelPath = os.path.relpath(rmsFile, start=DRRmsDir)
            rmsHistFileRelPath = os.path.relpath(rmsHistFile, start=DRRmsDir)
            os.system("ln -s %s %s" %(rmsFileRelPath, DRRmsDir))
            os.system("ln -s %s %s" %(rmsHistFileRelPath, DRRmsDir))
        # (only cross-match folders that have a best-matches table are linked into the data release)
        if _xmatchDirHasBestMatches(xmatchDir) and not os.path.exists(os.path.join(DRXmatchDir, os.path.basename(xmatchDir))):
            xmatchDirRelPath = os.path.relpath(xmatchDir, start=DRXmatchDir)
            os.system("ln -s %s %s" %(xmatchDirRelPath, DRXmatchDir))

    DRImages.sort('path')
    DRImagesOutFileName = DRDir+os.path.sep+'images_%s.fits' %dataRelease
    DRImages.write(DRImagesOutFileName, overwrite = True)
    print("\nWrote %s" % (DRImagesOutFileName))

    globalTabsDict = {band: None for band in includeBands}

    completenessDone = False
    completenessTypeCounts = {0: 0, 1: 0, 2: 0}
    tabFilesList=sorted(glob.glob(DRCatDir+os.path.sep+'*_bdsfcat.fits'))
    for t in tabFilesList:

        if t.find("srl_bdsfcat") == -1:
            tab=atpy.Table().read(t)
            freqGHz=tab.meta['FREQ0']/1e9
            bandKey=getBandKey(freqGHz)

            if bandKey not in includeBands:
                continue

            tab = tab[tab['Total_flux'] > 0.0]
            tab['freqGHz']=freqGHz
            tab['radCatPath'] = os.path.join("catalogs", os.path.basename(t))

            # Renaming source names from MSS to MWM
            tab['Source_name'] = [name.replace('MSS', 'MWM', 1) if name.startswith('MSS') else name for name in tab['Source_name']]
            # Creating resolved flag
            bmajCatDeg = float(tab.meta['BMAJ']) # in degrees
            bminCatDeg = float(tab.meta['BMIN']) # in degrees
            tab['beamMajDeg'] = bmajCatDeg
            tab['beamMinDeg'] = bminCatDeg
            resolved = (tab['Maj'] - bmajCatDeg) >= 2 * tab['E_Maj'] # Eq. 6 Heywood et al. 2022
            tab['resolvedFlag_axes'] = resolved.astype(int)
            # Size calculation
            tab['sizeMajAsec'] = np.sqrt(np.clip(tab['Maj']**2 - tab['beamMajDeg']**2, 0, None)) * 3600
            tab['sizeMinAsec'] = np.sqrt(np.clip(tab['Min']**2 - tab['beamMinDeg']**2, 0, None)) * 3600

            # assigning completeness
            tabFlux = tab['Total_flux'] # in Jy
            tabBaseName = os.path.basename(t).split('_bdsfcat.fits')[0]
            # (flux density and size if the 2D completeness is available, otherwise flux density only)
            completenessDir = _findSourceInjectionDir(os.path.join(DRDir, 'sinjection'), tabBaseName)
            sizeMajArcsec = _getSizeArcsec(tab, [('DC_Maj', 3600.0), ('sizeMajAsec', 1.0)])
            completenessSources, completenessType = assignRadioCompleteness(tabFlux, sizeMajArcsec,
                                                                            completenessDir, tabBaseName)
            tab['completeness_rad'] = completenessSources
            tab['completenessType_rad'] = np.int16(completenessType)   # 2: flux + size, 1: flux only, 0: not available
            completenessTypeCounts[completenessType] += 1
            if completenessType > 0:
                completenessDone = True
            else:
                print("\nCompleteness information is not available for %s." %(t))

            tab.meta.clear()
            if globalTabsDict[bandKey] is None:
                globalTabsDict[bandKey]=tab
            else:
                globalTabsDict[bandKey]=atpy.vstack([globalTabsDict[bandKey], tab])

    print("\nRadio completeness assigned using flux density and size for %d catalogues, flux density only for %d,"
          " and not available for %d." %(completenessTypeCounts[2], completenessTypeCounts[1], completenessTypeCounts[0]))
    if completenessDone is False:
        print("\nCompleteness information is not available for some catalogues in this data release."
              " Consider running the source injection task to compute completeness,"
               " and rerun builddr.")

    for bandKey in globalTabsDict.keys():
        if globalTabsDict[bandKey] is not None:
            outFileName=os.path.join(DRDir, "survey_catalog_%s_%s.fits" % (bandKey, dataRelease))
            globalTabsDict[bandKey].sort('DEC')
            globalTabsDict[bandKey].sort('RA')
            globalTabsDict[bandKey].meta['DR']=dataRelease
            globalTabsDict[bandKey].meta['BAND']=bandKey
            globalTabsDict[bandKey].meta['QTYVALS'] = includeQuality
            globalTabsDict[bandKey].meta['BERKVER']=__version__
            globalTabsDict[bandKey].meta['DATEMADE']=datetime.date.today().isoformat()
            globalTabsDict[bandKey].write(outFileName, overwrite = True)

            # Removing source-level duplicates (due to overlapping pointings)
            DRGlobalCatUnique = catalogs.removeDuplicateSources(globalTabsDict[bandKey], matchRadius_arcsec=6.0, raCol='RA', decCol='DEC',
                                                      fluxCol='Total_flux', fluxErrCol='E_Total_flux')
            DRGlobalCatUniqueFile = os.path.join(DRDir, "survey_catalog_%s_unique_%s.fits" % (bandKey, dataRelease))
            DRGlobalCatUnique.write(DRGlobalCatUniqueFile, overwrite=True)

            catalogs.catalog2DS9(globalTabsDict[bandKey], outFileName.replace(".fits", ".reg"),
                                 idKeyToUse = 'Source_name', RAKeyToUse = 'RA', decKeyToUse = 'DEC')
            print("\nWrote %s" % (outFileName))

    # Building cross-matched table and Luminosity table

    # cosmology for distance calculations
    cosmologyDR = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

    DRXmatchTab = None
    DRzmaxTab = None
    xmatchFilesList = sorted(glob.glob(DRXmatchDir+os.path.sep+'xmatch_*'+os.path.sep+'xmatchtable_bestmatches_*.fits'))

    filesToProcess = xmatchFilesList
    radCatNames = [os.path.basename(xmatchFileNow).split('xmatchtable_bestmatches_')[1].split('_DECaLSDR10_rband')[0]+'.fits'
                  for xmatchFileNow in filesToProcess]
    args = [(f, DRImages, DRDir, cosmologyDR, 'zmax', 'zphot_ifnot_spec_opt', radCatFileName) for f, radCatFileName in zip(filesToProcess, radCatNames)]

    doParallel = True

    if doParallel is True:
        nProcess = max(1, int(os.cpu_count()-1))
        with Pool(processes=nProcess) as pool:
            results = pool.map(processVmaxForXmatchFile, args)
    else:
        results = []
        for arg in args:
            result = processVmaxForXmatchFile(arg)
            results.append(result)

    DRXmatchTab = None
    DRzmaxTab = None
    for xmatchTab, zmaxTab in results:
        if 'ref_cat_opt' in zmaxTab.colnames:
            zmaxTab.remove_column('ref_cat_opt')
        if 'ref_cat_opt' in xmatchTab.colnames:
            xmatchTab.remove_column('ref_cat_opt')
        if xmatchTab is None:
            continue
        DRXmatchTab = xmatchTab if DRXmatchTab is None else atpy.vstack([DRXmatchTab, xmatchTab])
        DRzmaxTab   = zmaxTab if DRzmaxTab is None else atpy.vstack([DRzmaxTab, zmaxTab])

    if DRXmatchTab is not None:
        DRXmatchFile = os.path.join(DRDir, 'xmatchCat_DECaLSDR10_r_%s.fits' % dataRelease)
        DRXmatchTab.write(DRXmatchFile, overwrite=True)
        print("\nWrote %s" % (DRXmatchFile))

        # Unique version (duplicates from overlapping pointings removed)
        DRXmatchTabUnique = catalogs.removeDuplicateSources(DRXmatchTab, matchRadius_arcsec=6.0, raCol='RA_rad', decCol='DEC_rad',
                                                            fluxCol='Total_flux_rad', fluxErrCol='E_Total_flux_rad')
        DRXmatchUniqueFile = os.path.join(DRDir, 'xmatchCatUnique_DECaLSDR10_r_%s.fits' % dataRelease)
        DRXmatchTabUnique.write(DRXmatchUniqueFile, overwrite=True)
        print("\nWrote %s" % (DRXmatchUniqueFile))

    if DRzmaxTab is not None:
        DRzmaxFile = os.path.join(DRDir, 'zmaxCat_%s.fits' % dataRelease)
        DRzmaxTab.write(DRzmaxFile, overwrite=True)
        print("\nWrote %s" % (DRzmaxFile))

    # Removing source-level duplicates (due to overlapping pointings)

    DRzmaxTabUnique = catalogs.removeDuplicateSources(DRzmaxTab, matchRadius_arcsec=6.0, raCol='RA_rad', decCol='DEC_rad',
                                                      fluxCol='Total_flux_rad', fluxErrCol='E_Total_flux_rad')

    DRzmaxUniqueFile = os.path.join(DRDir, 'zmaxCatUnique_%s.fits' % dataRelease)
    DRzmaxTabUnique.write(DRzmaxUniqueFile, overwrite=True)

    # Removing pointings that are not well covered by opticalsurvey.
    # Visually checking the sky coverage of cross-matches using the RadOptSyPlot*.png in each crossmatch dir.
    # The pointings to avoid is in the directory xmatches_DECaLSDR10/xmatchPlotsToAvoid/

    DRzmaxTabUniqueBest = DRzmaxTabUnique.copy()
    xmatchPlotsToAvoidDir = os.path.join(DRXmatchDir, 'xmatchPlotsToAvoid')
    if os.path.exists(xmatchPlotsToAvoidDir) is True:
        plotNamesInAvoidDir = sorted(glob.glob(DRXmatchDir+os.path.sep+'xmatchPlotsToAvoid'+os.path.sep+'RadOptSkyPlot_*.png'))
        if len(plotNamesInAvoidDir) > 0:
            print("\nRemoving pointings that are not well covered by optical survey...\n")
            for plotName in plotNamesInAvoidDir:
                catAvoidName = os.path.basename(plotName).split('RadOptSkyPlot_')[1].replace('_DECaLSDR10_rband.png', '.fits')
                keep = np.array([_normaliseCatName(c) != _normaliseCatName(catAvoidName) for c in DRzmaxTabUniqueBest['radCatName']])
                DRzmaxTabUniqueBest = DRzmaxTabUniqueBest[keep]
                print("Removed %s from zmax table." %catAvoidName)

    DRzmaxUniqueBestFile = os.path.join(DRDir, 'zmaxCatUniqueBest_%s.fits' % dataRelease)
    DRzmaxTabUniqueBest.write(DRzmaxUniqueBestFile, overwrite=True)

    print("\n" + "═" * 40)
    print("Successfully built %s ║" %dataRelease)
    print("═" * 40 + "\n")

#------------------------------------------------------------------------------------------------------------
def builddb():
    """Build database...

    """

    # Build global catalog in each band (L-band, UHF)
    globalTabsDict={'L': None, 'UHF': None, 'S': None}

    # File to list the catalogues with wrapping issue
    catWrapIssueList = startup.config['productsDir']+os.path.sep+'catalogs_wrappingissue.txt'
    if os.path.exists(catWrapIssueList):
        os.remove(catWrapIssueList)

    tabFilesList=sorted(glob.glob(startup.config['productsDir']+os.path.sep+'catalogs'+os.path.sep+'*_bdsfcat.fits'))
    for t in tabFilesList:
        if t.find("srl_bdsfcat") == -1:
            tab=atpy.Table().read(t)
            # Fixing RA
            tabWrapAngle = catalogs.detectWrapAngle(tab['RA'])
            tab = catalogs.fixRA(tab, raCol='RA', wrapAngle=tabWrapAngle)
            catalogs.listCatalogInFile(t, catWrapIssueList)
            # if any(tab['RA'] < 0.0):
            #     # This pybdsf catalog has -180 to 180 wrapping. Need to change to 360 wrapping
            #     tab = catalogs.fixRA(tab, raCol='RA', wrapAngle=360)
            #     catalogs.listCatalogInFile(t, catWrapIssueList)
            tab = tab[tab['Total_flux'] > 0.0] # ignoring negative flux entries
            freqGHz=tab.meta['FREQ0']/1e9
            tab['freqGHz']=freqGHz # for scaling of the global catalog
            bandKey=getBandKey(freqGHz)

            # Renaming source names from MSS to MWM
            tab['Source_name'] = [name.replace('MSS', 'MWM', 1) if name.startswith('MSS') else name for name in tab['Source_name']]
            # Creating resolved flag
            bmajCatDeg = float(tab.meta['BMAJ']) # in degrees
            bminCatDeg = float(tab.meta['BMIN']) # in degrees
            tab['beamMajDeg'] = bmajCatDeg
            tab['beamMinDeg'] = bminCatDeg
            resolved = (tab['Maj'] - float(bmajCatDeg)) >= 2 * tab['E_Maj'] # Eq. 6 Heywood et al. 2022
            tab['resolvedFlag_axes'] = resolved.astype(int)

            #tab.meta=None # It'd be good to clear this... but the catalog matching stuff wants many things from here
            tab.meta.clear() #
            if globalTabsDict[bandKey] is None:
                globalTabsDict[bandKey]=tab
            else:
                globalTabsDict[bandKey]=atpy.vstack([globalTabsDict[bandKey], tab])

    # To implement: remove duplicates [we may want to make time series?]

    # Sorting, additional meta data, output
    for bandKey in globalTabsDict.keys():
        if globalTabsDict[bandKey] is not None:
            outFileName=startup.config['productsDir']+os.path.sep+"survey_catalog_%s.fits" % (bandKey)
            globalTabsDict[bandKey].sort('DEC')
            globalTabsDict[bandKey].sort('RA')
            globalTabsDict[bandKey].meta['BAND']=bandKey
            globalTabsDict[bandKey].meta['BERKVER']=__version__
            globalTabsDict[bandKey].meta['DATEMADE']=datetime.date.today().isoformat()
            globalTabsDict[bandKey].write(outFileName, overwrite = True)
            catalogs.catalog2DS9(globalTabsDict[bandKey], outFileName.replace(".fits", ".reg"),
                                 idKeyToUse = 'Source_name', RAKeyToUse = 'RA', decKeyToUse = 'DEC')
            print("\nWrote %s" % (outFileName))

    # Quality flags table
    qualFileName=startup.config['productsDir']+os.path.sep+"qualityFlags.csv"
    if os.path.exists(qualFileName):
        qualTab=atpy.Table().read(qualFileName)
    else:
        qualTab=None

    # Make image table - centre coords, radius [approx.], RMS, band, image path - UHF and L together.
    # Report command (when we make it) could load and dump some of that info

    outFileName=startup.config['productsDir']+os.path.sep+"images.fits"
    imageOutDir = startup.config['productsDir']+os.path.sep+"images"

    if os.path.exists(outFileName):
        existingImgTab = atpy.Table.read(outFileName)
        processedImages = set(existingImgTab['path'])
    else:
        existingImgTab = None
        processedImages = set()

    imgFilesList=sorted(glob.glob(startup.config['productsDir']+os.path.sep+"images"+os.path.sep+"pbcorr_*.fits"))

    # Removing images from images.fits that no longer exist in the productsDir
    if existingImgTab is not None:
        currentImages = set(os.path.join("images", os.path.basename(p)) for p in imgFilesList)
        mask = [p in currentImages for p in existingImgTab['path']]
        existingImgTab = existingImgTab[mask]
        processedImages = set(existingImgTab['path'])

    statsDictList=[]

    for imgFile in imgFilesList:
        pathName = imgFile.replace(startup.config['productsDir']+os.path.sep, '')

        imgFileName = pathName.split(os.path.sep)[-1]
        commonFilenamePart = imgFileName.split('.ms_pcalmask-MFS-image.fits')[0]
        commonFilenamePartForCat = commonFilenamePart.split(".")[0]
        radCatName = "catalogs"+os.path.sep+commonFilenamePartForCat+"_bdsfcat.fits"
        radCatPath = startup.config['productsDir']+os.path.sep+radCatName

        if not os.path.exists(radCatPath):
            nPybdsfSources = -99
            radCatNameToImages = '-'
        else:
            radCatNameToImages = radCatName
            radCatTab = atpy.Table.read(radCatPath)
            nPybdsfSources = len(radCatTab)

        if pathName in processedImages:
            continue #avoiding recomputing statistics

        statDict=images.getImagesStats(imgFile)
        if statDict is None:
            continue
        captureBlockId=os.path.split(statDict['path'])[-1].split('img_')[-1].split('_sdp')[0]
        statDict['captureBlockId']=captureBlockId
        statDict['path']=pathName
        statDict['band']=getBandKey(statDict['freqGHz'])
        statDict['radioCatPath']=radCatNameToImages
        statDict['nRadioSources']=nPybdsfSources
        statsDictList.append(statDict)

    if statsDictList:
        newImgTab = atpy.Table()
        for key in statsDictList[0].keys():
            newImgTab[key] = [s[key] for s in statsDictList]
    else:
        newImgTab = None

    # Combine old and new tables if both exist
    if existingImgTab is not None and newImgTab is not None:
        imgTab = existingImgTab.copy()
        imgTab = atpy.vstack([imgTab, newImgTab])
    elif existingImgTab is not None:
        imgTab = existingImgTab
    else:
        imgTab = newImgTab

    # plotting images
    for row in imgTab:
        imgPath = startup.config['productsDir'] + os.path.sep + row['path']
        plotName = os.path.basename(row['path']).replace('.fits', '.png')
        plotPath = os.path.join(imageOutDir, plotName)
        if not os.path.exists(plotPath):
            statDictRow = dict(zip(imgTab.colnames, row))
            images.plotImages(imgPath, imageOutDir, axLabelDeg = True, statsDict=statDictRow, overwrite=False)

    # Update quality flag column - we only use quality column from qualTab
    imgTab['quality']=99
    imgTab['comments']=Column(length=len(imgTab), dtype='U200')
    if qualTab is not None:
        for irow in imgTab:
            mask=irow['path'] == qualTab['path']
            if 'quality' in qualTab.keys():
               #pulling already existing quality and assigning 99 to newly added images
                if mask.any():
                    irow['quality']=qualTab[mask]['quality'][0]
                else:
                    irow['quality'] = 99
            if 'comments' in qualTab.keys():
                if mask.any():
                    irow['comments']=qualTab[mask]['comments'][0]
                else:
                    irow['comments'] = 'NA'

    # Output

    imgTab.sort('path')
    imgTab.meta['BERKVER']=__version__
    imgTab.meta['DATEMADE']=datetime.date.today().isoformat()
    imgTab.write(qualFileName, overwrite = True)
    print("\nWrote %s" % (qualFileName))
    imgTab.remove_column('comments')
    imgTab.write(outFileName, overwrite = True)
    print("\nWrote %s" % (outFileName))

    # Generate survey mask in some format - we'll use that to get total survey area

#------------------------------------------------------------------------------------------------------------
def _processXmatchOnRadioCatalogue(args):

    radCat, optSurvey, optBandToMatch, optPosErrAsecValue = args

    catalogName = radCat.split(os.path.sep)[-1]

    # making a subscript for this particular match
    outSubscript = '%s_%s_%sband' %(catalogName.replace('.fits',''), optSurvey, optBandToMatch)# , str(np.round(searchRadiusArcsec, 1)).replace(".","p"))

    xmatchDirPath = os.path.join(startup.config['productsDir'], 'xmatches_%s' %optSurvey)
    os.makedirs(xmatchDirPath, exist_ok = True)

    radCatTab = atpy.Table().read(radCat)
    freqGHz=radCatTab.meta['FREQ0']/1e9
    radBandName=getBandKey(freqGHz)

    # determining search radius based on radio positional uncertainty and optical positional uncertainty

    sigmaRadPosDeg = np.sqrt(radCatTab['E_RA']**2 + radCatTab['E_DEC']**2)
    sigmaRadPosMeanDeg = np.mean(sigmaRadPosDeg)

    # add optical positional uncertainty (e.g. 0.2 arcsec for DECaLS)
    sigmaOptDeg = optPosErrAsecValue/3600.

    # combined sigma
    sigmaTotalDeg = np.sqrt(sigmaRadPosMeanDeg**2 + sigmaOptDeg**2)

    # search radius as 3-sigma
    searchRadiusDegVal = 3 * sigmaTotalDeg
    searchRadiusArcsec = searchRadiusDegVal * 3600

    xmatchTab = crossmatch.xmatchRadioOptical(radioCatFilePath=radCat,
                                                radioBand=radBandName,
                                                xmatchDirPath=xmatchDirPath,
                                                optSurvey=optSurvey,
                                                optMagCol = optBandToMatch,
                                                searchRadiusArcsec = searchRadiusArcsec,
                                                radRACol='RA', radERACol='E_RA',
                                                radDecCol='DEC', radEDecCol='E_DEC',
                                                radEMajCol='E_Maj',
                                                radEMinCol='E_Min', radPACol='PA',
                                                outSubscript=outSubscript,
                                                optPosErrAsecValue=optPosErrAsecValue, nMagBins=15, beamSizeArcsecValue=6.0,
                                                skipIfExists=True
                                                )

    if xmatchTab is None:
        return None

    if 'radCatPath' not in xmatchTab.columns:
        radCatPathName = radCat.replace(startup.config['productsDir']+os.path.sep, '')
        xmatchTab.add_column(radCatPathName, name='radCatPath', index=0)

    return xmatchTab

#------------------------------------------------------------------------------------------------------------
def xmatch(optSurveyInput='decalsdr10'):
    """Does cross-matching...

    """

    optSurveyParamDict = {'decalsdr10': {'surveyTag': 'DECaLSDR10',
                                         'optPosErrAsecValue': 0.2},
                          'rubindp1': {'surveyTag': 'RubinDP1',
                                         'optPosErrAsecValue': 0.05}}

    optSurvey = optSurveyParamDict[optSurveyInput]['surveyTag']
    optPosErrAsecValue = optSurveyParamDict[optSurveyInput]['optPosErrAsecValue']
    optBandToMatch = 'r'
    #searchRadiusArcsec = 4.0

    print("\n" + "═" * 40)
    print("║ Cross-matching with %s ║" %optSurvey)
    print("═" * 40 + "\n")

    globalBestXmatchTabName=startup.config['productsDir']+os.path.sep+"xmatchCat_%s_%sband.fits" %(optSurvey, optBandToMatch)

    globalBestXmatchTab=None
    radCatFilesList=sorted([f for f in glob.glob(startup.config['productsDir']+os.path.sep+'catalogs'+os.path.sep+'*_bdsfcat.fits')
                            if 'srl_bdsfcat' not in os.path.basename(f)])

    args = [(radCat, optSurvey, optBandToMatch, optPosErrAsecValue) for radCat in radCatFilesList]

    doParallelXmatch = False
    if doParallelXmatch is True:
        nProcess = max(1, int(os.cpu_count()-4))
        with Pool(processes=nProcess) as pool:
            tables = pool.map(_processXmatchOnRadioCatalogue, args)
        tables = [t for t in tables if t is not None]
    else:
        tables = []
        for arg in args:
            table = _processXmatchOnRadioCatalogue(arg)
            if table is not None:
                tables.append(table)

    if len(tables) == 0:
        print("\nNo radio catalogs had successful cross-matching with %s!\n" % optSurvey)
        return

    for tab in tables:
        if tab is not None and 'ref_cat_opt' in tab.colnames:
            tab.remove_column('ref_cat_opt')

    globalBestXmatchTab = atpy.vstack(tables)

    globalBestXmatchTab.write(globalBestXmatchTabName, overwrite=True)

    print("\n" + "-" * 100)
    print("\nWrote %s" % (globalBestXmatchTabName))
    print("\n" + "-" * 100)

    # Unique version: a radio source detected in overlapping pointings is kept once (the detection with the higher
    # signal-to-noise ratio of its flux density), as for the survey catalogue
    globalBestXmatchUniqueTabName = startup.config['productsDir']+os.path.sep+"xmatchCatUnique_%s_%sband.fits" %(optSurvey, optBandToMatch)
    globalBestXmatchUniqueTab = catalogs.removeDuplicateSources(globalBestXmatchTab, matchRadius_arcsec=6.0,
                                                                raCol='RA_rad', decCol='DEC_rad',
                                                                fluxCol='Total_flux_rad', fluxErrCol='E_Total_flux_rad')
    globalBestXmatchUniqueTab.write(globalBestXmatchUniqueTabName, overwrite=True)

    print("\nWrote %s" % (globalBestXmatchUniqueTabName))
    print("\n" + "-" * 100)

#------------------------------------------------------------------------------------------------------------
def collect():
    """Collect...

    """

    print("Collecting processed data products...")
    if 'BERK_NODES_FILE' not in os.environ.keys():
        print("You need to set the BERK_NODES_FILE environment variable to use the 'collect' task.")
        sys.exit()

    try:
        stubs=[]
        with open(os.environ['BERK_NODES_FILE'], "r") as inFile:
            for line in inFile.readlines():
                stubs.append(line.strip())
    except:
        import urllib.request  # the lib that handles the url stuff
        stubs=[]
        for line in urllib.request.urlopen(os.environ['BERK_NODES_FILE']):
            l=line.decode('utf-8')
            if l[0] != "#" and len(l) > 3:
                stubs.append(l.strip())

    # Get images
    print("Collecting images...")
    toPath=startup.config['productsDir']+os.path.sep+"images"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        imgPath="processing/*/IMAGES/pbcorr*pcalmask-MFS-image.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+imgPath, toPath)
        os.system(cmd)

    # Get catalogs
    print("Collecting catalogs...")
    toPath=startup.config['productsDir']+os.path.sep+"catalogs"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        catPath="processing/*/IMAGES/pbcorr_trim_*_pybdsf/*_bdsfcat.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+catPath, toPath)
        os.system(cmd)

    # Get rms maps
    print("Collecting rms images...")
    toPath=startup.config['productsDir']+os.path.sep+"rms"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        rmsPath="processing/*/IMAGES/pbcorr_trim_*_pybdsf/*_rms.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+rmsPath, toPath)
        os.system(cmd)

    # Get mean maps
    print("Collecting mean images...")
    toPath=startup.config['productsDir']+os.path.sep+"mean"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        meanPath="processing/*/IMAGES/pbcorr_trim_*_pybdsf/*_mean.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+meanPath, toPath)
        os.system(cmd)

    # Get residual maps
    print("Collecting residual images...")
    toPath=startup.config['productsDir']+os.path.sep+"residual"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        residualPath="processing/*/IMAGES/pbcorr_trim_*_pybdsf/*_gaus_resid.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+residualPath, toPath)
        os.system(cmd)

    print("Finished!")
    sys.exit()


#------------------------------------------------------------------------------------------------------------
def setup(captureBlockId):
    """Set up the data for further processing with process1 and process2 tasks

    """

    # Forget staging dir, just do a symbolic link to the MSCache dir
    MSPath=os.environ['BERK_MSCACHE']+os.path.sep+captureBlockId+"_sdp_l0.ms"

    # Setup in processing dir
    MSProcessDir=startup.config['processingDir']+os.path.sep+captureBlockId
    if os.path.exists(MSProcessDir+os.path.sep+'project_info.json') is True:
        print("Observation already setup - re-run with process1 task to continue")
        sys.exit()
    os.makedirs(MSProcessDir)
    os.chdir(MSProcessDir)
    os.system("ln -s %s" % (os.path.abspath(MSPath)))
    oxdirs=['setups', 'tools', 'oxkat', 'data']
    for oxdir in oxdirs:
        os.system("ln -s %s" % (startup.config['oxkatDir']+os.path.sep+oxdir))

    # Run oxkat scripts
    # Setup
    cmd="python3 setups/0_GET_INFO.py %s" % (os.environ['BERK_PLATFORM'])
    os.system(cmd)
    if os.path.exists("submit_info_job.sh") is False:
        raise Exception("Failed to generate submit_info_job.sh")
    cmd="berk_chain %s submit_info_job.sh" % (startup.config['workloadManager'])
    jobID=jobs.submitJob(cmd, "SUBMIT_SETUP", dependentJobIDs = None,
                          workloadManager = startup.config['workloadManager'],
                          time = "00:05:00")
    print("Setup jobs submitted")
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def process1(captureBlockId):
    """1GC processing

    """

    MSPath=os.environ['BERK_MSCACHE']+os.path.sep+captureBlockId+"_sdp_l0.ms"
    MSProcessDir=startup.config['processingDir']+os.path.sep+captureBlockId
    os.chdir(MSProcessDir)

    if os.path.exists('project_info.json') is False:
        raise Exception("You need to run the 'setup' task on this observation first")

    # 1GC
    cmd="python3 setups/1GC.py %s" % (os.environ['BERK_PLATFORM'])
    os.system(cmd)
    if os.path.exists("submit_1GC_jobs.sh") is False:
        raise Exception("Failed to generate submit_1GC_jobs.sh")
    cmd="berk_chain %s submit_1GC_jobs.sh" % (startup.config['workloadManager'])
    jobID=jobs.submitJob(cmd, "SUBMIT_1GC", dependentJobIDs = None,
                          workloadManager = startup.config['workloadManager'],
                          time = "00:05:00")
    print("1GC jobs submitted")
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def process2(captureBlockId):
    """Imaging and 2GC processing

    """

    MSPath=os.environ['BERK_MSCACHE']+os.path.sep+captureBlockId+"_sdp_l0.ms"
    MSProcessDir=startup.config['processingDir']+os.path.sep+captureBlockId
    os.chdir(MSProcessDir)

    # FLAG and 2GC can be chained
    cmd="python3 setups/FLAG.py %s" % (os.environ['BERK_PLATFORM'])
    os.system(cmd)
    cmd="python3 setups/2GC.py %s" % (os.environ['BERK_PLATFORM'])
    os.system(cmd)
    if os.path.exists("submit_flag_jobs.sh") is False:
        raise Exception("Failed to generate submit_flag_jobs.sh")
    if os.path.exists("submit_2GC_jobs.sh") is False:
        raise Exception("Failed to generate submit_2GC_jobs.sh")
    cmd="berk_chain %s submit_flag_jobs.sh submit_2GC_jobs.sh" % (startup.config['workloadManager'])
    jobID=jobs.submitJob(cmd, "SUBMIT_2GC", dependentJobIDs = None,
                          workloadManager = startup.config['workloadManager'],
                          time = "00:05:00")
    print("FLAG and 2GC jobs submitted")
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def analyse(captureBlockId):
    """Analyse...

    """

    # Setup in processing dir
    MSProcessDir=startup.config['processingDir']+os.path.sep+captureBlockId
    if os.path.exists(MSProcessDir) == False:
        raise Exception("Processing directory %s does not exist - you need to process the data before the analyse task will run." % (MSProcessDir))
    os.chdir(MSProcessDir)
    os.system("ln -s %s" % (startup.config["catalogScriptsDir"]+os.path.sep+"sourcefinding.py"))
    #os.system("ln -s %s" % (startup.config["catalogScriptsDir"]+os.path.sep+"catalog_matching.py"))
    os.system("ln -s %s" % (startup.config["catalogScriptsDir"]+os.path.sep+"parsets"))

    # Source finding is fairly lightweight so we put everything in one job script
    # We will have issues with needing to see the internet to fetch cross match catalogs though
    # So we will need to cache NVSS catalogs for a given direction when doing 'stage'
    # OR forget using catalog_matching.py here and just do later with database/catalog scripts
    imgPaths=glob.glob("IMAGES/*pcalmask-MFS-image.fits")
    for i in imgPaths:
        if i.find("pbcorr_trim") != -1:
            continue
        imgPath=os.path.abspath(i)
        cmd="mkat_primary_beam_correct %s -T" % (imgPath)

        imgDir, imgFileName=os.path.split(imgPath)
        pbcorrImgPath=imgDir+os.path.sep+"pbcorr_trim_"+imgFileName
        cmd=cmd+"\npython3 sourcefinding.py c %s -o fits --survey MWM" % (pbcorrImgPath)

        # The bit below isn't going to work on compute nodes - so we may as well move this
        #label=imgFileName.split(".ms_")[0].split(".")[0]
        #catPath=imgDir+os.path.sep+"pbcorr_trim_"+label+"_pybdsf"+os.path.sep+"pbcorr_trim_"+label+"_bdsfcat.fits"
        #cmd=cmd+"\npython3 catalog_matching.py %s NVSS --astro --flux" % (catPath)

        jobID=jobs.submitJob(cmd, 'source-finding-%s' % (imgFileName), dependentJobIDs = None,
                              nodes = 1, tasks = 20, mem = 64000, time = "02:00:00",
                              cmdIsBatchScript = False,
                              workloadManager = startup.config['workloadManager'])
        print("Submitted source finding and analysis job %d" % (jobID))
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def summarize(dataBase='PARENT'):
    """Summarize the progress in berk processing.

    Prints a summary and makes plots (sky coverage, RMS area coverage, source counts). For the EDR it also writes the
    lists of pointings without or with only partial optical (DECaLS) coverage and the cross-match statistics (Q0,
    LR threshold, reliability, completeness) of every pointing, see _summarizeEDRXmatches.

    Args:
            dataBase (str): 'PARENT' for all the processed data, or 'EDR' for the early data release.

    """

    if dataBase == 'PARENT':
        dataBaseDir = startup.config['productsDir']
        subScript = ''
    elif dataBase == 'EDR':
        dataBaseDir = os.path.join(startup.config['rootDir'], dataBase)
        subScript = '_EDR'

    imagesFileName = dataBaseDir+os.path.sep+"images%s.fits" %subScript
    xmatchFileName = dataBaseDir+os.path.sep+"xmatchCat_DECaLSDR10_r%s.fits" %subScript
    rmsDirPath = dataBaseDir+os.path.sep+"rms"

    imagesTab = atpy.Table().read(imagesFileName)
    xmatchTab = atpy.Table().read(xmatchFileName)

    allBandColors = {'L': '#e35c1e', 'UHF': '#1e21e3', 'S': '#145a32'}
    bandsInDB = set(imagesTab['band'])
    bandColorDict = {band: color for band, color in allBandColors.items() if band in bandsInDB}
    orderedBands = list(bandColorDict.keys())
    bandTotalAreaDict = {}

    # Printing summary

    print("\n" + "═" * 50)
    print(("║ SUMMARY OF MEERKAT OBSERVATIONS IN %s DATABASE ║"%dataBase).center(50))
    print("═" * 50 + "\n")

    for band in orderedBands:
        bandMaskImages = imagesTab['band'] == band
        bandDataImages = imagesTab[bandMaskImages]
        bandCountImages = len(bandDataImages)
        bandTotalArea = bandDataImages['skyArea_sqDeg'].sum()
        bandTotalAreaDict[band] = bandTotalArea

        catFileName = os.path.join(dataBaseDir, "survey_catalog_%s%s.fits" % (band, subScript))
        catalogTab = atpy.Table().read(catFileName)
        nPybdsfSources = len(catalogTab)

        bandMaskXmatches = xmatchTab['band'] == band
        bandDataXmatches = xmatchTab[bandMaskXmatches]
        bandCountXmatches = len(bandDataXmatches)

        print("\n---------------- %s band --------------------" %band)
        print("\nNumber of pointings: %d \
                \nTotal area: %.3f sq. deg. \
                \nNumber of PyBDSF sources: %d \
              \nNumber of DECaLS cross-matches: %d\n" % (bandCountImages, bandTotalArea, nPybdsfSources, bandCountXmatches))

    # Plotting sky coverage

    skyPlotName = dataBaseDir+os.path.sep+'MeerKAT_pointings%s.png' %subScript

    summaryPlots.plotSkyCoverage(imagesTab, bandColorDict=bandColorDict, plotOutPath=skyPlotName, plotProjection='aitoff')


    # Plotting RMS area coverage

    # The RMS area calculations are slow. They are skipped if the text files written by the unique RMS area calculation
    # (MeerKAT_RMS_area_coverage_unique<DR>_<band>.txt, which the source counts read) already exist.
    # Delete them to calculate the RMS area coverage again.

    rmsAreaCoveragePlotName = dataBaseDir+os.path.sep+'MeerKAT_RMS_area_coverage%s.png' % subScript
    rmsAreaCoverageUniqPlotName = dataBaseDir+os.path.sep+'MeerKAT_RMS_area_coverage_unique%s.png' % subScript
    uniqTxtFiles = [rmsAreaCoverageUniqPlotName.replace('.png', '') + '_%s.txt' % band for band in orderedBands]

    if len(uniqTxtFiles) > 0 and all(os.path.exists(f) for f in uniqTxtFiles):
        print("\nRMS area coverage already calculated (found %s): skipping. Delete the file(s) to calculate it again."
              % ', '.join(os.path.basename(f) for f in uniqTxtFiles))
    else:
        summaryPlots.plotRMSAreaCoverageCumulative(rmsDirPath, rmsAreaCoveragePlotName, bandColorDict, nRMSBins=100)

        # Plotting RMS area coverage taking into account overlapping pointings

        # Finding overlapping pointings
        overlapGroupsByBand, isolatedIndicesByBand = images.findOverlappingPointings(imagesTab,  raCol='centre_RADeg', decCol='centre_decDeg', areaCol='skyArea_sqDeg')

        summaryPlots.plotRMSAreaCoverageCumulativeUnique(rmsDirPath, rmsAreaCoverageUniqPlotName, bandColorDict, imagesTab, overlapGroupsByBand, isolatedIndicesByBand, nRMSBins=100)

    # # Plotting sourcecount

    sourceCountPlotName = dataBaseDir+os.path.sep+'MeerKAT_sourcecount%s.png' % subScript
    summaryPlots.plotSourceCounts(dataBaseDir, catSubScript=subScript, fluxCol='Total_flux', fluxMin=1E-5, fluxMax=1.0, nFluxBins=39, bandColorDict=bandColorDict, plotOutPath=sourceCountPlotName)


    # Cross-match lists and statistics of the EDR pointings
    if dataBase == 'EDR':
        try:
            _summarizeEDRXmatches(dataBase, imagesTab)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print("\nCould not write the EDR cross-match lists and statistics: %s" % e)

#------------------------------------------------------------------------------------------------------------
def report():
    """Report...

    """

    # We'll tidy this all up later...

    # Get catalogs if we don't have them already...
    firstPath=startup.config['cacheDir']+os.path.sep+"first_14dec17.fits"
    if os.path.exists(firstPath) == False:
        print("Fetching and caching FIRST catalog [this is only done once]")
        topDir=os.getcwd()
        os.chdir(os.path.abspath(startup.config['cacheDir']))
        os.system("wget http://sundog.stsci.edu/first/catalogs/first_14dec17.fits.gz")
        os.system("gunzip first_14dec17.fits.gz")
        os.chdir(topDir)
    firstTab=atpy.Table().read(firstPath)

    # Compare L-band flux densities with FIRST [which only goes to dec -10 degrees]
    tab=atpy.Table().read(os.path.abspath(startup.config['productsDir'])+os.path.sep+"survey_catalog_L.fits")
    x_tab, x_firstTab, rDeg=catalogs.crossMatch(tab, firstTab, radiusArcmin = 6.0/60)

    # /home/matty/Documents/Rubin_x_MeerKAT/products/survey_catalog_L.fits


    # http://sundog.stsci.edu/first/catalogs/first_14dec17.fits.gz

    import IPython
    IPython.embed()
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def _getSizeArcsec(tab, candidates):
    """Returns the deconvolved major-axis size of each source in arcsec, or None if no usable column exists.

    Args:
        tab (astropy.table.Table): Table with the radio sources.
        candidates (list of tuple): (column name, factor that converts the column to arcsec), in order of preference.

    Returns:
        numpy.ndarray or None: Size in arcsec for each row of the table (NaN where undefined).
    """
    for colName, toArcsec in candidates:
        if colName in tab.colnames:
            sizes = np.ma.filled(np.ma.asarray(tab[colName], dtype=float), np.nan)
            return np.asarray(sizes, dtype=float) * toArcsec
    return None

#------------------------------------------------------------------------------------------------------------
def _loadCompleteness2D(path2D):
    """Loads a flux density x size completeness grid (written by sourceinjection) and builds an interpolator.

    Cells without enough injected sources are NaN in the file; they are filled with the value of the nearest
    valid cell. Queries outside the grid are clamped to the grid edge (as np.interp does for the 1D curve).

    Returns:
        dict with the interpolator and grid limits, or None if the grid has no valid cell.
    """
    from scipy.interpolate import RegularGridInterpolator
    from scipy.ndimage import distance_transform_edt

    with np.load(path2D) as grid:
        fluxEdges = np.asarray(grid['fluxEdgesJy'], dtype=float)
        sizeEdges = np.asarray(grid['sizeEdgesArcsec'], dtype=float)
        comp = np.array(grid['completeness'], dtype=float)
        spacing = str(grid['sizeSpacing']) if 'sizeSpacing' in grid.files else None

    if spacing is None:   # older files: infer the spacing of the size axis
        steps = np.diff(sizeEdges)
        spacing = 'linear' if np.allclose(steps, steps[0]) else 'log'
    logSize = (spacing == 'log')

    fluxCen = np.sqrt(fluxEdges[:-1] * fluxEdges[1:])
    sizeCen = np.sqrt(sizeEdges[:-1] * sizeEdges[1:]) if logSize else 0.5 * (sizeEdges[:-1] + sizeEdges[1:])

    valid = np.isfinite(comp)
    if not valid.any():
        return None
    if not valid.all():
        nearest = distance_transform_edt(~valid, return_distances=False, return_indices=True)
        comp = comp[tuple(nearest)]

    xFlux = np.log10(fluxCen)
    ySize = np.log10(sizeCen) if logSize else sizeCen
    interp = RegularGridInterpolator((xFlux, ySize), comp, method='linear', bounds_error=False, fill_value=None)

    return {'interp': interp, 'logSize': logSize,
            'xRange': (xFlux.min(), xFlux.max()), 'yRange': (ySize.min(), ySize.max())}

#------------------------------------------------------------------------------------------------------------
def _findSourceInjectionDir(sinjectRoot, baseName):
    """Finds the source injection folder of a radio catalogue, i.e. <sinjectRoot>/sinjected_<type>_<baseName>_<N>sources_<R>reps.

    Folders from injections with random sizes ('extended') are preferred over point-source-only ones
    ('unresolved'). Only folders that contain completeness files for this catalogue are considered. If there
    are several (e.g. runs with different numbers of sources or repetitions), the most recent one is used.

    Args:
        sinjectRoot (str): Directory with the source injection folders (e.g. <DRDir>/sinjection).
        baseName (str): Name of the radio catalogue without suffix.

    Returns:
        str or None: Path to the folder, or None if no suitable folder exists.
    """
    for sourceType in ('extended', 'unresolved'):
        pattern = os.path.join(glob.escape(sinjectRoot),
                               'sinjected_%s_%s_*sources_*reps' % (sourceType, glob.escape(baseName)))
        dirs = [d for d in glob.glob(pattern)
                if os.path.exists(os.path.join(d, 'completeness_%s.txt' % baseName))
                or os.path.exists(os.path.join(d, 'completeness2D_%s.npz' % baseName))]
        if len(dirs) > 0:
            dirs.sort(key=os.path.getmtime)
            if len(dirs) > 1:
                print("Several '%s' source injection folders found for %s; using the most recent: %s"
                      % (sourceType, baseName, os.path.basename(dirs[-1])))
            return dirs[-1]
    return None

#------------------------------------------------------------------------------------------------------------
def assignRadioCompleteness(fluxJy, sizeMajAsec, completenessDir, baseName):
    """Assigns the radio source detection completeness to each source.

    If the flux density x size completeness grid (completeness2D_<baseName>.npz) exists and the sizes are
    available, the completeness is interpolated in both flux density and deconvolved size. Otherwise the
    1D completeness as a function of flux density (completeness_<baseName>.txt) is used.

    Args:
        fluxJy (array): Flux density of each source in Jy.
        sizeMajAsec (array or None): Deconvolved major-axis FWHM of each source in arcsec (0 for unresolved
            sources). None if sizes are not available.
        completenessDir (str or None): Directory with the completeness files (None if there is none).
        baseName (str): Name of the radio catalogue without suffix (as used in the completeness file names).

    Returns:
        tuple: (completeness, completenessType). completeness is an array (-99 if no completeness information
            is available). completenessType is 2 for flux + size, 1 for flux only and 0 if not available.
    """

    fluxJy = np.asarray(fluxJy, dtype=float)
    if completenessDir is None:
        return np.full(len(fluxJy), -99.0), 0
    path1D = os.path.join(completenessDir, 'completeness_%s.txt' % baseName)
    path2D = os.path.join(completenessDir, 'completeness2D_%s.npz' % baseName)

    comp1D = None
    if os.path.exists(path1D):
        tab1D = atpy.Table.read(path1D, format='ascii')
        comp1D = np.interp(fluxJy, np.asarray(tab1D['FluxBinCentre_Jy'], dtype=float),
                           np.asarray(tab1D['MeanCompleteness'], dtype=float))

    if sizeMajAsec is not None and os.path.exists(path2D):
        grid = _loadCompleteness2D(path2D)
        if grid is not None:
            size = np.asarray(sizeMajAsec, dtype=float)
            sizeOK = np.isfinite(size)
            size = np.where(sizeOK, size, 0.0)

            x = np.log10(np.clip(fluxJy, 1e-30, None))
            y = np.log10(np.clip(size, 1e-3, None)) if grid['logSize'] else size   # size 0 -> smallest size bin
            x = np.clip(x, *grid['xRange'])
            y = np.clip(y, *grid['yRange'])

            comp2D = np.full(len(fluxJy), np.nan)
            ok = np.isfinite(fluxJy)
            comp2D[ok] = grid['interp'](np.column_stack([x[ok], y[ok]]))
            if comp1D is not None:   # sources without a size: flux only
                comp2D[~sizeOK] = comp1D[~sizeOK]
            return np.clip(comp2D, 0.0, 1.0), 2

    if comp1D is not None:
        return comp1D, 1

    return np.full(len(fluxJy), -99.0), 0

#------------------------------------------------------------------------------------------------------------
def _fitsImageIsComplete(path):
    """True if the FITS image exists and is not truncated, i.e. the file holds all the data its header declares."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    try:
        with fits.open(path, memmap=True) as hdul:
            hdu = hdul[0]
            hdr = hdu.header
            nBytes = abs(hdr['BITPIX']) // 8 * int(np.prod([hdr['NAXIS%d' % k] for k in range(1, hdr['NAXIS'] + 1)]))
            return os.path.getsize(path) >= hdu._data_offset + nBytes
    except Exception:
        return False

#------------------------------------------------------------------------------------------------------------
def _findEDRImage(baseName):
    """Finds the radio image of an EDR pointing in <rootDir>/EDR/images, given the name of its PyBDSF catalogue
    (without the _bdsfcat.fits suffix)."""
    imageDir = os.path.join(startup.config['rootDir'], 'EDR', 'images')
    matches = sorted(glob.glob(os.path.join(glob.escape(imageDir), glob.escape(baseName) + '.*image.fits')))
    if len(matches) == 0:
        matches = sorted(glob.glob(os.path.join(glob.escape(imageDir), glob.escape(baseName) + '*image.fits')))
    if len(matches) != 1:
        raise FileNotFoundError("Expected one image matching '%s*image.fits' in %s, found %d: %s"
                                % (baseName, imageDir, len(matches), [os.path.basename(m) for m in matches]))
    return matches[0]

#------------------------------------------------------------------------------------------------------------
def _ensurePyBDSFMaps(baseName, mapsDir):
    """Makes sure that complete residual, rms and mean maps exist for an image, next to its catalogue.

    The catalogue with the source names (<baseName>_bdsfcat.fits) and the three maps are taken from the products
    directory if they are all there and the maps are not truncated. Otherwise PyBDSF is run again on the image in
    <rootDir>/EDR/images (with the parameters in sourceinjection.pybdsfArgs) and the maps, together with the raw PyBDSF
    catalogue (<baseName>_srl_bdsfcat.fits, which has no source names and only serves to check that the extraction
    finished and to compare the number of sources), are written to mapsDir. sourceinjection.executeSingle looks in
    the current directory before the products directory, so estimateCompleteness runs with mapsDir as the current
    directory. The catalogue with the source names must exist in the products directory: it cannot be re-created here.

    Args:
        baseName (str): Name of the radio catalogue without the _bdsfcat.fits suffix.
        mapsDir (str): Directory for re-extracted files.
    """

    prodDir = startup.config['productsDir']
    namedCat = baseName + '_bdsfcat.fits'        # catalogue with source names (products directory)
    rawCat = baseName + '_srl_bdsfcat.fits'      # raw PyBDSF catalogue (written by the re-extraction)
    mapFiles = {'residual': baseName + '_gaus_resid.fits', 'rms': baseName + '_rms.fits', 'mean': baseName + '_mean.fits'}

    def badMaps(directories):
        """Names of the maps that are missing or truncated; directories maps each map type to a folder."""
        return [name for kind, name in mapFiles.items() if not _fitsImageIsComplete(os.path.join(directories[kind], name))]

    if not os.path.exists(os.path.join(prodDir, 'catalogs', namedCat)):
        raise FileNotFoundError("The catalogue with source names %s is missing in %s"
                                % (namedCat, os.path.join(prodDir, 'catalogs')))

    # maps re-extracted earlier (the raw catalogue is written last, so it shows that the extraction finished)
    badInMaps = badMaps({k: mapsDir for k in mapFiles}) + ([] if os.path.exists(os.path.join(mapsDir, rawCat)) else [rawCat])
    if len(badInMaps) == 0:
        print("Using the maps re-extracted earlier in %s" % mapsDir)
        return

    # leftovers of an interrupted extraction must not be picked up (files in the current directory are used first)
    for name in list(mapFiles.values()) + [rawCat]:
        if os.path.lexists(os.path.join(mapsDir, name)):
            os.remove(os.path.join(mapsDir, name))

    badInProducts = badMaps({k: os.path.join(prodDir, k) for k in mapFiles})
    if len(badInProducts) == 0:
        print("Complete catalogue and maps found in the products directory.")
        return

    print("Missing or truncated in the products directory: %s" % badInProducts)
    imagePath = _findEDRImage(baseName)
    linkPath = os.path.join(mapsDir, os.path.basename(imagePath))   # PyBDSF writes its log next to the image it is given
    if os.path.lexists(linkPath):
        os.remove(linkPath)
    os.symlink(os.path.realpath(imagePath), linkPath)

    print("Re-running PyBDSF on %s" % imagePath)
    sourceinjection.sourceExtract(imageFile=linkPath, pybdsfBaseName=baseName, saveRMSMeanResidualMaps=True, outDir=mapsDir)

    badAfter = badMaps({k: mapsDir for k in mapFiles}) + ([] if os.path.exists(os.path.join(mapsDir, rawCat)) else [rawCat])
    if len(badAfter) > 0:
        raise RuntimeError("PyBDSF re-extraction failed for %s (missing or truncated: %s). See the messages above."
                           % (baseName, badAfter))

    # the re-extracted sources should be those of the catalogue in the products directory
    nNew = len(atpy.Table.read(os.path.join(mapsDir, rawCat)))
    nOld = len(atpy.Table.read(os.path.join(prodDir, 'catalogs', namedCat)))
    print("Re-extracted catalogue has %d sources, the products catalogue has %d" % (nNew, nOld))
    if nOld > 0 and abs(nNew - nOld) > 0.01 * nOld:
        print("WARNING: the numbers differ by more than 1%. The PyBDSF parameters in sourceinjection.pybdsfArgs "
              "may differ from those used for the products catalogue.")

#------------------------------------------------------------------------------------------------------------
def estimateCompleteness(beginImage=0, endImage=None, nInjectionSources=5000, nRepetitions=50, minFluxJyInj=1e-5, maxFluxJyInj=1.0, radiusFactorToInject=0.95,
                         minSizeArcsec=1.0, maxSizeArcsec=60.0, minAxisRatio=0.3, crowdingFactor=5.0, sinjectionDirName='sinjection'):
    """Estimate completeness as a function of flux density and source size by source injection.

    Sources with random deconvolved sizes (log-uniform between minSizeArcsec and maxSizeArcsec) are injected
    into the residual image of each pointing and recovered with PyBDSF. See sourceinjection.executeSingle.

    """

    # the catalogues with the source names (not the raw PyBDSF ones, which have 'srl' in their names)
    tabFilesList=sorted([f for f in glob.glob(os.path.join(startup.config['rootDir'], 'EDR', 'catalogs', '*_bdsfcat.fits'))
                         if 'srl_bdsfcat' not in os.path.basename(f)])

    sinjectDirPath = os.path.join(startup.config['rootDir'], 'EDR', sinjectionDirName)
    os.makedirs(sinjectDirPath, exist_ok = True)

    if endImage is None:
        selectedTabs = tabFilesList[beginImage:]
    else:
        selectedTabs = tabFilesList[beginImage:endImage]

    # Re-extracted PyBDSF files go to mapsDir, which is also the current directory while the images are processed
    # (sourceinjection.executeSingle looks in the current directory first and in the products directory second)
    mapsDir = os.path.join(sinjectDirPath, 'pybdsf_maps')
    os.makedirs(mapsDir, exist_ok = True)

    originalDir = os.getcwd()
    os.chdir(mapsDir)
    try:
        for tab in selectedTabs:
            baseName = tab.split('_bdsfcat.fits')[0]
            print("\n" + "═" * 70)
            print("║ Estimating completeness for %s ║" %baseName)
            _ensurePyBDSFMaps(os.path.basename(baseName), mapsDir)
            sourceinjection.executeSingle(imageName=baseName,
                nInjectionSources=nInjectionSources,
                nRepetitions=nRepetitions,
                minFluxJyInj=minFluxJyInj,
                maxFluxJyInj=maxFluxJyInj,
                radiusFactorToInject=radiusFactorToInject,
                outDir=sinjectDirPath,
                injectUnresolvedOnly=False,
                minSizeArcsec=minSizeArcsec,
                maxSizeArcsec=maxSizeArcsec,
                minAxisRatio=minAxisRatio,
                crowdingFactor=crowdingFactor)
            print("\n" + "═" * 70)
    finally:
        os.chdir(originalDir)


#------------------------------------------------------------------------------------------------------------
def xmatchSummary(dataRelease='EDR', optSurveyTag='DECaLSDR10', optBand='r', outFileName=None):
    """Collects the likelihood-ratio cross-matching statistics of every cross-matched radio catalogue.

    For each <DRDir>/xmatches_<optSurveyTag>/xmatch_<radio catalogue>_<optSurveyTag>_<optBand>band folder it reads
    LR_THR, REL and COMP (the LR threshold at which reliability and completeness balance, and their values) from
    the header of xmatchtable_bestmatches_*.fits, and Q0 from Q0_*.txt, which is stored separately and is not in
    the header. The search radius (SEAR_RAD) and the number of radio sources with an accepted counterpart
    (rows of the best-matches table) are also recorded.

    Args:
        dataRelease (str): Data release, default is 'EDR'.
        optSurveyTag (str): Tag of the optical survey used in the folder names, default is 'DECaLSDR10'.
        optBand (str): Optical band used for the cross-match, default is 'r'.
        outFileName (str): Output file (FITS). Default: <DRDir>/xmatch_summary_<optSurveyTag>_<optBand>band.fits.
            A CSV copy with the same name is also written.

    Returns:
        astropy.table.Table: One row per cross-matched radio catalogue.
    """

    DRDir = os.path.join(startup.config['rootDir'], dataRelease)
    xmatchRoot = os.path.join(DRDir, 'xmatches_%s' % optSurveyTag)
    suffix = '_%s_%sband' % (optSurveyTag, optBand)

    xmatchDirs = sorted(d for d in glob.glob(os.path.join(glob.escape(xmatchRoot), 'xmatch_*' + suffix))
                        if os.path.isdir(d))   # symlinks to the product directory are followed
    print("\nFound %d cross-match folders in %s" % (len(xmatchDirs), xmatchRoot))

    rows = []
    for xmatchDir in xmatchDirs:
        outSubscript = os.path.basename(xmatchDir)[len('xmatch_'):]
        radCatName = outSubscript[:-len(suffix)] if outSubscript.endswith(suffix) else outSubscript

        bestFile = os.path.join(xmatchDir, 'xmatchtable_bestmatches_%s.fits' % outSubscript)
        q0File = os.path.join(xmatchDir, 'Q0_%s.txt' % outSubscript)

        lrThr = rel = comp = searchRadArcsec = np.nan
        nMatched = nRadio = nOptical = -1
        if os.path.exists(bestFile):
            hdr = fits.getheader(bestFile, ext=1)   # table meta is stored in the header of the binary table
            lrThr = float(hdr.get('LR_THR', np.nan))
            rel = float(hdr.get('REL', np.nan))
            comp = float(hdr.get('COMP', np.nan))
            if 'SEAR_RAD' in hdr:
                searchRadArcsec = float(str(hdr['SEAR_RAD']).split()[0])
            nMatched = int(hdr.get('NAXIS2', -1))
        else:
            print("No best-matches table in %s" % xmatchDir)

        q0 = np.nan
        if os.path.exists(q0File):
            q0 = float(np.atleast_1d(np.loadtxt(q0File))[0])
        else:
            print("No Q0 file in %s" % xmatchDir)

        # LRStats_*.ecsv (written by the cross-match) holds all of these for every pointing, also those without a
        # best-matches table; when it exists it takes precedence
        statsFile = os.path.join(xmatchDir, 'LRStats_%s.ecsv' % outSubscript)
        if os.path.exists(statsFile):
            stats = atpy.Table.read(statsFile, format='ascii.ecsv')
            lrThr, rel, comp, q0 = [float(stats[c][0]) for c in ('LR_THR', 'REL', 'COMP', 'Q0')]
            searchRadArcsec = float(stats['SEAR_RAD_arcsec'][0])
            nRadio, nOptical, nMatched = [int(stats[c][0]) for c in ('N_RADIO', 'N_OPT', 'N_BEST')]

        rows.append((radCatName, lrThr, rel, comp, q0, searchRadArcsec, nMatched, nRadio, nOptical))

    if len(rows) == 0:
        print("\nNo cross-matched catalogues found.")
        return None

    summaryTab = atpy.Table(rows=rows, names=('radCatName', 'LR_THR', 'REL', 'COMP', 'Q0', 'SEAR_RAD_arcsec', 'nMatched', 'nRadio', 'nOptical'),
                            dtype=('U200', 'f8', 'f8', 'f8', 'f8', 'f8', 'i8', 'i8', 'i8'))

    if outFileName is None:
        outFileName = os.path.join(DRDir, 'xmatch_summary_%s_%sband.fits' % (optSurveyTag, optBand))
    summaryTab.write(outFileName, overwrite=True)

    print("\nWrote %s" % outFileName)
    print("\n%-17s %9s %9s %9s %9s" % ('', 'min', 'median', 'max', 'N valid'))
    for col in ('LR_THR', 'REL', 'COMP', 'Q0', 'SEAR_RAD_arcsec'):
        vals = np.asarray(summaryTab[col], dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            print("%-17s %9.3f %9.3f %9.3f %9d" % (col, vals.min(), np.median(vals), vals.max(), len(vals)))
    nBadQ0 = int(np.sum(np.isfinite(summaryTab['Q0']) & ((summaryTab['Q0'] <= 0) | (summaryTab['Q0'] > 1))))
    if nBadQ0 > 0:
        print("\nWarning: Q0 is outside (0, 1] for %d catalogues. It comes from an unconstrained linear fit, so it can "
              "leave the physical range." % nBadQ0)

    return summaryTab


#------------------------------------------------------------------------------------------------------------
def _normaliseCatName(name):
    """Radio catalogue name without directory and without the .fits, _srl_bdsfcat or _bdsfcat suffixes."""
    name = os.path.basename(str(name).strip())
    for suffix in ('.fits', '_srl_bdsfcat', '_bdsfcat'):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name

def _readCatNameList(listFile):
    """Set of normalised catalogue names listed (one per line) in a text file; empty set if the file is missing."""
    if not os.path.exists(listFile):
        print("Warning: %s not found" % listFile)
        return set()
    with open(listFile, 'r', encoding='utf-8') as infile:
        return set(_normaliseCatName(line.split()[0]) for line in infile if line.strip() and not line.startswith('#'))

def _unitVectors(raDeg, decDeg):
    """Cartesian unit vectors of sky positions (no RA wrapping issues)."""
    ra, dec = np.deg2rad(raDeg), np.deg2rad(decDeg)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])

#------------------------------------------------------------------------------------------------------------
def xmatchCoverage(dataRelease='EDR', optSurveyTag='DECaLSDR10', optBand='r', radioBand='L',
                   maxSepArcsec=30.0, maxSearchArcsec=300.0, partialMinFraction=0.01, fullFraction=0.99,
                   nIdentified=None):
    """Finds which unique radio sources lie outside the optical (DECaLS) coverage, fully or partially.

    Coverage is decided from positions. The optical sources fetched for each pointing (the file
    <optSurveyTag>_sources_<radio catalogue>.fits in each xmatch folder, which covers about 4 times the area
    of the pointing) trace the optical footprint. A unique radio source is covered if the nearest optical source
    of ANY pointing's fetch is closer than maxSepArcsec. Positions, not pointings, are used, so this handles
    pointings that are partly covered, pointings without any optical source (listed in
    <DR>_no_<survey>_sources.txt, which have no xmatch folder), pointings with optical sources but no
    overlap (listed in <DR>_no_counterparts_<survey>_<band>band.txt), and overlapping pointings.

    The per-source table and the global numbers use the unique radio sources. The per-pointing fractions use ALL the
    sources of each pointing (survey_catalog_<band>_<DR>.fits), not only the unique ones: in the unique catalogue a
    source in the overlap of two pointings is kept for only one of them, so the unique sources of a pointing can miss
    exactly the part of it that is covered, and its coverage would be underestimated.

    Args:
        dataRelease (str): Data release, default is 'EDR'.
        optSurveyTag (str): Tag of the optical survey used in file and folder names, default is 'DECaLSDR10'.
        optBand (str): Optical band used for the cross-match, default is 'r'.
        radioBand (str): Radio band of the survey catalogue, default is 'L'.
        maxSepArcsec (float): A radio source is covered if an optical source lies within this distance.
            Inside the footprint the distance to the nearest optical source is a few arcsec, outside it grows with the
            distance to the footprint edge.
        maxSearchArcsec (float): Distances are only computed up to this value (larger means 'no optical source nearby').
        partialMinFraction, fullFraction (float): A pointing is called partially covered if the fraction of its
            sources (all of them, including those in overlaps) outside the coverage is between partialMinFraction and
            fullFraction, and not covered above that.
        nIdentified (int): Optional number of unique radio sources with an optical counterpart, to print the
            identification fraction among all and among covered sources.

    Returns:
        tuple: (table with one row per unique radio source, table with one row per pointing). The pointing table has
            nSources, nOutsideCoverage, fractionOutsideCoverage and coverageStatus (all sources of the pointing), and
            nUniqueSources, nUniqueOutsideCoverage (its unique sources).
    """
    from scipy.spatial import cKDTree

    DRDir = os.path.join(startup.config['rootDir'], dataRelease)
    xmatchRoot = os.path.join(DRDir, 'xmatches_%s' % optSurveyTag)
    suffix = '_%s_%sband' % (optSurveyTag, optBand)

    uniq = atpy.Table.read(os.path.join(DRDir, 'survey_catalog_%s_unique_%s.fits' % (radioBand, dataRelease)))
    nUniq = len(uniq)
    uniqXYZ = _unitVectors(np.asarray(uniq['RA'], dtype=float), np.asarray(uniq['DEC'], dtype=float))
    pointingOfSource = np.array([_normaliseCatName(pth) for pth in uniq['radCatPath']])

    allTab = atpy.Table.read(os.path.join(DRDir, 'survey_catalog_%s_%s.fits' % (radioBand, dataRelease)))
    nAll = len(allTab)
    allXYZ = _unitVectors(np.asarray(allTab['RA'], dtype=float), np.asarray(allTab['DEC'], dtype=float))
    pointingOfAll = np.array([_normaliseCatName(pth) for pth in allTab['radCatPath']])
    minSepAll = np.full(nAll, np.inf)

    noSourcesList = _readCatNameList(os.path.join(DRDir, '%s_no_%s_sources.txt' % (dataRelease, optSurveyTag)))
    noCounterpartsList = _readCatNameList(os.path.join(DRDir, '%s_no_counterparts_%s_%sband.txt'
                                                       % (dataRelease, optSurveyTag, optBand)))

    chordMax = 2.0 * np.sin(np.deg2rad(maxSearchArcsec / 3600.0) / 2.0)
    minSepArcsec = np.full(nUniq, np.inf)

    xmatchDirs = sorted(d for d in glob.glob(os.path.join(glob.escape(xmatchRoot), 'xmatch_*' + suffix))
                        if os.path.isdir(d))
    print("\nFound %d cross-match folders in %s" % (len(xmatchDirs), xmatchRoot))

    xmatchedPointings = set()
    for xmatchDir in xmatchDirs:
        outSubscript = os.path.basename(xmatchDir)[len('xmatch_'):]
        radCatName = outSubscript[:-len(suffix)] if outSubscript.endswith(suffix) else outSubscript
        xmatchedPointings.add(_normaliseCatName(radCatName))

        optFile = os.path.join(xmatchDir, '%s_sources_%s.fits' % (optSurveyTag, radCatName))
        if not os.path.exists(optFile):
            print("No %s sources file in %s" % (optSurveyTag, xmatchDir))
            continue

        with fits.open(optFile, memmap=True) as hdul:   # only the two position columns are read
            optRA = np.asarray(hdul[1].data['RADeg'], dtype=float)
            optDec = np.asarray(hdul[1].data['decDeg'], dtype=float)
        if len(optRA) == 0:
            continue

        tree = cKDTree(_unitVectors(optRA, optDec))
        chord, _ = tree.query(uniqXYZ, k=1, distance_upper_bound=chordMax, workers=-1)
        found = np.isfinite(chord)
        sepArcsec = np.full(nUniq, np.inf)
        sepArcsec[found] = np.degrees(2.0 * np.arcsin(chord[found] / 2.0)) * 3600.0
        minSepArcsec = np.minimum(minSepArcsec, sepArcsec)

        chordAll, _ = tree.query(allXYZ, k=1, distance_upper_bound=chordMax, workers=-1)
        foundAll = np.isfinite(chordAll)
        sepAll = np.full(nAll, np.inf)
        sepAll[foundAll] = np.degrees(2.0 * np.arcsin(chordAll[foundAll] / 2.0)) * 3600.0
        minSepAll = np.minimum(minSepAll, sepAll)

    covered = minSepArcsec <= maxSepArcsec
    coveredAll = minSepAll <= maxSepArcsec

    category = np.where(np.isin(pointingOfSource, list(noSourcesList)), 'noOpticalSources',
               np.where(np.isin(pointingOfSource, list(noCounterpartsList)), 'noCounterparts',
               np.where(np.isin(pointingOfSource, list(xmatchedPointings)), 'xmatched', 'notProcessed')))

    # ---- per-source table
    sourceTab = atpy.Table()
    sourceTab['radCatPath'] = uniq['radCatPath']
    sourceTab['RA'] = uniq['RA']
    sourceTab['DEC'] = uniq['DEC']
    sourceTab['minSepOptical_arcsec'] = np.where(np.isfinite(minSepArcsec), minSepArcsec, np.nan)
    sourceTab['hasOpticalCoverage'] = covered
    sourceTab['pointingCategory'] = category

    # ---- per-pointing table (all the sources of each pointing, see the docstring)
    names, inverse = np.unique(pointingOfAll, return_inverse=True)
    nPerPointing = np.bincount(inverse)
    nUncPerPointing = np.bincount(inverse, weights=(~coveredAll).astype(float)).astype(int)
    fracUnc = nUncPerPointing / np.maximum(nPerPointing, 1)
    uniqNames, inverseUniq = np.unique(pointingOfSource, return_inverse=True)
    nUniqPerName = dict(zip(uniqNames, np.bincount(inverseUniq)))
    nUniqUncPerName = dict(zip(uniqNames, np.bincount(inverseUniq, weights=(~covered).astype(float)).astype(int)))
    xmatchedSet, noSourcesSet, noCounterpartsSet = set(xmatchedPointings), set(noSourcesList), set(noCounterpartsList)
    catPerPointing = np.array(['noOpticalSources' if n in noSourcesSet else 'noCounterparts' if n in noCounterpartsSet
                               else 'xmatched' if n in xmatchedSet else 'notProcessed' for n in names])
    status = np.where(fracUnc >= fullFraction, 'notCovered',
             np.where(fracUnc >= partialMinFraction, 'partiallyCovered', 'covered'))
    pointingTab = atpy.Table({'pointing': names, 'category': catPerPointing, 'nSources': nPerPointing,
                              'nOutsideCoverage': nUncPerPointing, 'fractionOutsideCoverage': fracUnc,
                              'nUniqueSources': np.array([nUniqPerName.get(n, 0) for n in names]),
                              'nUniqueOutsideCoverage': np.array([nUniqUncPerName.get(n, 0) for n in names]),
                              'coverageStatus': status})

    outBase = os.path.join(DRDir, 'xmatch_coverage_%s_%sband' % (optSurveyTag, optBand))
    sourceTab.write(outBase + '.fits', overwrite=True)
    pointingTab.write(outBase + '_by_pointing.csv', format='ascii.csv', overwrite=True)

    # ---- report
    nUnc = int((~covered).sum())
    print("\nUnique radio sources                       : %d" % nUniq)
    print("  with %s coverage                    : %d (%.1f%%)" % (optSurveyTag, nUniq - nUnc, 100.0 * (nUniq - nUnc) / nUniq))
    print("  outside the coverage                     : %d (%.1f%%)" % (nUnc, 100.0 * nUnc / nUniq))
    print("\n%-17s %9s %9s %9s %9s" % ('pointing group', 'pointings', 'sources', 'outside', 'covered'))
    for cat in ('noOpticalSources', 'noCounterparts', 'xmatched', 'notProcessed'):
        inCat = category == cat
        if inCat.any():
            print("%-17s %9d %9d %9d %9d" % (cat, int((catPerPointing == cat).sum()), int(inCat.sum()),
                                             int((inCat & ~covered).sum()), int((inCat & covered).sum())))
    print("\nPointings by coverage (fraction of all their sources outside the coverage; partial: %.2f-%.2f):"
          % (partialMinFraction, fullFraction))
    for st in ('notCovered', 'partiallyCovered', 'covered'):
        print("  %-17s %d" % (st, int((status == st).sum())))
    if nIdentified is not None:
        print("\nIdentified sources: %d -> fraction of all unique sources %.1f%%, of covered sources %.1f%%"
              % (nIdentified, 100.0 * nIdentified / nUniq, 100.0 * nIdentified / max(nUniq - nUnc, 1)))
    print("\nWrote %s.fits and %s_by_pointing.csv" % (outBase, outBase))

    return sourceTab, pointingTab


#------------------------------------------------------------------------------------------------------------
def _edrIdentificationStats(dataBase, srcTab, zRange=(0.2, 1.6), radioBand='L'):
    """Numbers of unique radio sources that lie inside the optical (DECaLS) coverage, have an optical ID and a
    redshift. Written to <DR>_identification_statistics.csv in the data release directory (columns: quantity, number,
    percentOfAll, percentOfCovered).

    All counts are of unique radio sources (each physical source counted once, also where pointings overlap):
      N unique        survey_catalog_<band>_unique_<DR>.fits
      covered         a DECaLS source (any pointing's fetch) lies within maxSepArcsec of the source (xmatchCoverage)
      with an ID      xmatchCatUnique_DECaLSDR10_r_<DR>.fits (best matches above the LR threshold, duplicates removed)
      with a redshift zmaxCatUnique_<DR>.fits (spectroscopic if available, otherwise photometric)

    Args:
        dataBase (str): Name of the data release.
        srcTab (astropy.table.Table): Per-source table returned by xmatchCoverage (same order as the unique catalogue).
        zRange (tuple): Redshift range of the luminosity functions, to count the sources inside it.
        radioBand (str): Radio band of the survey catalogue.
    """
    DRDir = os.path.join(startup.config['rootDir'], dataBase)

    nUniq = len(atpy.Table.read(os.path.join(DRDir, 'survey_catalog_%s_unique_%s.fits' % (radioBand, dataBase))))
    covered = np.asarray(srcTab['hasOpticalCoverage'], dtype=bool)
    if len(covered) != nUniq:
        raise ValueError("The coverage table and the unique catalogue have different lengths (%d, %d): run xmatchCoverage "
                         "again after builddr" % (len(covered), nUniq))

    idFile = os.path.join(DRDir, 'xmatchCatUnique_DECaLSDR10_r_%s.fits' % dataBase)
    if os.path.exists(idFile):
        idTab = atpy.Table.read(idFile)
    else:
        print("No %s: removing the duplicates from xmatchCat_DECaLSDR10_r_%s.fits here (run builddr to make it)" % (idFile, dataBase))
        idTab = catalogs.removeDuplicateSources(atpy.Table.read(os.path.join(DRDir, 'xmatchCat_DECaLSDR10_r_%s.fits' % dataBase)),
                                                matchRadius_arcsec=6.0, raCol='RA_rad', decCol='DEC_rad',
                                                fluxCol='Total_flux_rad', fluxErrCol='E_Total_flux_rad')
    zTab = atpy.Table.read(os.path.join(DRDir, 'zmaxCatUnique_%s.fits' % dataBase))
    zBestFile = os.path.join(DRDir, 'zmaxCatUniqueBest_%s.fits' % dataBase)
    zBestTab = atpy.Table.read(zBestFile) if os.path.exists(zBestFile) else None

    zspec = np.asarray(zTab['zspec_opt'], dtype=float)
    hasSpec = np.isfinite(zspec) & (zspec != -99)
    zVals = np.asarray(zTab['zphot_ifnot_spec_opt'], dtype=float)

    nCov, nId, nZ = int(covered.sum()), len(idTab), len(zTab)
    nSpec = int(hasSpec.sum())

    def pct(a, b):
        return 100.0 * a / b if b > 0 else np.nan

    quantities = [
        ('unique radio sources', nUniq),
        ('inside the optical coverage', nCov),
        ('with an optical ID', nId),
        ('with a redshift', nZ),
        ('with a spectroscopic redshift', nSpec),
        ('with a photometric redshift only', nZ - nSpec),
        ('IDs without a redshift', nId - nZ),
        ('with a redshift in %.1f <= z <= %.1f' % zRange, int(np.sum((zVals >= zRange[0]) & (zVals <= zRange[1]))))]
    if zBestTab is not None:
        quantities.append(('in zmaxCatUniqueBest (luminosity-function sample)', len(zBestTab)))
    rows = [(name, n, pct(n, nUniq), np.nan if name == 'unique radio sources' else pct(n, nCov)) for name, n in quantities]

    if nId > nCov:
        print("WARNING: there are more IDs than sources inside the coverage: check maxSepArcsec in xmatchCoverage")

    outFile = os.path.join(DRDir, '%s_identification_statistics.csv' % dataBase)
    atpy.Table(rows=rows, names=('quantity', 'number', 'percentOfAll', 'percentOfCovered'),
               dtype=('U80', 'i8', 'f8', 'f8')).write(outFile, format='ascii.csv', overwrite=True)
    print("Wrote %s" % outFile)


#------------------------------------------------------------------------------------------------------------
def _summarizeEDRXmatches(dataBase, imagesTab):
    """Writes, for the pointings of a data release (e.g. 'EDR'), the lists of pointings without and with partial optical
    (DECaLS) coverage, and a table with the cross-match statistics of every pointing.

    The lists written by the cross-match in the products directory cover all processed pointings, so here they are
    cut down to the pointings of this data release. All files go to <rootDir>/<dataBase>/:

        <DR>_no_DECaLSDR10_sources.txt              pointings with no DECaLS source inside the pointing
        <DR>_no_counterparts_DECaLSDR10_rband.txt   pointings that DECaLS overlaps but with no match (no candidate
                                                    in the search radius, or none above the LR threshold)
        <DR>_partial_DECaLSDR10_coverage.txt        pointings where only part of the radio sources lie inside the
                                                    DECaLS footprint (decided by position, see xmatchCoverage)
        <DR>_xmatch_statistics.fits                 one row per pointing: status of the cross-match, Q0, LR threshold,
                                                    reliability, completeness, search radius, numbers of sources,
                                                    coverage
        <DR>_identification_statistics.csv          numbers of sources inside the coverage, with an ID and a redshift
        xmatch_summary_*.fits, xmatch_coverage_*    (written by xmatchSummary and xmatchCoverage)

    Existing lists with these names are kept as <name>.bak if their content changes.

    Args:
        dataBase (str): Name of the data release.
        imagesTab (astropy.table.Table): Images table of the data release (images_<DR>.fits).
    """
    import shutil

    DRDir = os.path.join(startup.config['rootDir'], dataBase)
    xmatchRoot = os.path.join(startup.config['productsDir'], 'xmatches_DECaLSDR10')

    names = [_normaliseCatName(pth) for pth in imagesTab['radioCatPath'] if str(pth).strip() not in ('-', '')]
    drNames = set(names)

    def writeList(fileName, listNames):
        path = os.path.join(DRDir, fileName)
        text = ''.join('%s_bdsfcat.fits\n' % n for n in sorted(listNames))
        if os.path.exists(path):
            with open(path, 'r', encoding='utf-8') as infile:
                old = infile.read()
            if old != text and not os.path.exists(path + '.bak'):
                shutil.copy2(path, path + '.bak')
        with open(path, 'w', encoding='utf-8') as outfile:
            outfile.write(text)
        print("Wrote %s (%d pointings)" % (path, len(listNames)))

    print("\n" + "═" * 50)
    print((" CROSS-MATCH LISTS AND STATISTICS OF %s " % dataBase).center(50, '═'))
    print("═" * 50)

    # ---- 1. lists restricted to this data release
    noSources = set(n for n in _readCatNameList(os.path.join(xmatchRoot, 'no_DECaLSDR10_sources.txt')) if n in drNames)
    noCounterparts = set(n for n in _readCatNameList(os.path.join(xmatchRoot, 'no_counterparts_DECaLSDR10_rband.txt'))
                         if n in drNames)
    writeList('%s_no_DECaLSDR10_sources.txt' % dataBase, noSources)
    writeList('%s_no_counterparts_DECaLSDR10_rband.txt' % dataBase, noCounterparts)

    # ---- 2. partial coverage, by position (xmatchCoverage reads the two lists written above)
    covTab = None
    srcTab = None
    try:
        srcTab, covTab = xmatchCoverage(dataRelease=dataBase)
    except Exception as e:
        print("\nThe coverage by position could not be computed: %s" % e)
    partial = set()
    if covTab is not None:
        partial = set(str(n) for n in covTab['pointing'][covTab['coverageStatus'] == 'partiallyCovered'])
        writeList('%s_partial_DECaLSDR10_coverage.txt' % dataBase, partial)

    # ---- 2b. identification and redshift statistics of the unique radio sources
    if srcTab is not None:
        try:
            _edrIdentificationStats(dataBase, srcTab)
        except Exception as e:
            import traceback
            traceback.print_exc()
            print("\nThe identification statistics could not be computed: %s" % e)

    # ---- 3. statistics of every pointing
    sumTab = xmatchSummary(dataRelease=dataBase)
    sumBy = {} if sumTab is None else {_normaliseCatName(r['radCatName']): r for r in sumTab}
    covBy = {} if covTab is None else {str(r['pointing']): r for r in covTab}

    rows = []
    for img in imagesTab:
        base = _normaliseCatName(img['radioCatPath'])
        if base in ('-', ''):
            continue
        if base in noSources:
            status = 'no_optical_sources'
        elif base in noCounterparts:
            status = 'no_counterparts'
        elif base in sumBy:
            status = 'matched'
        else:
            status = 'not_run'
        sm, cv = sumBy.get(base), covBy.get(base)
        rows.append((base, str(img['object']), str(img['captureBlockId']), float(img['centre_RADeg']),
                     float(img['centre_decDeg']), float(img['skyArea_sqDeg']), float(img['RMS_uJy/beam']),
                     int(img['nRadioSources']), status,
                     float(sm['Q0']) if sm is not None else np.nan, float(sm['LR_THR']) if sm is not None else np.nan,
                     float(sm['REL']) if sm is not None else np.nan, float(sm['COMP']) if sm is not None else np.nan,
                     float(sm['SEAR_RAD_arcsec']) if sm is not None else np.nan,
                     int(sm['nOptical']) if sm is not None else -1, int(sm['nMatched']) if sm is not None else -1,
                     int(cv['nUniqueSources']) if cv is not None else -1,
                     float(cv['fractionOutsideCoverage']) if cv is not None else np.nan,
                     str(cv['coverageStatus']) if cv is not None else ''))
    statNames = ('pointing', 'object', 'captureBlockId', 'centre_RADeg', 'centre_decDeg', 'skyArea_sqDeg', 'RMS_uJy/beam',
                 'nRadioSources', 'xmatchStatus', 'Q0', 'LR_THR', 'REL', 'COMP', 'SEAR_RAD_arcsec', 'nOptical', 'nMatched',
                 'nUniqueSources', 'fractionOutsideCoverage', 'coverageStatus')
    statDtypes = ('U200', 'U100', 'U30', 'f8', 'f8', 'f8', 'f8', 'i8', 'U30', 'f8', 'f8', 'f8', 'f8', 'f8', 'i8', 'i8', 'i8',
                  'f8', 'U30')
    statTab = atpy.Table(rows=rows, names=statNames, dtype=statDtypes)
    statBase = os.path.join(DRDir, '%s_xmatch_statistics' % dataBase)
    statTab.write(statBase + '.fits', overwrite=True)
    print("Wrote %s.fits (%d pointings)" % (statBase, len(statTab)))

    # ---- 5. summary
    print("\nCross-match status of the %d pointings of the %s:" % (len(statTab), dataBase))
    for status in ('matched', 'no_counterparts', 'no_optical_sources', 'not_run'):
        print("  %-20s %d" % (status, int(np.sum(statTab['xmatchStatus'] == status))))
    if covTab is not None:
        print("Coverage by position:")
        for status in ('covered', 'partiallyCovered', 'notCovered'):
            print("  %-20s %d" % (status, int(np.sum(statTab['coverageStatus'] == status))))
    for col in ('Q0', 'LR_THR', 'REL', 'COMP'):
        vals = np.asarray(statTab[col], dtype=float)
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            print("  %-7s median %.3f, range %.3f to %.3f (%d pointings)" % (col, np.median(vals), vals.min(), vals.max(), len(vals)))