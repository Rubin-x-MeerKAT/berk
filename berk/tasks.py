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

    zmaxFilePath = os.path.join(zMaxDirPath, 'zmax_%s.fits' %radCatName.split('_srl_bdsfcat')[0])

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

    rmsAreaFileName = radCatName.replace('_srl_bdsfcat.fits', '_rms_rmshist.txt')
    rmsAreaFile = os.path.join(DRDir, 'rms', rmsAreaFileName)
    rmsHistFromFile = np.loadtxt(rmsAreaFile)
    rmsBinsCentres, cumArea = rmsHistFromFile[:, 0], rmsHistFromFile[:, 4]
    zmaxTab['EffectiveAreaSqDeg'] = summaryPlots.getEffectiveAreaInFluxBinsfromRMS(zmaxTab['Total_flux_rad'].value,
                                                   rmsBinsCentres, cumArea, sigmaDetection=5.0)
    zmaxTab['WholeAreaSqDeg'] = np.max(cumArea)

    DRImageRow = DRImages[DRImages['radioCatPath'] == 'catalogs/'+radCatName]

    if len(DRImageRow) == 0:
        return None, None

    # Radio completeness: as a function of flux density and size if the 2D completeness is available,
    # otherwise as a function of flux density only
    radCatBaseName = radCatName.split('_srl_bdsfcat')[0]
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


        pybdsfCommonName = os.path.basename(row['radioCatPath']).split('_srl_bdsfcat')[0]

        catFile = os.path.join(startup.config['productsDir'], 'catalogs', pybdsfCommonName+'_bdsfcat.fits')
        srlCatFile = os.path.join(startup.config['productsDir'], 'catalogs', pybdsfCommonName+'_srl_bdsfcat.fits')
        rmsFile = os.path.join(startup.config['productsDir'], 'rms', pybdsfCommonName+'_rms.fits')
        rmsHistFile = os.path.join(startup.config['productsDir'], 'rms', pybdsfCommonName+'_rms_rmshist.txt')

        xmatchDir = os.path.join(startup.config['productsDir'],
                                 'xmatches_DECaLSDR10',
                                 'xmatch_%s_srl_bdsfcat_DECaLSDR10_rband' %pybdsfCommonName)

        if os.path.exists(fitsFile) and not os.path.exists(os.path.join(DRDir, row['path'])):
            fitsFileRelPath = os.path.relpath(fitsFile, start=DRImageDir)
            pngFileRelPath = os.path.relpath(pngFile, start=DRImageDir)
            os.system("ln -s %s %s" %(fitsFileRelPath, DRImageDir))
            os.system("ln -s %s %s" %(pngFileRelPath, DRImageDir))
        if os.path.exists(catFile) and not os.path.exists(os.path.join(DRDir, row['radioCatPath'].replace('_srl',''))):
            catFileRelPath = os.path.relpath(catFile, start=DRCatDir)
            srlCatFileRelPath = os.path.relpath(srlCatFile, start=DRCatDir)
            os.system("ln -s %s %s" %(catFileRelPath, DRCatDir))
            os.system("ln -s %s %s" %(srlCatFileRelPath, DRCatDir))
        if os.path.exists(rmsFile) and not os.path.exists(os.path.join(DRRmsDir, os.path.basename(rmsFile))):
            rmsFileRelPath = os.path.relpath(rmsFile, start=DRRmsDir)
            rmsHistFileRelPath = os.path.relpath(rmsHistFile, start=DRRmsDir)
            os.system("ln -s %s %s" %(rmsFileRelPath, DRRmsDir))
            os.system("ln -s %s %s" %(rmsHistFileRelPath, DRRmsDir))
        if os.path.exists(xmatchDir) and not os.path.exists(os.path.join(DRXmatchDir, os.path.basename(xmatchDir))):
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
                DRzmaxTabUniqueBest = DRzmaxTabUniqueBest[DRzmaxTabUniqueBest['radCatName'] != 'catalogs/'+catAvoidName]
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
        radCatName = "catalogs"+os.path.sep+commonFilenamePartForCat+"_srl_bdsfcat.fits"
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
    radCatFilesList=sorted(glob.glob(startup.config['productsDir']+os.path.sep+'catalogs'+os.path.sep+'*srl_bdsfcat.fits'))

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

    Args:
            dataRelease (str): The name of the data release to build, default is 'EDR'.

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

    rmsAreaCoveragePlotName = dataBaseDir+os.path.sep+'MeerKAT_RMS_area_coverage%s.png' % subScript

    summaryPlots.plotRMSAreaCoverageCumulative(rmsDirPath, rmsAreaCoveragePlotName, bandColorDict, nRMSBins=100)

    # Plotting RMS area coverage taking into account overlapping pointings

    # Finding overlapping pointings
    overlapGroupsByBand, isolatedIndicesByBand = images.findOverlappingPointings(imagesTab,  raCol='centre_RADeg', decCol='centre_decDeg', areaCol='skyArea_sqDeg')

    rmsAreaCoverageUniqPlotName = dataBaseDir+os.path.sep+'MeerKAT_RMS_area_coverage_unique%s.png' % subScript

    summaryPlots.plotRMSAreaCoverageCumulativeUnique(rmsDirPath, rmsAreaCoverageUniqPlotName, bandColorDict, imagesTab, overlapGroupsByBand, isolatedIndicesByBand, nRMSBins=100)

    # # Plotting sourcecount

    sourceCountPlotName = dataBaseDir+os.path.sep+'MeerKAT_sourcecount%s.png' % subScript
    summaryPlots.plotSourceCounts(dataBaseDir, catSubScript=subScript, fluxCol='Total_flux', fluxMin=1E-5, fluxMax=1.0, nFluxBins=39, bandColorDict=bandColorDict, plotOutPath=sourceCountPlotName)

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
def estimateCompleteness(beginImage=0, endImage=None, nInjectionSources=5000, nRepetitions=50, minFluxJyInj=1e-5, maxFluxJyInj=1.0, radiusFactorToInject=0.95,
                         minSizeArcsec=1.0, maxSizeArcsec=60.0, minAxisRatio=0.3, crowdingFactor=5.0, sinjectionDirName='sinjection'):
    """Estimate completeness as a function of flux density and source size by source injection.

    Sources with random deconvolved sizes (log-uniform between minSizeArcsec and maxSizeArcsec) are injected
    into the residual image of each pointing and recovered with PyBDSF. See sourceinjection.executeSingle.

    """

    tabFilesList=sorted(glob.glob(os.path.join(startup.config['rootDir'], 'EDR', 'catalogs', '*srl_bdsfcat.fits')))

    sinjectDirPath = os.path.join(startup.config['rootDir'], 'EDR', sinjectionDirName)
    os.makedirs(sinjectDirPath, exist_ok = True)

    if endImage is None:
        selectedTabs = tabFilesList[beginImage:]
    else:
        selectedTabs = tabFilesList[beginImage:endImage]

    for tab in selectedTabs:
        baseName = tab.split('_srl_bdsfcat.fits')[0]
        print("\n" + "═" * 70)
        print("║ Estimating completeness for %s ║" %baseName)
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
        nMatched = -1
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

        rows.append((radCatName, lrThr, rel, comp, q0, searchRadArcsec, nMatched))

    if len(rows) == 0:
        print("\nNo cross-matched catalogues found.")
        return None

    summaryTab = atpy.Table(rows=rows, names=('radCatName', 'LR_THR', 'REL', 'COMP', 'Q0', 'SEAR_RAD_arcsec', 'nMatched'),
                            dtype=('U200', 'f8', 'f8', 'f8', 'f8', 'f8', 'i8'))

    if outFileName is None:
        outFileName = os.path.join(DRDir, 'xmatch_summary_%s_%sband.fits' % (optSurveyTag, optBand))
    summaryTab.write(outFileName, overwrite=True)
    summaryTab.write(os.path.splitext(outFileName)[0] + '.csv', format='ascii.csv', overwrite=True)

    print("\nWrote %s (and .csv)" % outFileName)
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

    Args:
        dataRelease (str): Data release, default is 'EDR'.
        optSurveyTag (str): Tag of the optical survey used in file and folder names, default is 'DECaLSDR10'.
        optBand (str): Optical band used for the cross-match, default is 'r'.
        radioBand (str): Radio band of the survey catalogue, default is 'L'.
        maxSepArcsec (float): A radio source is covered if an optical source lies within this distance.
            Check the distribution of nearest-neighbour separations in the plot written by this function: inside the
            footprint it peaks at a few arcsec, outside it grows with the distance to the footprint edge.
        maxSearchArcsec (float): Distances are only computed up to this value (larger means 'no optical source nearby').
        partialMinFraction, fullFraction (float): A pointing is called partially covered if the fraction of its unique
            sources outside the coverage is between partialMinFraction and fullFraction, and not covered above that.
        nIdentified (int): Optional number of unique radio sources with an optical counterpart, to print the
            identification fraction among all and among covered sources.

    Returns:
        tuple: (table with one row per unique radio source, table with one row per pointing).
    """
    from scipy.spatial import cKDTree
    from matplotlib.figure import Figure

    DRDir = os.path.join(startup.config['rootDir'], dataRelease)
    xmatchRoot = os.path.join(DRDir, 'xmatches_%s' % optSurveyTag)
    suffix = '_%s_%sband' % (optSurveyTag, optBand)

    uniq = atpy.Table.read(os.path.join(DRDir, 'survey_catalog_%s_unique_%s.fits' % (radioBand, dataRelease)))
    nUniq = len(uniq)
    uniqXYZ = _unitVectors(np.asarray(uniq['RA'], dtype=float), np.asarray(uniq['DEC'], dtype=float))
    pointingOfSource = np.array([_normaliseCatName(pth) for pth in uniq['radCatPath']])

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

    covered = minSepArcsec <= maxSepArcsec

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

    # ---- per-pointing table
    names, inverse = np.unique(pointingOfSource, return_inverse=True)
    nPerPointing = np.bincount(inverse)
    nUncPerPointing = np.bincount(inverse, weights=(~covered).astype(float)).astype(int)
    fracUnc = nUncPerPointing / np.maximum(nPerPointing, 1)
    catPerPointing = np.array([category[np.where(inverse == i)[0][0]] for i in range(len(names))])
    status = np.where(fracUnc >= fullFraction, 'notCovered',
             np.where(fracUnc >= partialMinFraction, 'partiallyCovered', 'covered'))
    pointingTab = atpy.Table({'pointing': names, 'category': catPerPointing, 'nUniqueSources': nPerPointing,
                              'nOutsideCoverage': nUncPerPointing, 'fractionOutsideCoverage': fracUnc,
                              'coverageStatus': status})

    outBase = os.path.join(DRDir, 'xmatch_coverage_%s_%sband' % (optSurveyTag, optBand))
    sourceTab.write(outBase + '.fits', overwrite=True)
    pointingTab.write(outBase + '_by_pointing.csv', format='ascii.csv', overwrite=True)

    # ---- histogram of the nearest-neighbour separations (to choose maxSepArcsec)
    fig = Figure(figsize=(7, 4.5))
    ax = fig.subplots()
    sepForPlot = np.where(np.isfinite(minSepArcsec), minSepArcsec, maxSearchArcsec * 1.05)
    ax.hist(sepForPlot, bins=np.linspace(0, maxSearchArcsec * 1.05, 120), color='steelblue')
    ax.axvline(maxSepArcsec, color='r', ls='--', label='maxSepArcsec = %.0f"' % maxSepArcsec)
    ax.set_yscale('log')
    ax.set_xlabel('Distance to the nearest %s source / arcsec (last bin: none within %.0f")' % (optSurveyTag, maxSearchArcsec))
    ax.set_ylabel('Number of unique radio sources')
    ax.legend(frameon=False)
    fig.savefig(outBase + '_nnsep.png', dpi=150, bbox_inches='tight')

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
    print("\nPointings by coverage (fraction of their unique sources outside the coverage; partial: %.2f-%.2f):"
          % (partialMinFraction, fullFraction))
    for st in ('notCovered', 'partiallyCovered', 'covered'):
        print("  %-17s %d" % (st, int((status == st).sum())))
    if nIdentified is not None:
        print("\nIdentified sources: %d -> fraction of all unique sources %.1f%%, of covered sources %.1f%%"
              % (nIdentified, 100.0 * nIdentified / nUniq, 100.0 * nIdentified / max(nUniq - nUnc, 1)))
    print("\nWrote %s.fits, %s_by_pointing.csv and %s_nnsep.png" % (outBase, outBase, outBase))
    print("Check the histogram: the threshold should separate the peak of covered sources from the outside tail.")

    return sourceTab, pointingTab