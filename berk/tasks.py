"""

Tasks that can be performed by Berk

"""

import os
import sys
import subprocess
import glob
import datetime
import astropy.table as atpy
from . import startup, jobs, catalogs, images,  __version__, crossmatch, summaryPlots
import shlex
import matplotlib.pyplot as plt
import numpy as np

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

#------------------------------------------------------------------------------------------------------------
def builddb():
    """Build database...

    """

    # Build global catalog in each band (L-band, UHF)
    globalTabsDict={'L': None, 'UHF': None, 'S': None}

    # Fixing RA
    tabFilesList=glob.glob(startup.config['productsDir']+os.path.sep+'catalogs'+os.path.sep+'*_bdsfcat.fits')
    for t in tabFilesList:
        if t.find("srl_bdsfcat") == -1:
            tab=atpy.Table().read(t)
            if any(tab['RA'] < 0.0):
                tab = catalogs.fixRA(tab, raCol='RA', wrapAngle=360)
            freqGHz=tab.meta['FREQ0']/1e9
            bandKey=getBandKey(freqGHz)
            #tab.meta=None # It'd be good to clear this... but the catalog matching stuff wants many things from here
            #tab.meta.clear() #
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

    imgFilesList=glob.glob(startup.config['productsDir']+os.path.sep+"images"+os.path.sep+"pbcorr_*.fits")
    statsDictList=[]

    for imgFile in imgFilesList:
        pathName = imgFile.replace(startup.config['productsDir']+os.path.sep, '')

        if (pathName in processedImages):
            continue #avoiding recomputing statistics

        statDict=images.getImagesStats(imgFile)
        captureBlockId=os.path.split(statDict['path'])[-1].split('img_')[-1].split('_sdp')[0]
        statDict['captureBlockId']=captureBlockId
        statDict['path']=pathName
        statDict['band']=getBandKey(statDict['freqGHz'])
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
    if qualTab is not None:
        for irow in imgTab:
            mask=irow['path'] == qualTab['path']
            if 'quality' in qualTab.keys():
               #pulling already existing quality and assigning 99 to newly added images
                if mask.any():
                    irow['quality']=qualTab[mask]['quality'][0]
                else:
                    irow['quality'] = 99

    # Output
    imgTab.meta['BERKVER']=__version__
    imgTab.meta['DATEMADE']=datetime.date.today().isoformat()
    imgTab.write(outFileName, overwrite = True)
    print("\nWrote %s" % (outFileName))
    imgTab.write(qualFileName, overwrite = True)

    # Generate survey mask in some format - we'll use that to get total survey area

def xmatch():
    """Does cross-matching...

    """

    optSurvey = 'DECaLS'
    optSurveyDR = 'DR10'
    optBandToMatch = 'r'
    searchRadiusArcsec = 4.0

    globalBestXmatchTabName=startup.config['productsDir']+os.path.sep+"xmatchCat_%s%s_%s_%sasec.fits" %(optSurvey,
                                                                                                    optSurveyDR,
                                                                                                    optBandToMatch,
                                                                                                    str(searchRadiusArcsec).replace('.', 'p'))


    globalBestXmatchTab=None
    radCatFilesList=sorted(glob.glob(startup.config['productsDir']+os.path.sep+'catalogs'+os.path.sep+'*srl_bdsfcat.fits'))

    for radCat in radCatFilesList:
        catalogName = radCat.split(os.path.sep)[-1]
        captureBlockId = catalogName.split('_')[3]

        # making a subscript for this particular match
        outSubscript = '%s_%s%s_%sband_%sasec' %(catalogName.replace('.fits',''), optSurvey, optSurveyDR, optBandToMatch, str(searchRadiusArcsec).replace(".","p"))

        xmatchDirPath = os.path.join(startup.config['productsDir'], 'xmatches')
        os.makedirs(xmatchDirPath, exist_ok = True)

        radCatTab = atpy.Table().read(radCat)
        freqGHz=radCatTab.meta['FREQ0']/1e9
        radBandName=getBandKey(freqGHz)

        xmatchTab = crossmatch.xmatchRadioOptical(radioCatFilePath=radCat,
                                                    radioBand=radBandName,
                                                    xmatchDirPath=xmatchDirPath,
                                                    optSurvey=optSurvey,
                                                    optSurveyDR=optSurveyDR,
                                                    optMagCol = optBandToMatch,
                                                    searchRadiusArcsec = searchRadiusArcsec,
                                                    makePlots=True,
                                                    radRACol='RA', radERACol='E_RA',
                                                    radDecCol='DEC', radEDecCol='E_DEC',
                                                    radEMajCol='E_Maj',
                                                    radEMinCol='E_Min', radPACol='PA',
                                                    outSubscript=outSubscript,
                                                    optPosErrAsecValue=0.2, nMagBins=15, beamSizeArcsecValue=6.0,
                                                    saveFiles = True, skipIfExists=True
                                                    )

        if xmatchTab:
            if globalBestXmatchTab is None:
                globalBestXmatchTab = xmatchTab
            else:
                globalBestXmatchTab = atpy.vstack([globalBestXmatchTab, xmatchTab])

    globalBestXmatchTab.write(globalBestXmatchTabName, overwrite = True)
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
                stubs.append(line)
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

    # Get rms
    print("Collecting rms images...")
    toPath=startup.config['productsDir']+os.path.sep+"rms"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        rmsPath="processing/*/IMAGES/pbcorr_trim_*_pybdsf/*_rms.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+rmsPath, toPath)
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
        cmd=cmd+"\npython3 sourcefinding.py c %s -o fits --survey MSS" % (pbcorrImgPath)

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
def summarize():
    """Summarize the progress in berk processing.

    """
    print("\n" + "═" * 50)
    print("║ SUMMARY OF MEERKAT DATA PROCESSING ║".center(50))
    print("═" * 50 + "\n")


    imagesFileName = startup.config['productsDir']+os.path.sep+"images.fits"

    imagesTab = atpy.Table().read(imagesFileName)

    bandColorDict={'L': '#e35c1e', 'UHF': '#1e21e3', 'S': '#145a32'}

    # Plotting sky coverage

    skyPlotName = startup.config['productsDir']+os.path.sep+'MeerKAT_pointings.png'

    summaryPlots.plotSkyCoverage(imagesTab, bandColorDict=bandColorDict, plotOutPath=skyPlotName, plotProjection='aitoff')

    # Plotting sourcecount

    sourceCountPlotName = startup.config['productsDir']+os.path.sep+'MeerKAT_sourcecount.png'
    summaryPlots.plotSourceCounts(imagesTab, fluxCol='Total_flux', bandColorDict=bandColorDict, nFluxBins=39, plotOutPath=sourceCountPlotName)

    # Plotting RMS coverage
    #rmsCoveragePlotName = startup.config['productsDir']+os.path.sep+'MeerKAT_RMS_coverage.png'
    #summaryPlots.plotRMSCoverage(rmsCoveragePlotName)

    # Plotting RMS area coverage

    rmsAreaCoveragePlotName = startup.config['productsDir']+os.path.sep+'MeerKAT_RMS_area.png'

    summaryPlots.plotRMSAreaCoverage(rmsAreaCoveragePlotName, bandColorDict)

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
