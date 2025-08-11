"""

Routines for cross matching, using likelihood ratio method.

"""

import os
import numpy as np
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.units import Quantity
from astropy import units as u
from astropy.table import Table, hstack
import matplotlib.pyplot as plt
from zCluster import retrievers
from . import startup, catalogs
from collections import defaultdict
from astropy.units import Quantity


def _writeNoOpticalCatFile(catalogName, outFileName):
    """
    Append or create a text file to record the name of a radio catalogue with no optical counterpart.

    Args:
        catalogName (:obj:`str`): Name or identifier of the radio catalogue without an optical match.
        outFileName (:obj:`str`): Path to the output text file where the catalogue name will be recorded.

    Returns:
        :obj:`int`: Returns 0 upon successful write operation.
    """

    if os.path.exists(outFileName):
        with open(outFileName, 'a', encoding='utf8') as outFile:
            outFile.write("%s\n" %(catalogName))
    else:
        with open(outFileName, 'w', encoding='utf8') as outFile:
            outFile.write("%s\n" %(catalogName))
    return 0

def filterBadRows(table, columnsToCheck, badValues=[99., 999., -99., -999.]):
    """
    Filters rows in an Astropy table where any of the specified columns contain bad placeholder values.

    Args:
        table (astropy.table.Table): The input table to filter.
        columnsToCheck (list of str): List of column names to check for bad values.
        badValues (list of float, optional): Values that should be considered invalid. Defaults to [99., 999., -99., -999.].

    Returns:
        astropy.table.Table: A filtered table with bad rows removed.
    """

    # Start with all rows as valid
    mask = np.ones(len(table), dtype=bool)

    for col in columnsToCheck:
        values = table[col]
        # Convert to raw numeric values if unit is present
        try:
            values = values.value
        except AttributeError:
            pass
        mask &= ~np.isin(values, badValues)

    nRemoved = len(table) - np.sum(mask)
    print("\nRemoved %d bad rows out of %d!" % (nRemoved, len(table)))

    return table[mask]

def _getUnitlessValues(col):
    """
    Returns unitless numeric values from an astropy Table Column.
    If the column has a unit (not None), returns .value.
    Otherwise, returns the column as-is.
    """
    if hasattr(col, 'unit') and col.unit is not None:
        return col.value
    else:
        return col

def makeMagBins(magnitudes, nBins):
    """ Compute custom magnitude bin edges for a given array of magnitudes.

    Args:
        magnitudes (:obj:`array_like`): Array of magnitudes to be binned.
        nBins (:obj:`int`): Desired number of bins.

    Returns:
        magBins (:obj:`np.ndarray`): Array of bin edges.
    """
    minMag = np.min(magnitudes)
    maxMag = np.max(magnitudes)

    magBins = np.linspace(minMag, maxMag, nBins + 1)

    return magBins


def getMagHist(optMag, nBins):
    """Compute histogram of optical magnitudes.

    Args:
        optMag (:obj:`array_like`): Array of optical magnitudes (e.g., from parent optical catalogue).
        nBins (:obj:`int`): Number of bins to use for the histogram.

    Returns:
        tuple:
            - optMagHist (:obj:`np.ndarray`): Counts in each magnitude bin.
            - optMagBins (:obj:`np.ndarray`): Edges of the magnitude bins.
    """

    optMagBins = makeMagBins(optMag, nBins)

    optMagHist, _ = np.histogram(optMag, bins=optMagBins)

    return optMagHist, optMagBins

def getNM(optMag, nBins, areaSqArcsec):
    """Estimate the magnitude distribution n(m) from an optical catalogue.

    Args:
        optMag (:obj:`array_like`): Array of optical magnitudes (e.g., from parent optical catalogue).
        nBins (:obj:`int`): Number of bins to use for the magnitude histogram.
        areaSqArcsec (:obj:`astropy.units.Quantity`): Survey area in square arcseconds (with units attached).

    Returns:
        :obj:`np.ndarray`: The magnitude distribution n(m), i.e., number of sources per square arcsecond per bin.
    """

    optMagHist, optMagBins = getMagHist(optMag, nBins)

    if hasattr(areaSqArcsec, 'unit'):
        areaSqArcsec = areaSqArcsec.value

    binWidth = np.diff(optMagBins)[0]
    distNm = optMagHist / (areaSqArcsec * binWidth)
    return distNm

def randomPointsInCircle(centerRA, centerDec, radiusDeg, nPoints):
    """
    Generate uniformly distributed random points within a circular area on the sky.

    Args:
        centerRA (:obj:`float`): Right Ascension of the circle's center (degrees).
        centerDec (:obj:`float`): Declination of the circle's center (degrees).
        radiusDeg (:obj:`float`): Radius of the circular area (degrees).
        nPoints (:obj:`int`): Number of random points to generate.

    Returns:
        tuple:
            - randRA (array): Array of RA values (degrees) for the random points.
            - randDec (array): Array of Dec values (degrees) for the random points.
    """

    centerCoord = SkyCoord(ra=centerRA * u.deg, dec=centerDec * u.deg, frame='icrs')

    randRad = np.sqrt(np.random.uniform(0, radiusDeg**2, nPoints))  # Radius
    randTheta = np.random.uniform(0, 2 * np.pi, nPoints)  # Angle

    randRA = centerRA + randRad * np.cos(randTheta) / np.cos(np.radians(centerDec))
    randDec = centerDec + randRad * np.sin(randTheta)

    randCoords = SkyCoord(ra=randRA * u.deg, dec=randDec * u.deg, frame='icrs')
    toCentreDistances = centerCoord.separation(randCoords)
    withinCircle = toCentreDistances < radiusDeg * u.deg
    return randRA[withinCircle], randDec[withinCircle]


def makeRandomCat(centerRA, centerDec, radiusDeg, nRandomPoints, radRACol, radDecCol):
    """
    Create a random catalogue of sources uniformly distributed within a circular sky region.

    Args:
        centerRA (:obj:`float`): RA of the circle center (degrees).
        centerDec (:obj:`float`): Dec of the circle center (degrees).
        radiusDeg (:obj:`float`): Radius of the circle (degrees).
        nRandomPoints (:obj:`int`): Total number of random points to generate.
        radRACol (:obj:`str`): Column name for RA in the output table.
        radDecCol (:obj:`str`): Column name for Dec in the output table.

    Returns:
        :obj:`~astropy.table.Table`: Astropy Table with columns `[radRACol, radDecCol]` containing the random positions.
    """

    randRA, randDec = randomPointsInCircle(centerRA, centerDec, radiusDeg, nRandomPoints)
    randomCat = Table([randRA, randDec], names=(radRACol, radDecCol))
    return randomCat


def countBlanks(skyCatCoords1, skyCatCoords2, searchRadRsArcsec):
    """
    Counts the number of sources in skyCatCoords1 that have no match within searchRadRsArcsec in skyCatCoords2.

    Parameters:
    - skyCatCoords1 : SkyCoord of first catalogue.
    - skyCatCoords2 : SkyCoord of second catalogue.
    - searchRadRsArcsec: search radius in arcseconds

    Returns:
    - nBlanks: number of skyCatCoords1 positions with no match in skyCatCoords2
    """

    idx2, sep2d, _ = skyCatCoords1.match_to_catalog_sky(skyCatCoords2)

    nBlanks = np.sum(sep2d.arcsec > searchRadRsArcsec)

    return nBlanks

def getQ0(radCatCoords, randRadCatCoords, optCatCoords, searchRadRsArcsec, sigmaRad):
    """
    Estimate Q0, the fraction of radio sources with real counterparts, using the blank fields method.

    Args:
        radCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the radio catalogue positions.
        randRadCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the randomised radio catalogue positions.
        optCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the optical/IR catalogue positions.
        searchRadRsArcsec (:obj:`float`): Search radius in arcseconds within which counterparts are considered.
        sigmaRad (:obj:`float`): Typical positional uncertainty (sigma) of radio sources in arcseconds.

    Returns:
        float: Estimated Q0 value representing the fraction of radio sources with true counterparts.

    Raises:
        ValueError: If the number of blank fields in the random catalogue (nBlankRand) is zero, preventing division by zero.
    """

    nBlankReal = countBlanks(radCatCoords, optCatCoords, searchRadRsArcsec)

    nBlankRand = countBlanks(randRadCatCoords, optCatCoords, searchRadRsArcsec)

    if nBlankRand == 0:
        print("nBlankRand is zero — cannot divide by zero when computing Q0. Setting Q0 = 1.0")
        return 1.0

    Frs = 1 - np.exp( -0.5 * (searchRadRsArcsec**2 / sigmaRad**2))

    if Frs == 0:
        print("Frs is zero — cannot divide by zero when computing Q0. Setting Q0 = 1")
        return 1.0

    Q0 = (1 - nBlankReal/nBlankRand)/Frs

    print("\nQ0 = ", Q0)
    return Q0

def getQM(radCatCoords, nRadio, optCatCoords, optMagList, searchRadRmaxArcsec, nM, Q0):

    """
    Estimate q(m): the magnitude distribution of real counterparts to radio sources.

    Args:
        radCatCoords (SkyCoord): SkyCoord array of radio source positions.
        nRadio (int): Number of radio sources.
        optCatCoords (SkyCoord): SkyCoord array of optical catalogue positions.
        optMagList (array_like): List or array of optical magnitudes corresponding to optCatCoords.
        searchRadRmaxArcsec (float): Maximum search radius in arcseconds.
        nM (array_like): Surface density of background sources per magnitude bin, in sources/arcsec².
        Q0 (float): Fraction of radio sources that have real optical counterparts (0 < Q0 < 1).

    Returns:
        np.ndarray: Estimated q(m), the probability distribution of magnitudes for real counterparts.
    """
    nMagBins = len(nM)

    # getting total(m)

    searchRadRmaxDeg = searchRadRmaxArcsec / 3600.

    # Find all pairs within r_max
    _, idxOpt, _, _ = search_around_sky(radCatCoords, optCatCoords, searchRadRmaxDeg * u.deg)

    if len(idxOpt) == 0:
        print("\nNo cross-matches...!")
        return None

    # idxOpt are indices of optical sources within r_max of any radio source
    matchedMags = optMagList[idxOpt]

    totalM, _ = getMagHist(optMag=matchedMags, nBins=nMagBins)

    areaPerSource = np.pi * searchRadRmaxArcsec**2  # in sq.arcsec
    backgroundCounts = nM * nRadio * areaPerSource

    realM = totalM - backgroundCounts
    realM[realM < 0] = 0.  # prevent negative values

    if realM.sum() > 0:
        qM = (realM / realM.sum()) * Q0
    else:
        qM = np.zeros_like(realM)

    return qM

def getValueFromMagBins(magnitude, magBins, binValues):
    """
    Returns the corresponding bin value for a given magnitude based on magnitude bins.

    Parameters
    ----------
    magnitude : float
        The magnitude value for which the bin value is to be returned.
    magBins : array-like
        The edges of the magnitude bins. Length must be N+1 for N binValues.
    binValues : array-like
        The value associated with each bin (e.g. q(m) or n(m)). Length must be N.

    Returns
    -------
    float
        The value corresponding to the bin in which the magnitude falls.
        Returns np.nan if the value is outside the bin range or if inputs are inconsistent.
    """
    if len(binValues) != len(magBins) - 1:
        raise ValueError("Length of bin_values must be one less than length of bin_edges.")

    binIndex = np.digitize(magnitude, magBins, right=True) - 1
    if 0 <= binIndex < len(binValues):
        return binValues[binIndex]
    return np.nan

def getCentreRadius(radioCat, radRACol, radDecCol):

    """Calculate the approximate center and maximum radius of a radio source catalogue footprint.

    Args:
        radioCat (:obj:`~astropy.table.Table`): Astropy Table containing radio source catalogue.
        radRACol (:obj:`str`): Name of the column for right ascension (degrees) in `radioCat`.
        radDecCol (:obj:`str`): Name of the column for declination (degrees) in `radioCat`.

    Returns:
        tuple:
            - center_ra (:obj:`float`): RA of the bounding-box center (degrees).
            - center_dec (:obj:`float`): Dec of the bounding-box center (degrees).
    """

    raColEntries = radioCat[radRACol]
    decColEntries = radioCat[radDecCol]

    # Check if RA and Dec has units; if yes, extract raw values
    if isinstance(raColEntries, Quantity):
        raList = raColEntries.value
    else:
        raList = raColEntries

    if isinstance(decColEntries, Quantity):
        decList = decColEntries.value
    else:
        decList = decColEntries

    raList = _getUnitlessValues(radioCat[radRACol])
    decList = _getUnitlessValues(radioCat[radDecCol])
    catCoords = SkyCoord(ra=raList * u.deg, dec=decList * u.deg, frame='icrs')

    # Centre from bounding box
    centerRA = (np.min(raList) + np.max(raList)) / 2
    centerDec = (np.min(decList) + np.max(decList)) / 2
    centerCoord = SkyCoord(ra=centerRA * u.deg, dec=centerDec * u.deg, frame='icrs')

    # True maximum angular distance to any source
    toCentreseparations = centerCoord.separation(catCoords)
    radiusDegMax = np.max(toCentreseparations).deg

    return centerRA, centerDec, radiusDegMax

def retrieveDECaLS(centerRA, centerDec, radiusDeg, DR='DR10'):
    """  Retrieve DECaLS sources within a circular region around a given sky position.

    Args:
        centerRA (:obj:`float`): Right ascension (RA) of the centre of the search region, in degrees.
        centerDec (:obj:`float`): Declination (Dec) of the centre of the search region, in degrees.
        radiusDeg (:obj:`float`): Search radius around the central position, in degrees.
        DR (:obj:`str`, optional): DECaLS data release to use. Currently only 'DR10' is supported. Default is 'DR10'.

    Returns:
        :obj:`astropy.table.Table`: Table of DECaLS sources within the specified region.

    Raises:
        :obj:`Exception`: If a data release other than 'DR10' is requested.
    """

    print("\nRetrieving DECaLS %s sources with RA_central=%.2f deg, Dec_central=%.2f deg, and radius=%.2f deg" \
          % (DR, centerRA, centerDec, radiusDeg))
    if DR == 'DR10':
        decalsCat = Table(retrievers.DL_DECaLSDR10Retriever(centerRA, centerDec,
                                                             halfBoxSizeDeg = radiusDeg,
                                                             DR = None))

    else:
        raise Exception("DR to be DR10")

    return decalsCat

def getFR(rOffsetArcsec, radioSource, opticalSource, radRACol, radDecCol, radEMajCol, radEMinCol, radPACol, optRACol, optDecCol, optPosErrCol, sigmaAst=0.6):

    """
    Calculate the probability distribution f(r) of offset r between radio and a potential counterpart.

    Args:
        rOffsetArcsec (:obj:`float` or :obj:`np.ndarray`): Angular offset (separation) between radio and optical positions,
            in arcsec.
        radioSource (:obj:`~astropy.table.Row`): A single row from the radio catalogue table.
        opticalSource (:obj:`~astropy.table.Row`): A single row from the optical catalogue table.
        radRACol (:obj:`str`): Key for right ascension of the radio source in `radioSource`.
        radDecCol (:obj:`str`): Key for declination of the radio source in `radioSource`.
        radEMajCol (:obj:`str`): Key for major axis FWHM error of the radio source Gaussian fit.
        radEMinCol (:obj:`str`): Key for minor axis FWHM error of the radio source Gaussian fit.
        radPACol (:obj:`str`): Key for position angle (PA, degrees east of north) of the radio source major axis.
        optRACol (:obj:`str`): Key for right ascension of the optical source in `opticalSource`.
        optDecCol (:obj:`str`): Key for declination of the optical source in `opticalSource`.
        optPosErrCol (:obj:`str`): Key for positional uncertainty of the optical source.
        sigmaAst (:obj:`float`, optional): Astrometric uncertainty between radio and optical surveys (default 0.6 arcsec).

    Returns:
        float or np.ndarray:
            The value of the positional probability distribution f(r) for the given offset(s).

    """

    deltaMaj = radioSource[radEMajCol]
    deltaMin = radioSource[radEMinCol]
    sigmaMajRad = deltaMaj/np.sqrt(4 * np.log(2))
    sigmaMinRad = deltaMin/np.sqrt(4 * np.log(2))

    # getting vector joining radio to optical counterpart

    RARad = radioSource[radRACol]
    decRad = radioSource[radDecCol]
    RAOpt = opticalSource[optRACol]
    decOpt = opticalSource[optDecCol]

    dRA = (RAOpt - RARad) * np.cos(np.radians(decRad))
    dDec = decOpt - decRad

    thetaDir = np.arctan2(dDec, dRA) # in radians

    # angle between major axis and radio-optical vector

    positionAngle = np.radians(radioSource[radPACol])
    thetaPADir = thetaDir - positionAngle

    sigmaDirRad = np.sqrt(
        (sigmaMajRad * np.cos(thetaPADir))**2 +
        (sigmaMinRad * np.sin(thetaPADir))** 2
    )

    # calculating sigmas for optical

    sigmaOpt = opticalSource[optPosErrCol]
    sigmaRAOpt = sigmaOpt
    sigmaDecOpt = sigmaOpt # Only one position uncertainty for optical
    sigmaMajOpt = np.sqrt(
        (sigmaRAOpt * np.sin(positionAngle))**2 +
        (sigmaDecOpt * np.cos(positionAngle))**2
    )
    positionAngleMin = positionAngle + np.pi/2.
    sigmaMinOpt = np.sqrt(
        (sigmaRAOpt * np.sin(positionAngleMin))**2 +
        (sigmaDecOpt * np.cos(positionAngleMin))**2
    )

    sigmaDirOpt = np.sqrt(
        (sigmaRAOpt * np.cos(thetaDir))**2 +
        (sigmaDecOpt * np.sin(thetaDir))**2
    )

    sigmaMaj = np.sqrt(sigmaMajRad**2 + sigmaMajOpt**2 + sigmaAst**2)
    sigmaMin = np.sqrt(sigmaMinRad**2 + sigmaMinOpt**2 + sigmaAst**2)
    sigmaDir = np.sqrt(sigmaDirRad**2 + sigmaDirOpt**2 + sigmaAst**2)

    probDistR = (1 / (2 * np.pi * sigmaMaj * sigmaMin)) * np.exp(-0.5 * (rOffsetArcsec**2 / sigmaDir**2))

    return probDistR  # f(r)

def computeRelCompl(LRTab, Q0, NRadio, LRThreshold):
    """
    Compute the completeness and reliability for a given likelihood ratio (LR) threshold.

    Args:
        LRTab (:obj:`astropy.table.Table`): Table containing likelihood ratio values for cross-matched sources.
            Must include a column named 'LR'.
        Q0 (:obj:`float`): The fraction of true counterparts expected among all radio sources (i.e., the prior).
        NRadio (:obj:`int`): Total number of radio sources in the catalogue.
        LRThreshold (:obj:`float`): Likelihood ratio threshold above which identifications are considered reliable.

    Returns:
        :obj:`tuple`: Tuple containing:
            - completeness (:obj:`float`): Fraction of real identifications above the LR threshold.
            - reliability (:obj:`float`): Fraction of accepted identifications that are correct.
    """


    LRVals = LRTab['LR'].value

    # Completeness: sum over LR_i < L_thr
    thresholdMaskCompl = LRVals < LRThreshold
    complSum = np.sum((Q0 * LRVals[thresholdMaskCompl]) / (Q0 * LRVals[thresholdMaskCompl] + (1 - Q0)))

    completeness = 1 - complSum / (Q0 * NRadio)

    # Reliability: sum over LR_i >= L_thr
    thresholdMaskRel = LRVals >= LRThreshold
    relSum = np.sum((1 - Q0) / (Q0 * LRVals[thresholdMaskRel] + (1 - Q0)))

    reliability = 1 - relSum / (Q0 * NRadio)

    return completeness, reliability

def computeLR(radioCat, opticalCat, searchRadiusArcsec, optMagCol, magBins, qMList, nMList, radRACol, radDecCol, optRACol, optDecCol, radEMajCol, radEMinCol, radPACol, optPosErrCol):
    """
    Compute the Likelihood Ratio (LR) for matches between radio and optical sources.

    For each radio source, finds optical candidates within the search radius and calculates
    the LR = (q(m)/n(m)) * f(r) for each candidate, where:
      - f(r): positional probability density based on offsets and positional errors,
      - q(m): magnitude distribution of true counterparts,
      - n(m): magnitude distribution of background objects.

    Returns a merged Astropy Table containing columns from both catalogs for each matched pair,
    plus additional columns 'f_r', 'q_m', 'n_m', and 'LR'.

    Args:
        radioCat (:obj:`~astropy.table.Table`): Radio source catalogue.
        opticalCat (:obj:`~astropy.table.Table`): Optical source catalogue.
        searchRadiusArcsec (:obj:`float`): Search radius around radio sources in arcseconds.
        optMagCol (:obj:`str`): Column name for optical magnitudes in opticalCat.
        magBins (:obj:`array_like`): Bin edges used for q(m) and n(m) calculation.
        qMList (:obj:`array_like`): q(m) values for magnitude bins.
        nMList (:obj:`array_like`): n(m) values for magnitude bins.
        radRACol (:obj:`str`): RA column name in radioCat.
        radDecCol (:obj:`str`): Dec column name in radioCat.
        optRACol (:obj:`str`): RA column name in opticalCat.
        optDecCol (:obj:`str`): Dec column name in opticalCat.
        radEMajCol (:obj:`str`): Major axis error column name in radioCat.
        radEMinCol (:obj:`str`): Minor axis error column name in radioCat.
        radPACol (:obj:`str`): Position angle column name in radioCat.
        optPosErrCol (:obj:`str`): Positional error column name in opticalCat.

    Returns:
        astropy.table.Table: Merged table with one row per matched pair, containing all columns
        from both input tables (optical and radio columns postfixed with '_opt' and '_rad' resp.) and columns 'f_r', 'q_m', 'n_m', and 'LR'.
    """

    radRAList = _getUnitlessValues(radioCat[radRACol])
    radDecList = _getUnitlessValues(radioCat[radDecCol])

    optRAList = _getUnitlessValues(opticalCat[optRACol])
    optDecList = _getUnitlessValues(opticalCat[optDecCol])

    radCatCoords = SkyCoord(ra= radRAList * u.deg,
                               dec=radDecList * u.deg)
    optCatCoords = SkyCoord(ra= optRAList * u.deg,
                               dec= optDecList * u.deg)

    # Find all pairs within r_max
    idxRadio, idxOpt, _, _ = search_around_sky(radCatCoords, optCatCoords, searchRadiusArcsec * u.arcsec)

    if len(idxOpt) == 0:
        print("\nNo cross-matches...!")
        return None

    rowsRadio = []
    rowsOptical = []
    radOptSeparation = []
    fRVals = []
    qMVals = []
    nMVals = []
    LRVals = []

    for rIdx, oIdx in zip(idxRadio, idxOpt):
        radioSource = radioCat[rIdx]
        opticalSource = opticalCat[oIdx]



        radioCoord = SkyCoord(ra=radRAList[rIdx] * u.deg,
                               dec=radDecList[rIdx] * u.deg)
        opticalCoord = SkyCoord(ra=optRAList[oIdx] * u.deg,
                                     dec=optDecList[oIdx] * u.deg)

        radOptOffsetArcsec = radioCoord.separation(opticalCoord).to(u.arcsec).value

        fRPair = getFR(radOptOffsetArcsec, radioSource, opticalSource, radRACol, radDecCol, radEMajCol, radEMinCol, radPACol, optRACol, optDecCol, optPosErrCol, sigmaAst=0.6)

        opticalMagnitude = opticalSource[optMagCol]

        qMPair = getValueFromMagBins(opticalMagnitude, magBins, qMList)

        nMPair = getValueFromMagBins(opticalMagnitude, magBins, nMList)

        if nMPair == 0 or np.isnan(nMPair) or np.isnan(qMPair):
            # Avoid division by zero or NaNs
            LRPair = np.nan
        else:
            LRPair = (qMPair / nMPair) * fRPair

        # Collect for merged table
        rowsRadio.append(dict(radioSource))
        rowsOptical.append(dict(opticalSource))
        radOptSeparation.append(radOptOffsetArcsec)
        fRVals.append(fRPair)
        qMVals.append(qMPair)
        nMVals.append(nMPair)
        LRVals.append(LRPair)

    # Create tables for matched pairs
    radioMatchesTab = Table(rowsRadio)
    opticalMatchesTab = Table(rowsOptical)

    # Rename columns to avoid clashes
    for col in opticalMatchesTab.colnames:
        opticalMatchesTab.rename_column(col, f"{col}_opt")
    for col in radioMatchesTab.colnames:
        radioMatchesTab.rename_column(col, f"{col}_rad")

    # Merge radio + optical tables horizontally
    radOptMergedTab = hstack([radioMatchesTab, opticalMatchesTab])

    radOptMergedTab['rad_opt_sep_asec'] = radOptSeparation
    radOptMergedTab['f_r'] = fRVals
    radOptMergedTab['f_r'] = fRVals
    radOptMergedTab['q_m'] = qMVals
    radOptMergedTab['n_m'] = nMVals
    radOptMergedTab['LR'] = LRVals

    return radOptMergedTab

def xmatchRadioOptical(radioCatFilePath, radioBand, xmatchDirPath, optSurvey, optSurveyDR, optMagCol, searchRadiusArcsec, makePlots, radRACol, radDecCol, radERACol, radEDecCol, radEMajCol, radEMinCol, radPACol, outSubscript, optPosErrAsecValue=0.2, nMagBins=15, beamSizeArcsecValue=6.0, saveFiles = True, skipIfExists=True):
    """
    Perform likelihood ratio crossmatching between a radio source catalog and an optical survey.

    For each radio source, this function finds candidate optical counterparts within the
    search radius and computes the Likelihood Ratio (LR) for each candidate, incorporating:
      - positional uncertainties of radio and optical sources,
      - magnitude distributions of true counterparts (q(m)),
      - background magnitude distributions (n(m)),
      - and spatial probability density f(r) based on offsets.

    Args:
        radioCatFilePath (:obj:`str`): Path to the input radio catalog FITS file.
        radioBand (:obj:`str`): Radio frequency band identifier ('L', 'UHF', 'S').
        xmatchDirPath (:obj:`str`): Directory to save all resulting output files and plots.
        optSurvey (:obj:`str`): Optical survey name (e.g., 'DECaLS').
        optSurveyDR (:obj:`str`): Optical survey data release identifier (e.g., 'DR10').
        optMagCol (:obj:`str`): Column name for optical magnitudes used in LR calculation.
        searchRadiusArcsec (:obj:`float`): Search radius around radio positions in arcseconds.
        makePlots (:obj:`bool`): Whether to generate diagnostic plots.
        radRACol (:obj:`str`): Radio catalog Right Ascension column name.
        radDecCol (:obj:`str`): Radio catalog Declination column name.
        radERACol (:obj:`str`): Radio catalog RA positional error column name.
        radEDecCol (:obj:`str`): Radio catalog Dec positional error column name.
        radEMajCol (:obj:`str`): Radio catalog major axis positional error column name.
        radEMinCol (:obj:`str`): Radio catalog minor axis positional error column name.
        radPACol (:obj:`str`): Radio catalog position angle column name.
        outSubscript (:obj:`str`): String appended to output filenames.
        optPosErrAsecValue (:obj:`float`, optional): Assumed optical positional error in arcseconds. Default is 0.2.
        nMagBins (:obj:`int`, optional): Number of magnitude bins used for magnitude distribution estimation. Default is 15.
        beamSizeArcsecValue (:obj:`float`, optional): Radio beam size in arcseconds. Default is 6.0.
        saveFiles (:obj:`bool`, optional): Whether to save output files. Default is True.
        skipIfExists (:obj:`bool`, optional): Skip processing if output files exist. Default is True.

    Returns:
        astropy.table.Table: Table of best crossmatched sources with LR values and Reliability and Completeness in meta.
    """

    radCatName = radioCatFilePath.split(os.path.sep)[-1]
    captureBlockId = radCatName.split('_')[3]
    targetName = (radCatName.split('_1024ch_')[1]).split('_srl_')[0]
    xmatchIndividualDirPath = os.path.join(xmatchDirPath, 'xmatch_%s' %outSubscript)
    xmatchTabName = xmatchIndividualDirPath+os.path.sep+"xmatchtable_%s" %outSubscript+".fits"
    xmatchBestMatchTabName = xmatchIndividualDirPath+os.path.sep+"xmatchtable_bestmatches_%s" %outSubscript+".fits"
    os.makedirs(xmatchIndividualDirPath, exist_ok = True)

    print("\n" + "═" * 100)
    print("║ Cross-matching %s ║" %radCatName)
    print("═" * 100 + "\n")

    if skipIfExists is True:
        doesItExist = os.path.exists(xmatchBestMatchTabName) and os.path.exists(xmatchTabName)
        if doesItExist:
            print("%s cross-match table exists! Skipping to next radio catalogue.\n" %radCatName)
            xmatchTable = Table.read(xmatchBestMatchTabName)
            return xmatchTable

    # checking if this catalog is listed as having no optical sources in the area
    noOptSourcesFilename = xmatchDirPath+os.path.sep+'no_%s_%s_sources.txt' %(optSurvey, optSurveyDR)

    if os.path.exists(noOptSourcesFilename):
        with open(noOptSourcesFilename, 'r', encoding="utf-8") as infile:
            noOptSourcesCatNames = [line.strip() for line in infile]
        if radCatName in noOptSourcesCatNames:
            return None

    # checking if this catalog is listed as having no optical counterparts in the optical
    noOptCounterpartsFilename = xmatchDirPath+os.path.sep+'no_counterparts_%s%s_%sband_%sasec.txt' \
                           %(optSurvey, optSurveyDR, optMagCol, str(searchRadiusArcsec).replace(".","p"))

    if os.path.exists(noOptCounterpartsFilename):
        with open(noOptCounterpartsFilename, 'r', encoding="utf-8") as infile:
            noOptCounterpartsCatNames = [line.strip() for line in infile]
        if radCatName in noOptCounterpartsCatNames:
            return None

    if optSurvey.lower() != 'decals':
        print("\nERROR: Code not working for %s survey." %optSurvey)
        # Currently configured only for DECaLS...
        return None

    optPosErrValueDeg = (optPosErrAsecValue*u.arcsec).to(u.deg).value

    radioSources = Table.read(radioCatFilePath, format='fits', hdu=1)

    nRadio = len(radioSources)
    print("\nNumber of radio sources: %d" %nRadio)

    radColumns = [radRACol, radDecCol, radERACol, radEDecCol, radEMajCol, radEMinCol, radPACol]

    if not all(col in radioSources.colnames for col in radColumns):
        print("Required columns of radio catalogue is not well set. Exiting!")
        return

    if any(radioSources[radRACol] < 0.0):
        radioSources = catalogs.fixRA(radioSources, raCol=radRACol, wrapAngle=360)

    radRAValDegList = _getUnitlessValues(radioSources[radRACol])
    radDecValDegList = _getUnitlessValues(radioSources[radDecCol])

    radioSourcesCoords = SkyCoord(ra=radRAValDegList * u.deg, dec=radDecValDegList * u.deg, frame='icrs')

    centerRA, centerDec, radiusDeg = getCentreRadius(radioSources, radRACol, radDecCol)

    skyAreaSqDeg = np.pi*radiusDeg**2
    skyAreaSqArcsec = skyAreaSqDeg*3600**2

    sigmaRadPos = np.sqrt(radioSources[radERACol]**2 + radioSources[radEDecCol]**2)
    sigmaRadPosMean = np.mean(sigmaRadPos)
    sigmaRadPosMeanArcsec = sigmaRadPosMean * 3600.

    # Collecting optical sources
    optRACol, optDecCol = 'RADeg', 'decDeg'
    optPosErrCol = 'pos_err'

    optCatFileName = xmatchIndividualDirPath+os.path.sep+'decals_sources_%s.fits' %outSubscript
    if os.path.exists(optCatFileName):
        opticalSources = Table.read(optCatFileName, format='fits', hdu=1)
    else:
        opticalSourcesRaw = retrieveDECaLS(centerRA, centerDec, radiusDeg, DR=optSurveyDR)
        if len(opticalSourcesRaw) == 0:
            print("\n%s: No optical sources found in %s database...!" %(radCatName, optSurvey))
            _writeNoOpticalCatFile(radCatName, noOptCounterpartsFilename)
            return

        # filtering optical catalogue
        opticalSources = filterBadRows(opticalSourcesRaw, [optMagCol])

    if len(opticalSources) == 0:
        print("\n%s: No optical sources with reliable %s-magnitude found in %s database...!" %(radCatName, optMagCol, optSurvey))
        _writeNoOpticalCatFile(radCatName, noOptCounterpartsFilename)
        return

    print("\nNumber of %s%s sources with reliable %s-magnitude in the sky region: %d" % (optSurvey, optSurveyDR, optMagCol, len(opticalSources)))

    if optPosErrCol not in opticalSources.colnames:
        opticalSources[optPosErrCol] = optPosErrValueDeg

    optMagList = opticalSources[optMagCol]
    optMagBins = makeMagBins(optMagList, nMagBins)
    optMagBinCenters = 0.5 * (optMagBins[1:] + optMagBins[:-1])

    optRAValDegList = _getUnitlessValues(opticalSources[optRACol])
    optDecValDegList = _getUnitlessValues(opticalSources[optDecCol])

    optSourcesCoords = SkyCoord(ra=optRAValDegList * u.deg, dec= optDecValDegList * u.deg, frame='icrs')

    nRandomRadSources = 1 * nRadio # TODO make sure about number of randoms. need to normalize somewhere if different?
    randomRadioSources=makeRandomCat(centerRA, centerDec, radiusDeg, nRandomRadSources, radRACol, radDecCol)

    randRadRAValDegList = _getUnitlessValues(randomRadioSources[radRACol])
    randRadDecValDegList = _getUnitlessValues(randomRadioSources[radDecCol])

    randomRadioSourcesCoords = SkyCoord(ra= randRadRAValDegList * u.deg, dec= randRadDecValDegList * u.deg, frame='icrs')

    #Q0 = 0.7024573416919023 #TODO

    Q0 = getQ0(radCatCoords=radioSourcesCoords, randRadCatCoords=randomRadioSourcesCoords, optCatCoords=optSourcesCoords, searchRadRsArcsec=beamSizeArcsecValue, sigmaRad=sigmaRadPosMeanArcsec)

    # Finding n(m)

    nM = getNM(optMag=optMagList, nBins=nMagBins, areaSqArcsec=skyAreaSqArcsec)

    # Finding q(m)

    qM = getQM(radCatCoords=radioSourcesCoords, nRadio=nRadio, optCatCoords=optSourcesCoords, optMagList=optMagList, searchRadRmaxArcsec=beamSizeArcsecValue, nM=nM, Q0=Q0)
    if qM is None:
        _writeNoOpticalCatFile(radCatName, noOptCounterpartsFilename)
        return None

    print("\nComputing LR ...")

    xmatchTable = computeLR(radioCat=radioSources, opticalCat=opticalSources, searchRadiusArcsec=searchRadiusArcsec, optMagCol=optMagCol, magBins=optMagBins, qMList=qM, nMList=nM, radRACol=radRACol, radDecCol=radDecCol, optRACol=optRACol, optDecCol=optDecCol, radEMajCol=radEMajCol, radEMinCol=radEMinCol, radPACol=radPACol, optPosErrCol=optPosErrCol)

    if xmatchTable is None:
        print("\n%s: No cross-matched objects...!" %radCatName)
        _writeNoOpticalCatFile(radCatName, noOptCounterpartsFilename)
        return

    if 'captureBlockId' not in xmatchTable.columns:
        xmatchTable.add_column(captureBlockId, name='captureBlockId', index=0)
    if 'object' not in xmatchTable.columns:
        xmatchTable.add_column(targetName, name='object', index=1)
    if 'band' not in xmatchTable.columns:
        xmatchTable.add_column(radioBand, name='band')

    # Finding reliability and completeness

    print("\nComputing Completeness and Reliability ...")

    LRThreshValues = []
    completenessValues = []
    reliabilityValues = []

    for LRThresh in np.arange(0.02, 10.5, 0.02):
        completeness, reliability = computeRelCompl(LRTab=xmatchTable, Q0=Q0, NRadio=nRadio, LRThreshold=LRThresh)

        LRThreshValues.append(LRThresh)
        completenessValues.append(completeness)
        reliabilityValues.append(reliability)

    deltaCR = np.abs(np.array(completenessValues) - np.array(reliabilityValues))
    CRBalanceLRThresholdIndex = np.argmin(deltaCR)
    CRBalanceLRThreshold = LRThreshValues[CRBalanceLRThresholdIndex]
    CRBalanceRel = reliabilityValues[CRBalanceLRThresholdIndex]
    CRBalanceComp = completenessValues[CRBalanceLRThresholdIndex]

    print("\nLR threshold with almost equal reliability (%0.2f) and completeness (%0.2f) = %0.2f"
          %(CRBalanceRel, CRBalanceComp, CRBalanceLRThreshold))

    xmatchLRThresholdTable = xmatchTable[xmatchTable['LR'] >= CRBalanceLRThreshold]

    if len(xmatchLRThresholdTable) == 0:
        print("\nNo cross-matched sources above LR threshold!")
        _writeNoOpticalCatFile(radCatName, noOptCounterpartsFilename)
        return None


    # Filters the full xmatch table to keep only the matches with the maximum 'LR' value for each radio source
    xmatchLRThresholdTable.sort(['Source_id_rad', 'LR'])
    xmatchLRThresholdTable.reverse()
    groupedxmatchLRThresholdTable = xmatchLRThresholdTable.group_by('Source_id_rad')
    xmatchBestMatchTable = groupedxmatchLRThresholdTable.groups.aggregate(lambda rows: rows[0])

    xmatchBestMatchTable.meta['OPT_SUR']='%s%s' %(optSurvey, optSurveyDR)
    xmatchBestMatchTable.meta['SEAR_RAD']='%f arcsec' %searchRadiusArcsec
    xmatchBestMatchTable.meta['LR_THR']=CRBalanceLRThreshold
    xmatchBestMatchTable.meta['REL']=CRBalanceRel
    xmatchBestMatchTable.meta['COMP']=CRBalanceComp

    if saveFiles is True:

        opticalSources.write(optCatFileName, format='fits', overwrite=True)

        randomRadioSources.write(xmatchIndividualDirPath+os.path.sep+'Randoms_%s.fits' %outSubscript, format='fits', overwrite=True)

        np.savetxt(xmatchIndividualDirPath+os.path.sep+'Q0_%s.txt' %outSubscript, [Q0], fmt='%f')

    if makePlots is True:
        # Plotting optical and radio sources
        RADecplotOutName = "%s/RadOptSkyPlot_%s.png" %(xmatchIndividualDirPath, outSubscript)
        plt.figure(figsize=(6, 6))
        plt.scatter(opticalSources[optRACol], opticalSources[optDecCol], s=1,
                    c='#8AD5F1', label='%s (N=%d)'%(optSurvey, len(opticalSources)))
        plt.scatter(radioSources[radRACol], radioSources[radDecCol], s=2, c='#87340D',
                    label='MeerKAT (N=%d)' %len(radioSources))
        #plt.scatter(xmatchTable[radRACol+'_rad'], xmatchTable[radDecCol+'_rad'],
        #            marker='o', facecolor='None', linewidth=0.5, s=15,
        #            edgecolor='#06471D',
        #            label='MeerKATx%s (N=%d)' %(optSurvey, len(xmatchTable)))
        plt.scatter(xmatchBestMatchTable[radRACol+'_rad'], xmatchBestMatchTable[radDecCol+'_rad'],
                    marker='o', facecolor='None', linewidth=0.5, s=15,
                    edgecolor='#06471D',
                    label='MeerKATx%s (Best matches; N=%d)' %(optSurvey, len(xmatchBestMatchTable)))
        plt.title(radCatName)
        plt.title("%s\nSearch radius = %0.1f asec, %s band, Q0=%0.2f" \
                    % (radCatName, searchRadiusArcsec, optMagCol, Q0))
        plt.xlabel("RA (deg; J2000)")
        plt.ylabel("Dec (deg; J2000)")
        plt.legend(loc="lower left", scatterpoints=1, fontsize=10)
        plt.savefig(RADecplotOutName , dpi=300, bbox_inches = 'tight')
        plt.close()
        print("\nPlotted sky coverage of radio and optical sources!")

        # Plotting qm/nm as a function of magnitude

        qmNmPlotOutName = "%s/QmNmPlot_%s.png" %(xmatchIndividualDirPath, outSubscript)
        plt.figure(figsize=(8, 5))

        # avoiding error due to nM = 0
        qMnMRatio = np.full_like(qM, np.nan, dtype=float)
        mask = nM != 0
        qMnMRatio[mask] = qM[mask] / nM[mask]

        plt.plot(optMagBinCenters, qMnMRatio)
        plt.xlabel(optMagCol+'-band magnitude')
        plt.ylabel(r'$q(m)/n(m)$')
        plt.tight_layout()
        plt.savefig(qmNmPlotOutName, dpi=300, bbox_inches = 'tight')
        plt.close()
        print("\nPlotted q(m)/n(m) vs magnitude!")

        # Plotting LR and Rel
        LRRelPlotOutName = "%s/LRRelMagPlot_%s.png" %(xmatchIndividualDirPath, outSubscript)
        plt.figure(figsize=(8, 5))
        plt.plot(LRThreshValues, completenessValues, label='Completeness', color='blue')
        plt.plot(LRThreshValues, reliabilityValues, label='Reliability', color='green')
        plt.gca().axvline(x=CRBalanceLRThreshold, linestyle='dashed', color='k', label='Choosen LR Threshold')
        plt.xlabel('LR Threshold')
        plt.ylabel('Completeness/Reliability')
        plt.legend()
        plt.tight_layout()
        plt.savefig(LRRelPlotOutName , dpi=300, bbox_inches = 'tight')
        plt.close()
        print("\nPlotted Reliability/Completeness vs LR_threshold!")


    xmatchTable.write(xmatchTabName, format='fits', overwrite=True)
    print("\nWrote full cross-matched table %s." %xmatchTabName)

    xmatchBestMatchTable.write(xmatchBestMatchTabName, format='fits', overwrite=True)
    print("\nWrote cross-matched table with Reliability %0.2f and completeness %0.2f: %s." % (CRBalanceRel, CRBalanceComp, xmatchBestMatchTabName))

    return xmatchBestMatchTable
