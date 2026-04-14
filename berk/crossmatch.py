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
import requests
import pyvo as vo
from dl import authClient as ac
from getpass import getpass
from astropy_healpix import HEALPix

def DLLogin(username=None, password=None, max_attempts=3):
    """
    Logs into Data Lab and returns a token.
    Prompts user if username/password are not provided.
    Retries up to max_attempts times if login fails.
    """
    username = username or os.getenv("DL_USERNAME")
    password = password or os.getenv("DL_PASSWORD")

    if username is None:
        username = input("Enter Data Lab username: ")

    attempt = 0
    while attempt < max_attempts:
        if password is None:
            password = getpass("Enter Data Lab password: ")

        try:
            token = ac.login(username, password)
            print("\nLogin to %s successful!" %ac.whoAmI())
            return token
        except Exception as e:
            print("\nLogin failed:", e)
            password = None  # prompt again on next loop
            attempt += 1

    # If we get here, all attempts failed
    raise RuntimeError("Failed to login to Data Lab after %d attempts." % max_attempts)

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

def getNM(optMag, nBins, areaSqDeg):
    """Estimate the magnitude distribution n(m) from an optical catalogue.

    Args:
        optMag (:obj:`array_like`): Array of optical magnitudes (e.g., from parent optical catalogue).
        nBins (:obj:`int`): Number of bins to use for the magnitude histogram.
        areaSqDeg (:obj:`astropy.units.Quantity`): Survey area in square degrees (with units attached).

    Returns:
        :obj:`np.ndarray`: The magnitude distribution n(m), i.e., number of sources per square degree per bin.
    """

    optMagHist, optMagBins = getMagHist(optMag, nBins)

    if hasattr(areaSqDeg, 'unit'):
        areaSqDeg = areaSqDeg.value

    binWidth = np.diff(optMagBins)[0]
    distNm = optMagHist / (areaSqDeg * binWidth)
    return distNm

def randomPointsInCircleExactNPoints(centerRA, centerDec, radiusDeg, nPoints, rng=None):
    """
    Generate uniformly distributed random points within a circular area on the sky.

    Args:
        centerRA (:obj:`float`): Right Ascension of the circle's center (degrees).
        centerDec (:obj:`float`): Declination of the circle's center (degrees).
        radiusDeg (:obj:`float`): Radius of the circular area (degrees).
        nPoints (:obj:`int`): Number of random points to generate.
        rng (:obj:`numpy.random.Generator`, optional): A NumPy random number generator instance for reproducibility. 

    Returns:
        tuple:
            - randRA (array): Array of RA values (degrees) for the random points.
            - randDec (array): Array of Dec values (degrees) for the random points.
    """

    if rng is None:
        rng = np.random.default_rng()

    center = SkyCoord(centerRA*u.deg, centerDec*u.deg, frame='icrs')

    # Proper spherical random distribution
    randR = radiusDeg * np.sqrt(rng.uniform(0, 1, nPoints))
    randTheta = rng.uniform(0, 2*np.pi, nPoints)

    # Offset using spherical geometry
    coords = center.directional_offset_by(randTheta * u.rad,
                                          randR * u.deg)

    return coords.ra.deg, coords.dec.deg

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

    randRA, randDec = randomPointsInCircleExactNPoints(centerRA, centerDec, radiusDeg, nRandomPoints)
    randomCat = Table([randRA, randDec], names=(radRACol, radDecCol))
    return randomCat


def countBlanks(skyCatCoords1, skyCatCoords2, searchRadRsDeg):
    """
    Counts the number of sources in skyCatCoords1 that have no match within searchRadRsDeg in skyCatCoords2.

    Parameters:
    - skyCatCoords1 : SkyCoord of first catalogue.
    - skyCatCoords2 : SkyCoord of second catalogue.
    - searchRadRsDeg: search radius in degrees

    Returns:
    - nBlanks: number of skyCatCoords1 positions with no match in skyCatCoords2
    """

    idx2, sep2d, _ = skyCatCoords1.match_to_catalog_sky(skyCatCoords2)

    nBlanks = np.sum(sep2d.deg > searchRadRsDeg)

    return nBlanks

def getQ0FitForRs(radCatCoords, randRadCatCoords, optCatCoords, searchRadiusDeg, sigmaRadDeg):
    """
    Estimate Q0, the fraction of radio sources with real counterparts, using the blank fields method
    by fitting for a series of radii.

    Args:
        radCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the radio catalogue positions.
        randRadCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the randomised radio catalogue positions.
        optCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the optical/IR catalogue positions.
        searchRadRsDeg (:obj:`float`): Search radius in deg within which counterparts are considered (used to get maximum radii).
        sigmaRadDeg (:obj:`float`): Typical positional uncertainty (sigma) of radio sources in deg.

    Returns:
        tuple:
            - Q0 (float) : Value of Q0 after fitting
            - radiiDeg (array): Array of radii used for ratio estimation and fitting
            - UobsByUrandArray (array): Array of Uobs/Urandom ratio for the radii
            - FrsArray (array): Array of F(r) obtained for the radii

    Raises:
        ValueError: If the number of blank fields in the random catalogue (nBlankRand) is zero, preventing division by zero.
    """

    radiiDeg = np.linspace(0, searchRadiusDeg, 20)   # 20 radii up to the search limit

    UobsByUrandArray = []
    FrsArray = []

    for rs in radiiDeg:

        nBlankReal = countBlanks(radCatCoords, optCatCoords, rs)
        nBlankRand = countBlanks(randRadCatCoords, optCatCoords, rs)

        if nBlankRand == 0:
            print("nBlankRand is zero — cannot divide by zero when computing Q0. Setting Q0 = 1.0")
            continue

        Frs = 1 - np.exp( -0.5 * (rs**2 / sigmaRadDeg**2))

        UobsByUrandArray.append(nBlankReal / nBlankRand)
        FrsArray.append(Frs)

    UobsByUrandArray = np.array(UobsByUrandArray)
    FrsArray = np.array(FrsArray)

    # Fit: Y = 1 - Q0 * X => (1 - Y) = Q0 * X
    Q0, _ = np.polyfit(FrsArray, 1 - UobsByUrandArray, 1)

    return Q0, radiiDeg, UobsByUrandArray, FrsArray

def getQ0SingleRs(radCatCoords, randRadCatCoords, optCatCoords, searchRadRsDeg, sigmaRadDeg):
    """
    Estimate Q0, the fraction of radio sources with real counterparts, using the blank fields method.

    Args:
        radCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the radio catalogue positions.
        randRadCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the randomised radio catalogue positions.
        optCatCoords (:obj:`~astropy.coordinates.SkyCoord`): SkyCoord object for the optical/IR catalogue positions.
        searchRadRsDeg (:obj:`float`): Search radius in deg within which counterparts are considered.
        sigmaRadDeg (:obj:`float`): Typical positional uncertainty (sigma) of radio sources in deg.

    Returns:
        float: Estimated Q0 value representing the fraction of radio sources with true counterparts.

    Raises:
        ValueError: If the number of blank fields in the random catalogue (nBlankRand) is zero, preventing division by zero.
    """

    nBlankReal = countBlanks(radCatCoords, optCatCoords, searchRadRsDeg)

    nBlankRand = countBlanks(randRadCatCoords, optCatCoords, searchRadRsDeg)

    if nBlankRand == 0:
        print("nBlankRand is zero — cannot divide by zero when computing Q0. Setting Q0 = 1.0")
        return 1.0

    Frs = 1 - np.exp( -0.5 * (searchRadRsDeg**2 / sigmaRadDeg**2))

    if Frs == 0:
        print("Frs is zero — cannot divide by zero when computing Q0. Setting Q0 = 1")
        return 1.0

    Q0 = (1 - nBlankReal/nBlankRand)/Frs

    print("\nQ0 = ", Q0)
    return Q0

def getQM(radCatCoords, nRadio, optCatCoords, optMagList, searchRadRmaxDegVal, nM, Q0):

    """
    Estimate q(m): the magnitude distribution of real counterparts to radio sources.

    Args:
        radCatCoords (SkyCoord): SkyCoord array of radio source positions.
        nRadio (int): Number of radio sources.
        optCatCoords (SkyCoord): SkyCoord array of optical catalogue positions.
        optMagList (array_like): List or array of optical magnitudes corresponding to optCatCoords.
        searchRadRmaxDeg (float): Maximum search radius in deg.
        nM (array_like): Surface density of background sources per magnitude bin, in sources/sq.deg.
        Q0 (float): Fraction of radio sources that have real optical counterparts (0 < Q0 < 1).

    Returns:
        np.ndarray: Estimated q(m), the probability distribution of magnitudes for real counterparts.
    """
    nMagBins = len(nM)

    # getting total(m)

    # Find all pairs within r_max
    _, idxOpt, _, _ = search_around_sky(radCatCoords, optCatCoords, searchRadRmaxDegVal * u.deg)

    if len(idxOpt) == 0:
        print("\nNo cross-matches...!")
        return None

    # idxOpt are indices of optical sources within r_max of any radio source
    matchedMags = optMagList[idxOpt]

    totalM, _ = getMagHist(optMag=matchedMags, nBins=nMagBins)

    areaSearchRadiusSqDeg = np.pi * searchRadRmaxDegVal**2  # in sq.deg.
    backgroundCounts = nM * nRadio * areaSearchRadiusSqDeg

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

def getCentreRadiusFromImagesTab(radioCatFilePath):
    """
    Retrieve the pointing centre coordinates and radius
    for a given radio catalogue from the 'images.fits' table.

    Args:
        radioCatFilePath (:obj:`str`): Full path to the radio source catalogue
            (e.g., the PyBDSF output file). The function uses the catalogue
            name to locate the corresponding image entry in the global
            `images.fits` table.

    Returns:
        tuple:
            - fieldRACentre (:obj:`float`): Right ascension of the field centre (degrees).
            - fieldDecCentre (:obj:`float`): Declination of the field centre (degrees).
            - bandRadiusDeg (:obj:`float`): Approximate primary beam radius (degrees),
              determined by the observing band (L, UHF, or S).
    """

    bandRadiusDict={'L': 0.8, 'UHF': 1.2, 'S': 0.6}

    catalogName = radioCatFilePath.split(os.path.sep)[-1]
    commonFilenamePart = catalogName.split('_srl')[0]

    globalImagesFileName = startup.config['productsDir']+os.path.sep+"images.fits"
    globalImagesTab = Table.read(globalImagesFileName)

    fieldMask = np.array([commonFilenamePart in p for p in globalImagesTab['path']])

    fieldRACentre = globalImagesTab['centre_RADeg'][fieldMask][0]
    fieldDecCentre = globalImagesTab['centre_decDeg'][fieldMask][0]

    bandName = globalImagesTab['band'][fieldMask][0]
    bandRadiusDeg = bandRadiusDict[bandName]

    return fieldRACentre, fieldDecCentre, bandRadiusDeg

def getCentreRadiusFromCatalog(radioCat, radRACol, radDecCol):
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

def getDECaLSSkyAreaSqDeg(decalsCat, healPixNSide=4096, healPixCountCol='nest4096'):
    """
    Calculates the effective sky area using a HEALPix-based method.

    Args:
        decalsCat (:obj:`~astropy.table.Table`): DECaLS catalogue in Astropy Table format.
        healPixNSide (int): HEALPix Nside resolution of the index column (default: 4096 for DECaLS DR10).
        healPixCountCol (str): Column name containing the HEALPix index in the NESTED scheme (e.g., 'nest4096').

    Returns:
        float:
            Total sky area in square degrees spanned by the catalogue based on the unique HEALPix footprint.
    """

    # Check for the HEALPix column
    if healPixCountCol not in decalsCat.colnames:
        print("ERROR: %s column not found in the table. Cannot calculate HEALPix area." %healPixCountCol)
        return None

    # Get all unique HEALPix indices (the survey footprint)
    uniqueHealpixIndices = np.unique(decalsCat[healPixCountCol])
    numUniquePixels = len(uniqueHealpixIndices)

    # tractor cat documentation says:
    # best4096: HEALPIX index (Nsides 4096, Nest scheme)
    hp = HEALPix(nside=healPixNSide, order='nested', frame='icrs')

    # Get the area of one pixel
    areaPerPixelSqDegVal = hp.pixel_area.to(u.deg**2).value

    # Total HEALPix Area
    totalAreaSqDeg = numUniquePixels * areaPerPixelSqDegVal

    return totalAreaSqDeg

def retrieveDECaLSDR10(centerRA, centerDec, radiusDeg):
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

    token = DLLogin()

    print("\nRetrieving DECaLS DR10 sources with RA_central=%.2f deg, Dec_central=%.2f deg, and radius=%.2f deg" \
          % (centerRA, centerDec, radiusDeg))
    zClusterCacheDir = os.path.join(os.environ['ZCLUSTER_CACHE'], "zCluster", "cache")
    resultRetrieve = retrievers.DL_DECaLSDR10RetrieverPhotoZ(centerRA, centerDec,
                                                             halfBoxSizeDeg = radiusDeg,
                                                             DR = None,
                                                             optionsDict={'altCacheDir': zClusterCacheDir,
                                                                          'token': token})

    if resultRetrieve is not None:
        decalsCat = Table(resultRetrieve)
    else:
        return None

    return decalsCat

def retrieveRubinDP1(centerRA, centerDec, radiusDeg):
    """  Retrieve Rubin DP1 sources within a circular region around a given sky position.

    Args:
        centerRA (:obj:`float`): Right ascension (RA) of the centre of the search region, in degrees.
        centerDec (:obj:`float`): Declination (Dec) of the centre of the search region, in degrees.
        radiusDeg (:obj:`float`): Search radius around the central position, in degrees.

    Returns:
        :obj:`astropy.table.Table`: Table of Rubin sources within the specified region.

    """

    print("\nRetrieving Rubin sources with RA_central=%.2f deg, Dec_central=%.2f deg, and radius=%.2f deg" \
          % (centerRA, centerDec, radiusDeg))
    zClusterCacheDir = os.environ['ZCLUSTER_CACHE']+os.path.sep+"zCluster"+os.path.sep+"cache"

    RubinURL = os.environ['RUBIN_URL']
    RubinToken = os.environ['RUBIN_TOKEN']

    tapSession = requests.Session()
    tapSession.headers['Authorization'] = 'Bearer %s' %RubinToken

    tapService = vo.dal.TAPService(RubinURL, session=tapSession)

    resultRetrieve = retrievers.RubinDP1Retriever(centerRA, centerDec,
                                                halfBoxSizeDeg = radiusDeg,
                                                optionsDict={
                                                    'altCacheDir': zClusterCacheDir,
                                                    'TAP': tapService})

    if resultRetrieve is not None:
        RubinCat = Table(resultRetrieve)
    else:
        return None

    return RubinCat

def getFRSymErr(rOffsetDeg, radioSource, opticalSource, radERACol, radEDecCol, optPosErrCol):

    """
    Calculate the probability distribution f(r) of offset r between radio and a potential counterpart.

    Args:
        rOffsetDeg (:obj:`float` or :obj:`np.ndarray`): Angular offset (separation) between radio and optical positions,
            in degrees.
        radioSource (:obj:`~astropy.table.Row`): A single row from the radio catalogue table.
        opticalSource (:obj:`~astropy.table.Row`): A single row from the optical catalogue table.
        radERACol (:obj:`str`): Key for error in right ascension of the radio source in `radioSource`.
        radEDecCol (:obj:`str`): Key for error in declination of the radio source in `radioSource`.
        optPosErrCol (:obj:`str`): Key for positional uncertainty of the optical source.

    Returns:
        float or np.ndarray:
            The value of the positional probability distribution f(r) for the given offset(s).

    """

    # sigmaPos as in Eqn.4 of McAlpine+2012
    sigmaRARad = radioSource[radERACol]
    sigmaDecRad = radioSource[radEDecCol]

    sigmaRad = np.sqrt(sigmaRARad**2 + sigmaDecRad**2)
    sigmaOpt = opticalSource[optPosErrCol]

    sigmaPos = np.sqrt(sigmaRad**2 + sigmaOpt**2)

    # There is a missing -r in Eqn.4 of McAlpine+2012
    probDistR = (1 / (2 * np.pi * sigmaPos**2)) * np.exp(-0.5 * (rOffsetDeg**2 / sigmaPos**2))

    return probDistR  # f(r)

def getFRRadioAsymErrOptSymErr(rOffsetDeg, radioSource, opticalSource, radRACol, radDecCol, radEMajCol, radEMinCol, radPACol, optRACol, optDecCol, optPosErrCol):

    """
    Calculate the probability distribution f(r) of offset r between radio and a potential counterpart,

    Args:
        rOffsetDeg (:obj:`float` or :obj:`np.ndarray`): Angular offset (separation) between radio and optical positions,
            in deg.
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
        sigmaAstArcsecVal (:obj:`float`, optional): Astrometric uncertainty between radio and optical surveys (default 0.6 arcsec).

    Returns:
        float or np.ndarray:
            The value of the positional probability distribution f(r) for the given offset(s).

    """

    deltaMaj = radioSource[radEMajCol]
    deltaMin = radioSource[radEMinCol]
    sigmaMajRad = deltaMaj/np.sqrt(4 * np.log(2)) # William et al. 2019
    sigmaMinRad = deltaMin/np.sqrt(4 * np.log(2))

    # getting vector joining radio to optical counterpart

    RARad = radioSource[radRACol]
    decRad = radioSource[radDecCol]
    RAOpt = opticalSource[optRACol]
    decOpt = opticalSource[optDecCol]

    dRA = (RAOpt - RARad) * np.cos(np.radians(decRad))
    dDec = decOpt - decRad

    # angle (in radians) of the vector from the radio source to the
    # optical candidate, measured from the RA (east) direction
    thetaDir = np.arctan2(dDec, dRA)

    positionAngle = np.radians(radioSource[radPACol])

    # angle between major axis and radio-optical vector
    thetaPADir = thetaDir - positionAngle

    sigmaRad = np.sqrt(
        (sigmaMajRad * np.cos(thetaPADir))**2 +
        (sigmaMinRad * np.sin(thetaPADir))**2
    )

    sigmaOpt = opticalSource[optPosErrCol] # sigma for optical

    sigmaPos = np.sqrt(sigmaRad**2 + sigmaOpt**2)

    probDistR = (1 / (2 * np.pi * sigmaPos**2)) * np.exp(-0.5 * (rOffsetDeg**2 / sigmaPos**2))

    return probDistR  # f(r)


def getFRAsymErr(rOffsetDeg, radioSource, opticalSource, radRACol, radDecCol, radEMajCol, radEMinCol, radPACol, optRACol, optDecCol, optPosErrCol, sigmaAstArcsecVal=0.6):

    """
    Calculate the probability distribution f(r) of offset r between radio and a potential counterpart, based on William et al. 2019

    Args:
        rOffsetDeg (:obj:`float` or :obj:`np.ndarray`): Angular offset (separation) between radio and optical positions,
            in deg.
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
        sigmaAstArcsecVal (:obj:`float`, optional): Astrometric uncertainty between radio and optical surveys (default 0.6 arcsec).

    Returns:
        float or np.ndarray:
            The value of the positional probability distribution f(r) for the given offset(s).

    """

    sigmaAstDegVal = sigmaAstArcsecVal/3600.0

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

    # angle (in radians) of the vector from the radio source to the
    # optical candidate, measured from the RA (east) direction
    thetaDir = np.arctan2(dDec, dRA)

    positionAngle = np.radians(radioSource[radPACol])

    # angle between major axis and radio-optical vector
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

    sigmaMaj = np.sqrt(sigmaMajRad**2 + sigmaMajOpt**2 + sigmaAstDegVal**2)
    sigmaMin = np.sqrt(sigmaMinRad**2 + sigmaMinOpt**2 + sigmaAstDegVal**2)
    sigmaDir = np.sqrt(sigmaDirRad**2 + sigmaDirOpt**2 + sigmaAstDegVal**2)

    probDistR = (1 / (2 * np.pi * sigmaMaj * sigmaMin)) * np.exp(-0.5 * (rOffsetDeg**2 / sigmaDir**2))

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

def computeLR(radioCat, opticalCat, searchRadiusDegVal, optMagCol, magBins, qMList, nMList, radRACol, radDecCol, radERACol, radEDecCol, optRACol, optDecCol, radEMajCol, radEMinCol, radPACol, optPosErrCol, dofRSymErr=True, sigmaAstArcsecVal=0.6):
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
        searchRadiusDegVal (:obj:`float`): Search radius around radio sources in degrees.
        optMagCol (:obj:`str`): Column name for optical magnitudes in opticalCat.
        magBins (:obj:`array_like`): Bin edges used for q(m) and n(m) calculation.
        qMList (:obj:`array_like`): q(m) values for magnitude bins.
        nMList (:obj:`array_like`): n(m) values for magnitude bins.
        radRACol (:obj:`str`): RA column name in radioCat.
        radDecCol (:obj:`str`): Dec column name in radioCat.
        radERACol (:obj:`str`): Radio catalog RA positional error column name.
        radEDecCol (:obj:`str`): Radio catalog Dec positional error column name.
        optRACol (:obj:`str`): RA column name in opticalCat.
        optDecCol (:obj:`str`): Dec column name in opticalCat.
        radEMajCol (:obj:`str`): Major axis error column name in radioCat.
        radEMinCol (:obj:`str`): Minor axis error column name in radioCat.
        radPACol (:obj:`str`): Position angle column name in radioCat.
        optPosErrCol (:obj:`str`): Positional error column name in opticalCat.
        dofRSymErr (:obj:`bool`, optional): Uses symmetric error approach to calculate f(r). Default is True.
        sigmaAstArcsecVal (:obj:`float`, optional): Astrometric uncertainty between radio and optical surveys (default 0.6 arcsec).

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
    idxRadio, idxOpt, _, _ = search_around_sky(radCatCoords, optCatCoords, searchRadiusDegVal * u.deg)

    if len(idxOpt) == 0:
        print("\nNo optical candidates within the search radius...!")
        return None

    rowsRadio = []
    rowsOptical = []
    radOptSeparationDeg = []
    fRVals = []
    qMVals = []
    nMVals = []
    LRVals = []

    for rIdx, oIdx in zip(idxRadio, idxOpt):
        radioSource = radioCat[rIdx]
        opticalSource = opticalCat[oIdx]

        radioCoord = SkyCoord(ra=radRAList[rIdx] * u.deg, dec=radDecList[rIdx] * u.deg)
        opticalCoord = SkyCoord(ra=optRAList[oIdx] * u.deg, dec=optDecList[oIdx] * u.deg)

        radOptOffsetDeg = radioCoord.separation(opticalCoord).deg

        if dofRSymErr is True:
            fRPair = getFRRadioAsymErrOptSymErr(radOptOffsetDeg, radioSource, opticalSource, radRACol, radDecCol, radEMajCol, radEMinCol, radPACol, optRACol, optDecCol, optPosErrCol)
            #fRPair = getFRSymErr(radOptOffsetDeg, radioSource, opticalSource, radERACol, radEDecCol, optPosErrCol)
        else:
            fRPair = getFRAsymErr(radOptOffsetDeg, radioSource, opticalSource, radRACol, radDecCol, radEMajCol, radEMinCol, radPACol, optRACol, optDecCol, optPosErrCol, sigmaAstArcsecVal=sigmaAstArcsecVal)

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
        radOptSeparationDeg.append(radOptOffsetDeg)
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

    radOptMergedTab['rad_opt_sep_deg'] = radOptSeparationDeg
    radOptMergedTab['f_r'] = fRVals
    radOptMergedTab['f_r'] = fRVals
    radOptMergedTab['q_m'] = qMVals
    radOptMergedTab['n_m'] = nMVals
    radOptMergedTab['LR'] = LRVals

    return radOptMergedTab

def xmatchRadioOptical(radioCatFilePath, radioBand, xmatchDirPath, optSurvey, optMagCol, searchRadiusArcsec, radRACol, radDecCol, radERACol, radEDecCol, radEMajCol, radEMinCol, radPACol, outSubscript, optPosErrAsecValue, nMagBins=15, beamSizeArcsecValue=6.0, nAreaTimeOptFetch=2.0, skipIfExists=True, dofRSymErr=True):
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
        optSurvey (:obj:`str`): Optical survey name (e.g., 'DECaLSDR10', 'RubinDP1').
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
        nAreaTimeOptFetch (:obj:`float`, optional): Factor of extra optical area to be fetched. Default is 3.0.
        skipIfExists (:obj:`bool`, optional): Skip processing if output files exist. Default is True.
        dofRSymErr (:obj:`bool`, optional): Uses symmetric error approach to calculate f(r). Default is True.

    Returns:
        astropy.table.Table: Table of best crossmatched sources with LR values and Reliability and Completeness in meta.
    """

    # Quantities to degree from arcsec.
    beamSizeDegValue = beamSizeArcsecValue/3600.0
    searchRadiusDegVal = searchRadiusArcsec/3600.0
    optPosErrValueDeg = optPosErrAsecValue/3600.0

    radCatName = radioCatFilePath.split(os.path.sep)[-1]
    captureBlockId = radCatName.split('_')[3]
    targetName = (radCatName.split('_1024ch_')[1]).split('_srl_')[0]
    xmatchIndividualDirPath = os.path.join(xmatchDirPath, 'xmatch_%s' %outSubscript)
    xmatchTabName = xmatchIndividualDirPath+os.path.sep+"xmatchtable_%s" %outSubscript+".fits"
    xmatchBestMatchTabName = xmatchIndividualDirPath+os.path.sep+"xmatchtable_bestmatches_%s" %outSubscript+".fits"

    if skipIfExists is True:
        doesItExist = os.path.exists(xmatchBestMatchTabName) and os.path.exists(xmatchTabName)
        if doesItExist:
            xmatchTable = Table.read(xmatchBestMatchTabName)
            return xmatchTable

    # checking if this catalog is listed as having no optical sources in the area
    noOptSourcesFilename = xmatchDirPath+os.path.sep+'no_%s_sources.txt' %optSurvey

    if os.path.exists(noOptSourcesFilename):
        with open(noOptSourcesFilename, 'r', encoding="utf-8") as infile:
            noOptSourcesCatNames = [line.strip() for line in infile]
        if radCatName in noOptSourcesCatNames:
            return None

    # checking if this catalog is listed as having no optical counterparts in the optical
    noOptCounterpartsFilename = xmatchDirPath+os.path.sep+'no_counterparts_%s_%sband_%sasec.txt' \
                           %(optSurvey, optMagCol, str(searchRadiusArcsec).replace(".","p"))

    if os.path.exists(noOptCounterpartsFilename):
        with open(noOptCounterpartsFilename, 'r', encoding="utf-8") as infile:
            noOptCounterpartsCatNames = [line.strip() for line in infile]
        if radCatName in noOptCounterpartsCatNames:
            return None

    print("\n" + "-" * 100)
    print("║ Radio catalogue: %s ║" %radCatName)
    print("-" * 100 + "\n")

    radioSources = Table.read(radioCatFilePath, format='fits', hdu=1)

    nRadio = len(radioSources)
    print("\nNumber of radio sources: %d" %nRadio)

    radColumns = [radRACol, radDecCol, radERACol, radEDecCol, radEMajCol, radEMinCol, radPACol]

    if not all(col in radioSources.colnames for col in radColumns):
        print("\nERROR: Required columns of radio catalogue is not well set!")
        return None

    if any(radioSources[radRACol] < 0.0):
        # This pybdsf catalog has -180 to 180 wrapping. Need to change to 360 wrapping
        radioSources = catalogs.fixRA(radioSources, raCol=radRACol, wrapAngle=360)

    radRAValDegList = _getUnitlessValues(radioSources[radRACol])
    radDecValDegList = _getUnitlessValues(radioSources[radDecCol])

    radioSourcesCoords = SkyCoord(ra=radRAValDegList * u.deg, dec=radDecValDegList * u.deg, frame='icrs')

    centerRA, centerDec, radiusDeg = getCentreRadiusFromImagesTab(radioCatFilePath)

    sigmaRadPosDeg = np.sqrt(radioSources[radERACol]**2 + radioSources[radEDecCol]**2)
    sigmaRadPosMeanDeg = np.mean(sigmaRadPosDeg)

    # Collecting optical sources
    optRACol, optDecCol = 'RADeg', 'decDeg'
    optPosErrCol = 'pos_err_deg'

    optCatFileName = xmatchIndividualDirPath+os.path.sep+'%s_sources_%s.fits' %(optSurvey, outSubscript)

    # Radius for which optical catalogue is to be fetched
    radiusDegOptFetch = radiusDeg*nAreaTimeOptFetch

    if optSurvey == 'DECaLSDR10':
        opticalSourcesRaw = retrieveDECaLSDR10(centerRA, centerDec, radiusDegOptFetch)
    elif optSurvey == 'RubinDP1':
        opticalSourcesRaw = retrieveRubinDP1(centerRA, centerDec, radiusDegOptFetch)
    else:
        print("\nERROR: Berk is currently setup for DECaLSDR10 and RubinDP1 only.")
        return None

    if opticalSourcesRaw is None:
        print("\nRetrieval process unsuccessfull.")
        catalogs.listCatalogInFile(radCatName, noOptSourcesFilename)
        return None

    # filtering optical catalogue
    opticalSources = filterBadRows(opticalSourcesRaw, [optMagCol])

    if len(opticalSources) == 0:
        print("\n%s: No optical sources with reliable %s-magnitude found in %s database...!" %(radCatName, optMagCol, optSurvey))
        catalogs.listCatalogInFile(radCatName, noOptSourcesFilename)
        return None

    # getting sky area covered by the optical catalogue (DECaLS for now!) #TODO

    opticalSkyAreaSqDeg = getDECaLSSkyAreaSqDeg(opticalSources, healPixNSide=4096, healPixCountCol='nest4096')

    # writing and plotting the optical sources in the field
    os.makedirs(xmatchIndividualDirPath, exist_ok = True)
    print("\nNumber of %s sources with reliable %s-magnitude in the sky region: %d" % (optSurvey, optMagCol, len(opticalSources)))

    opticalSources.meta['AREASQDEG'] = opticalSkyAreaSqDeg
    opticalSources.write(optCatFileName, format='fits', overwrite=True)

    RADecplotOutName = "%s/RadOptSkyPlot_%s.png" %(xmatchIndividualDirPath, outSubscript)

    plt.figure(figsize=(6, 6))
    plt.scatter(opticalSources[optRACol], opticalSources[optDecCol], s=1,
                c='#8AD5F1', label='%s (N=%d)'%(optSurvey, len(opticalSources)))
    plt.scatter(radioSources[radRACol], radioSources[radDecCol], s=2, c='#87340D',
                label='MeerKAT (N=%d)' %len(radioSources))
    plt.title(radCatName)
    plt.xlabel("RA (deg; J2000)")
    plt.ylabel("Dec (deg; J2000)")
    plt.legend(loc="lower left", scatterpoints=1, fontsize=10)
    plt.savefig(RADecplotOutName, dpi=300, bbox_inches = 'tight')
    plt.close()

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
    randomRadioSources.write(xmatchIndividualDirPath+os.path.sep+'Randoms_%s.fits' %outSubscript, format='fits', overwrite=True)

    randRadRAValDegList = _getUnitlessValues(randomRadioSources[radRACol])
    randRadDecValDegList = _getUnitlessValues(randomRadioSources[radDecCol])

    randomRadioSourcesCoords = SkyCoord(ra= randRadRAValDegList * u.deg, dec= randRadDecValDegList * u.deg, frame='icrs')

    #Q0 = getQ0SingleRs(radCatCoords=radioSourcesCoords, randRadCatCoords=randomRadioSourcesCoords, optCatCoords=optSourcesCoords, #searchRadRsDeg=beamSizeDegValue, sigmaRadDeg=sigmaRadPosMeanDeg)

    Q0, radiiDeg, UobsByUrandArray, FrsArray = getQ0FitForRs(radCatCoords=radioSourcesCoords, randRadCatCoords=randomRadioSourcesCoords, optCatCoords=optSourcesCoords,
                       searchRadiusDeg=searchRadiusDegVal, sigmaRadDeg=sigmaRadPosMeanDeg)
    np.savetxt(xmatchIndividualDirPath+os.path.sep+'Q0_%s.txt' %outSubscript, [Q0], fmt='%f')

    # Finding n(m)

    nM = getNM(optMag=optMagList, nBins=nMagBins, areaSqDeg=opticalSkyAreaSqDeg)

    # Finding q(m)

    qM = getQM(radCatCoords=radioSourcesCoords, nRadio=nRadio, optCatCoords=optSourcesCoords, optMagList=optMagList, searchRadRmaxDegVal=beamSizeDegValue, nM=nM, Q0=Q0)
    if qM is None:
        catalogs.listCatalogInFile(radCatName, noOptCounterpartsFilename)
        return None

    print("\nComputing LR ...")

    xmatchTable = computeLR(radioCat=radioSources, opticalCat=opticalSources, searchRadiusDegVal=searchRadiusDegVal, optMagCol=optMagCol, magBins=optMagBins, qMList=qM, nMList=nM, radRACol=radRACol, radDecCol=radDecCol, radERACol=radERACol, radEDecCol=radEDecCol, optRACol=optRACol, optDecCol=optDecCol, radEMajCol=radEMajCol, radEMinCol=radEMinCol, radPACol=radPACol, optPosErrCol=optPosErrCol, dofRSymErr=dofRSymErr)

    if xmatchTable is None:
        print("\n%s: No cross-matched objects...!" %radCatName)
        catalogs.listCatalogInFile(radCatName, noOptCounterpartsFilename)
        return None

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
        catalogs.listCatalogInFile(radCatName, noOptCounterpartsFilename)
        return None


    # Filters the full xmatch table to keep only the matches with the maximum 'LR' value for each radio source
    xmatchLRThresholdTable.sort(['Source_id_rad', 'LR'])
    xmatchLRThresholdTable.reverse()
    groupedxmatchLRThresholdTable = xmatchLRThresholdTable.group_by('Source_id_rad')
    xmatchBestMatchTable = groupedxmatchLRThresholdTable.groups.aggregate(lambda rows: rows[0])

    xmatchBestMatchTable.meta['OPT_SUR']='%s' %optSurvey
    xmatchBestMatchTable.meta['SEAR_RAD']='%f arcsec' %searchRadiusArcsec
    xmatchBestMatchTable.meta['LR_THR']=CRBalanceLRThreshold
    xmatchBestMatchTable.meta['REL']=CRBalanceRel
    xmatchBestMatchTable.meta['COMP']=CRBalanceComp

    # Saving the results

    # Plotting optical and radio sources
    RADecplotOutName = "%s/RadOptSkyPlot_%s.png" %(xmatchIndividualDirPath, outSubscript)
    if os.path.exists(RADecplotOutName):
        os.remove(RADecplotOutName)
    plt.figure(figsize=(6, 6))
    plt.scatter(opticalSources[optRACol], opticalSources[optDecCol], s=1,
                c='#8AD5F1', label='%s (N=%d)'%(optSurvey, len(opticalSources)))
    plt.scatter(radioSources[radRACol], radioSources[radDecCol], s=2, c='#87340D',
                label='MeerKAT (N=%d)' %len(radioSources))
    plt.scatter(xmatchBestMatchTable[radRACol+'_rad'], xmatchBestMatchTable[radDecCol+'_rad'],
                marker='o', facecolor='None', linewidth=0.5, s=15,
                edgecolor='#06471D',
                label='MeerKATx%s (Best matches; N=%d)' %(optSurvey, len(xmatchBestMatchTable)))
    plt.title("%s\nSearch radius = %0.1f asec, %s band, Q0=%0.2f" \
                % (radCatName, searchRadiusArcsec, optMagCol, Q0))
    plt.xlabel("RA (deg; J2000)")
    plt.ylabel("Dec (deg; J2000)")
    plt.legend(loc="lower left", scatterpoints=1, fontsize=10)
    plt.savefig(RADecplotOutName, dpi=300, bbox_inches = 'tight')
    plt.close()
    print("\nPlotted sky coverage of radio and optical sources!")

    # Plotting optical and radio sources
    Q0RsplotOutName = "%s/Q0_rs_%s.png" %(xmatchIndividualDirPath, outSubscript)
    if os.path.exists(Q0RsplotOutName):
        os.remove(Q0RsplotOutName)
    plt.figure(figsize=(6, 6))
    plt.scatter(radiiDeg*3600, UobsByUrandArray)
    yFit = 1 - Q0 * FrsArray
    plt.plot(radiiDeg*3600, yFit, label='Fit: $1 - Q_0 F(r)$')
    plt.xlabel("radius (arcsec)")
    plt.ylabel(r"$1-Q_0 F(r)$")
    plt.legend(loc="lower left", scatterpoints=1, fontsize=10)
    plt.savefig(Q0RsplotOutName, dpi=300, bbox_inches = 'tight')
    plt.close()
    print("\nPlotted Q0 vs rs plot!")

    # Plotting q(m)/n(m) and n(m) as a function of magnitude

    qmNmPlotOutName = "%s/QmNmPlot_%s.png" % (xmatchIndividualDirPath, outSubscript)
    plt.figure(figsize=(8, 5))

    # Avoid division by zero
    qMnMRatio = np.full_like(qM, np.nan, dtype=float)
    mask = nM != 0
    qMnMRatio[mask] = qM[mask] / nM[mask]

    # Left axis for q(m)/n(m)
    fig, ax1 = plt.subplots(figsize=(8, 5))
    color1 = 'tab:blue'
    ax1.plot(optMagBinCenters, qMnMRatio, color=color1, label=r'$q(m)/n(m)$')
    ax1.set_xlabel(optMagCol + '-band magnitude')
    ax1.set_ylabel(r'$q(m)/n(m)$', color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)

    # Right axis for n(m)
    ax2 = ax1.twinx()
    color2 = 'tab:red'
    ax2.plot(optMagBinCenters, nM, color=color2, label=r'$n(m)$')
    ax2.set_ylabel(r'$n(m)$', color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)

    # Optional: add legends
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')

    plt.tight_layout()
    plt.savefig(qmNmPlotOutName, dpi=300, bbox_inches='tight')
    plt.close()
    print("\nPlotted q(m)/n(m) and n(m) vs magnitude!")

    # Plotting LR and Rel
    LRRelPlotOutName = "%s/LRRelMagPlot_%s.png" %(xmatchIndividualDirPath, outSubscript)
    plt.figure(figsize=(8, 5))
    plt.plot(LRThreshValues, completenessValues, label='Completeness', color='blue')
    plt.plot(LRThreshValues, reliabilityValues, label='Reliability', color='green')
    plt.gca().axvline(x=CRBalanceLRThreshold, linestyle='dashed', color='k', label='Choosen LR Threshold = %.2f' %CRBalanceLRThreshold)
    plt.xlabel('LR Threshold')
    plt.ylabel('Completeness or Reliability')
    plt.legend()
    plt.tight_layout()
    plt.savefig(LRRelPlotOutName , dpi=300, bbox_inches = 'tight')
    plt.close()
    print("\nPlotted Reliability/Completeness vs LR_threshold!")

    # Saving crossmatch tables

    xmatchTable.write(xmatchTabName, format='fits', overwrite=True)
    print("\nWrote full cross-matched table %s." %xmatchTabName)

    xmatchBestMatchTable.write(xmatchBestMatchTabName, format='fits', overwrite=True)
    print("\nWrote cross-matched table with Reliability %0.2f and completeness %0.2f: %s." % (CRBalanceRel, CRBalanceComp, xmatchBestMatchTabName))

    return xmatchBestMatchTable
