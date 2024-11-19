import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, hstack
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
import os
import argparse
from zCluster import retrievers
from . import startup, tasks

from matplotlib import font_manager as fm
this_platform = os.environ['THIS_PLATFORM']
if(this_platform == 'hp455'): #krishna's laptop
    fpath = "/home/krishna/Dropbox/fonts/cmunss.ttf"
    prop = fm.FontProperties(fname=fpath,size=12,math_fontfamily='stixsans')
else:
    prop = fm.FontProperties(size=12,math_fontfamily='stixsans')

def _write_no_optical_catFile(catalogName, outfilename):

	if os.path.exists(outfilename):
		with open(outfilename, 'a') as file:
			file.write("%s\n" %(catalogName))
	else:	
		with open(outfilename, 'w') as file:	
			file.write("%s\n" %(catalogName))
			
	return 0


def get_mag_hist(opt_mag, nbins):
	opt_mag_hist, opt_mag_bins = np.histogram(opt_mag, bins=nbins)
	return opt_mag_hist, opt_mag_bins

def get_n_m(opt_mag, nbins, area_asec):
	opt_mag_hist = get_mag_hist(opt_mag, nbins)[0]
	n_m = opt_mag_hist/area_asec.value # TODO : check if it is whole area of survey or area of the overlapping regions
	return n_m

def get_matches(radio_cat_table, optical_cat_table, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol):

	radio_coords = SkyCoord(ra=radio_cat_table[radRACol].value * u.deg, dec=radio_cat_table[radDecCol].value * u.deg)
	optical_coords = SkyCoord(ra=optical_cat_table[optRACol].value * u.deg, dec=optical_cat_table[optDecCol].value * u.deg)

	radio_cartesian = np.vstack(radio_coords.cartesian.xyz).T  
	optical_cartesian = np.vstack(optical_coords.cartesian.xyz).T

	optical_tree = cKDTree(optical_cartesian)
	search_radius_rad = search_radius_asec.to(u.rad).value 

	matches = optical_tree.query_ball_point(radio_cartesian, r=search_radius_rad)

	return matches

def get_separation(radio_cat_table, optical_cat_table, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol):

	radio_coords = SkyCoord(ra=radio_cat_table[radRACol].value * u.deg, dec=radio_cat_table[radDecCol].value * u.deg)
	optical_coords = SkyCoord(ra=optical_cat_table[optRACol].value * u.deg, dec=optical_cat_table[optDecCol].value * u.deg)
	matches = get_matches(radio_cat_table, optical_cat_table, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol)

	separation_list = []
	no_match_radio_idx = [] # List of radio sources with no match
	match_radio_idx = [] # List of radio sources with atleast one match
	match_optical_idx = [] # List of ALL optical sources that are matched to atleast one radio source

	for radio_idx, optical_indices in enumerate(matches):
		if len(optical_indices) == 0: # radio sources with no optical matches within the search radius
			no_match_radio_idx.append(radio_idx)
			continue 
		match_radio_idx.append(radio_idx)
		match_optical_idx.extend(optical_indices)
		radio_coord = radio_coords[radio_idx]
		optical_matches_coords = optical_coords[optical_indices]
		separations = radio_coord.separation(optical_matches_coords)
		for opt_idx, separation in zip(optical_indices, separations):
			separation_list.append((radio_idx, opt_idx, separation.arcsec))
		    
	return separation_list, match_radio_idx, match_optical_idx, no_match_radio_idx

 
def combined_error(sigma_rad, sigma_opt):
	return np.sqrt(sigma_rad**2 + sigma_opt**2)

def f_r(radio_coord, sigma_rad, optical_coord, sigma_opt):
	r = radio_coord.separation(optical_coord).value
	sigma_pos = combined_error(sigma_rad, sigma_opt)

	return (1 / (2 * np.pi * sigma_pos**2)) * np.exp(-0.5 * (r**2 / sigma_pos**2))
    
    
def random_points_in_circle(center_ra, center_dec, radius_deg, n_points):

	center_coord = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg, frame='icrs')

	rand_r = np.sqrt(np.random.uniform(0, radius_deg**2, n_points))  # Radius
	rand_theta = np.random.uniform(0, 2 * np.pi, n_points)  # Angle

	rand_ra = center_ra + rand_r * np.cos(rand_theta) / np.cos(center_dec * u.deg.to(u.rad))
	rand_dec = center_dec + rand_r * np.sin(rand_theta)

	coords = SkyCoord(ra=rand_ra * u.deg, dec=rand_dec * u.deg, frame='icrs')
	distances = center_coord.separation(coords)
	within_circle = distances < radius_deg * u.deg
	return rand_ra[within_circle], rand_dec[within_circle]

def make_random_cat(center_ra, center_dec, radius_deg, n_random_points, radRACol, radDecCol):
	rand_ra, rand_dec = random_points_in_circle(center_ra, center_dec, radius_deg, n_random_points)
	random_cat = Table([rand_ra, rand_dec], names=(radRACol, radDecCol))
	return random_cat

def get_Uobs_Urand_ratio(radio_sources, rand_radio_sources, optical_sources, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol):
	Uobs_r = len(get_separation(radio_sources, optical_sources, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol)[3])
	Urand_r = len(get_separation(rand_radio_sources, optical_sources, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol)[3])
	return Uobs_r/Urand_r

def Uobs_Urand_ratio_model(r , Q0, sigma):
	return 1. - (Q0*(1. - np.exp(-1*np.power(r,2.)/2.*np.power(sigma,2.))))


def get_q0(radio_sources, rand_radio_sources, optical_sources, beam_size, radRACol, radDecCol, optRACol, optDecCol):
	radii = np.arange( 1., beam_size, 0.5 )
	Uobs_Urand_ratios = []
	for radius in radii:
		search_radius_asec = float(radius)*u.arcsec
		Uobs_Urand_ratios.append(get_Uobs_Urand_ratio(radio_sources, rand_radio_sources, optical_sources, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol))
	params_fit, params_cov = curve_fit(Uobs_Urand_ratio_model, radii, Uobs_Urand_ratios)
	Q0, sig = params_fit
	Q0_err, sig_err = np.sqrt(np.diag(params_cov))
	return Q0, Q0_err

def get_q_m(radio_sources, rand_radio_sources, optical_sources, n_m, search_radius_asec, opt_mag_col, Q0, n_mag_bins, radRACol, radDecCol, optRACol, optDecCol):
	N_radio = len(radio_sources)
	separation_list, match_radio_idx, match_optical_idx, no_match_radio_idx = get_separation(radio_sources, optical_sources, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol)
	opt_match_mags = optical_sources[opt_mag_col][match_optical_idx]
	mag_bins = get_mag_hist(optical_sources[opt_mag_col], n_mag_bins)[1]
	total_m = np.histogram(opt_match_mags, bins=mag_bins)[0]
	real_m = total_m - (n_m * N_radio * np.pi * search_radius_asec.value**2)
	q_m = (real_m*Q0)/np.sum(real_m)
	return q_m
    
def get_qm_nm(mag, mag_bins, qm_nm_values):
	bin_index = np.digitize(mag, mag_bins, right=True) - 1
	if 0 <= bin_index < len(qm_nm_values):
		return qm_nm_values[bin_index]
	else:
		#print("Issue with getting qm_nm") 
		return np.nan  


def make_mag_bins(magnitudes, nbins):
	binwidth = (np.max(magnitudes)-np.min(magnitudes))/nbins
	mag_bins = np.arange( np.floor(np.min(magnitudes)), np.ceil(np.max(magnitudes)), binwidth )
	if np.ceil(np.max(magnitudes)) > np.max( mag_bins ):
		mag_bins = np.append( mag_bins, np.max(mag_bins)+0.4 )
	if np.floor(np.min(magnitudes)) < np.min( mag_bins ):
		mag_bins = np.insert( mag_bins, 0, np.min( mag_bins )-0.4 )
	return( mag_bins )
    
def add_postfix_to_columns(table, postfix):
	new_names = {colname: colname + postfix for colname in table.colnames}
	table.rename_columns(list(new_names.keys()), list(new_names.values()))
	return table
	
def get_centre_radius(radio_cat, radRACol, radDecCol):
	ra_list=list(radio_cat[radRACol])
	dec_list=list(radio_cat[radDecCol])
	coords = SkyCoord(ra=ra_list*u.degree, dec=dec_list*u.degree, frame='icrs')
	
	min_ra, max_ra = np.min(ra_list), np.max(ra_list)
	min_dec, max_dec = np.min(dec_list), np.max(dec_list)
	
	rad_ra = (max_ra - min_ra) / 2
	rad_dec = (max_dec - min_dec) / 2

	center_ra = (min_ra + max_ra) / 2
	center_dec = (min_dec + max_dec) / 2

	radius_deg_max = max(rad_ra, rad_dec)
	
	return center_ra, center_dec, rad_ra, rad_dec, radius_deg_max
	
def retrieve_DECaLS(center_ra, center_dec, radius_deg, DR='DR10'):
	print("\nRetrieving DECaLS %s sources with RAcen=%f, Deccen=%f, and radius=%f deg" %(DR, center_ra, center_dec, radius_deg)) 
	if(DR == 'DR10'):
		DECaLS_cat = Table(retrievers.DL_DECaLSDR10Retriever(center_ra, center_dec, halfBoxSizeDeg = radius_deg, DR = None))
	elif(DR == 'DR8'):
		DECaLS_cat = Table(retrievers.DL_DECaLSDR8Retriever(center_ra, center_dec, halfBoxSizeDeg = radius_deg, DR = None))
	return DECaLS_cat
	
def compute_LR_Rel(radio_sources, optical_sources, search_radius_asec, opt_mag_col, mag_bins, Q0, qm_nm_values, radRACol, radDecCol, optRACol, optDecCol, eRadRACol, eRadDecCol, optPosErrCol):

	matches = get_matches(radio_sources, optical_sources, search_radius_asec, radRACol, radDecCol, optRACol, optDecCol)
	
	LR_full = []
	radio_idx_list = []
	optical_indices_list = []

	for radio_idx, optical_indices in enumerate(matches):
		radio_idx_list.append(radio_idx)
		optical_indices_list.append(optical_indices)
		n_cand = len(optical_indices)
		LR_radio_i = []
		if n_cand > 0:
			for optical_idx in optical_indices:
				
            
				radio_coord = SkyCoord(ra=radio_sources[radRACol][radio_idx] * u.deg, dec=radio_sources[radDecCol][radio_idx] * u.deg)
				optical_coord = SkyCoord(ra=optical_sources[optRACol][optical_idx] * u.deg, dec=optical_sources[optDecCol][optical_idx] * u.deg)

				sigma_radio = np.sqrt(radio_sources[eRadRACol][radio_idx]**2 + radio_sources[eRadDecCol][radio_idx]**2)
				sigma_optical = optical_sources[optPosErrCol][optical_idx]

				F_R = f_r(radio_coord=radio_coord, sigma_rad=sigma_radio, optical_coord=optical_coord, sigma_opt=sigma_optical)

				QM_NM = get_qm_nm(optical_sources[opt_mag_col][optical_idx], mag_bins, qm_nm_values)

				LR = QM_NM * F_R

				LR_radio_i.append(LR)

		LR_full.append(LR_radio_i)
   
	# Computing reliability

	Rel_full = []
	for radio_idx, LRs in enumerate(LR_full):
		sum_Li = np.sum(LRs)
		Rel_radio_i = []
		for LR in LRs:
			Rel = LR/(sum_Li + (1-Q0))
			Rel_radio_i.append(Rel)
		Rel_full.append(Rel_radio_i)
	
	# Creating table of all matches within search radius with LR and Reliability

	xmatch_optical_idx_list = []
	xmatch_radio_idx_list = []
	xmatch_LR_list = []
	xmatch_Rel_list = []

	for radio_idx, optical_indices in enumerate(matches):
		for i, optical_idx in enumerate(optical_indices):
			xmatch_optical_idx_list.append(optical_idx)
			xmatch_radio_idx_list.append(radio_idx)
			xmatch_LR_list.append(LR_full[radio_idx][i])
			xmatch_Rel_list.append(Rel_full[radio_idx][i])
			
	radio_match_rows = radio_sources[xmatch_radio_idx_list] 
	optical_match_rows = optical_sources[xmatch_optical_idx_list] 
	
	radio_match_rows = add_postfix_to_columns(radio_match_rows, '_rad')
	optical_match_rows = add_postfix_to_columns(optical_match_rows, '_opt')


	# creating the cross-matched table
	xmatch_table = hstack([radio_match_rows, optical_match_rows])  # Stack radio and optical data horizontally
	xmatch_table['LRs'] = xmatch_LR_list  
	xmatch_table['Rels'] = xmatch_Rel_list  
	
	return xmatch_table
	
def xmatchOpt(catalog, opt_survey, opt_survey_dr, opt_mag_col, search_radius_asec, outPath, makePlots, radRACol, radDecCol, eRadRACol, eRadDecCol, outSubscript):

	catalogName = catalog.split(os.path.sep)[-1].replace(".fits","")

	search_radius_asec = search_radius_asec * u.arcsec
	n_mag_bins = 15
	beam_size = 7.0
	opt_pos_err = (0.2*u.arcsec).to(u.deg).value
	
	radio_sources = Table.read(catalog, format='fits', hdu=1)
	
	rad_pos_columns = [radRACol, radDecCol, eRadRACol, eRadDecCol]
	if not all(col in radio_sources.colnames for col in rad_pos_columns):
    		print("RA and Dec columns of radio catalogue is not well set. Exiting!")
    		return
	
	if any(radio_sources[radRACol] < 0.0):
            	print("\nCrossmatching - Fixing RA for %s" %catalogName)
            	radio_sources = tasks.fixRA(radio_sources, racol=radRACol, wrap_angle=360)
            	
	center_ra, center_dec, rad_ra, rad_dec, radius_deg = get_centre_radius(radio_sources, radRACol, radDecCol)
	N_radio = len(radio_sources)
	sky_area_sqdeg = np.pi*rad_ra*rad_dec*u.deg**2
	sky_area_sqasec = sky_area_sqdeg.to(u.arcsec**2)
	
	# Collecting optical sources
	if(opt_survey == 'DECaLS'):
		optRACol, optDecCol = 'RADeg', 'decDeg'
		optPosErrCol = 'pos_err'
		optical_sources = retrieve_DECaLS(center_ra, center_dec, radius_deg, DR=opt_survey_dr) 
		if(len(optical_sources) == 0):
			raise RuntimeError("\n%s: No optical sources found in %s database...!" %(catalogName, opt_survey))
		print("\n%d %s sources found in this region..." %(len(optical_sources), opt_survey))
	else:
		raise RuntimeError("\n%s: Optical survey not mentioned....!" %catalogName)
		
	if optPosErrCol not in optical_sources.colnames:
	
		optical_sources[optPosErrCol] = opt_pos_err

	rand_radio_sources=make_random_cat(center_ra, center_dec, radius_deg, N_radio, radRACol, radDecCol)

	print("\nFinding Q0...")
	Q0, Q0_err = get_q0(radio_sources, rand_radio_sources, optical_sources, beam_size, radRACol, radDecCol, optRACol, optDecCol)
	print("Q0=", Q0)
	n_m = get_n_m(opt_mag=list(optical_sources[opt_mag_col]), nbins=n_mag_bins, area_asec=sky_area_sqasec)
	q_m = get_q_m(radio_sources, rand_radio_sources, optical_sources, n_m, search_radius_asec, opt_mag_col, Q0, n_mag_bins, radRACol, radDecCol, optRACol, optDecCol)
	qm_nm_values = q_m/n_m
	mag_bins = make_mag_bins(optical_sources[opt_mag_col], n_mag_bins)

	print("\nComputing LR and Reliabilities...")
	xmatch_table = compute_LR_Rel(radio_sources, optical_sources, search_radius_asec, opt_mag_col, mag_bins, Q0, qm_nm_values, radRACol, radDecCol, optRACol, optDecCol, eRadRACol, eRadDecCol, optPosErrCol)
	
	# computing cross-match separation
	
	coords_rad = SkyCoord(ra=xmatch_table[radRACol+'_rad'].value*u.deg, dec=xmatch_table[radDecCol+'_rad'].value*u.deg)
	coords_opt = SkyCoord(ra=xmatch_table[optRACol+'_opt'].value*u.deg, dec=xmatch_table[optDecCol+'_opt'].value*u.deg)
	separation_rad_opt = coords_rad.separation(coords_opt).deg
	
	xmatch_table['rad_opt_sep_deg'] = separation_rad_opt
	
	if(len(xmatch_table) == 0):
		raise RuntimeError("\n%s: No cross-matched objects...!" %catalogName)
	
	if(makePlots):
		# Plotting optical and radio sources
		plotOutName = "%s/RadOptSkyPlot_%s.png" %(outPath, outSubscript)
		if not os.path.exists(plotOutName):
			plt.figure(figsize=(6, 6))	
			plt.scatter(optical_sources[optRACol], optical_sources[optDecCol], s=1, c='#8AD5F1', label='%s (N=%d)'%(opt_survey, len(optical_sources)))
			plt.scatter(radio_sources[radRACol], radio_sources[radDecCol], s=2, c='#87340D', label='MeerKAT (N=%d)' %len(radio_sources))
			plt.scatter(xmatch_table[radRACol+'_rad'], xmatch_table[radDecCol+'_rad'], marker='o', facecolor='None', linewidth=0.5, s=15, edgecolor='#06471D', label='MeerKATx%s (N=%d)' %(opt_survey, len(xmatch_table)))
			plt.title(catalogName)
			plt.title("%s\nSearch radius = %0.1f asec, %s band, Q0=%0.2f" %(catalogName, search_radius_asec.value, opt_mag_col, Q0), fontproperties=prop)
			plt.xlabel("RA (deg; J2000)", fontproperties=prop)
			plt.ylabel("Dec (deg; J2000)", fontproperties=prop)
			lgnd = plt.legend(loc="lower left", scatterpoints=1, fontsize=10, prop=prop)
			lgnd.legend_handles[0]._sizes = [30]
			lgnd.legend_handles[1]._sizes = [30]
			lgnd.legend_handles[2]._sizes = [30]
			plt.savefig(plotOutName , dpi=300, bbox_inches = 'tight')
			plt.close()
		
		# Plotting LR and Rel
		
		LRRelplotOutName = "%s/LRRelMagPlot_%s.png" %(outPath, outSubscript)
		if not os.path.exists(LRRelplotOutName):
			fig,ax=plt.subplots(nrows=2,ncols=1,sharex=True, sharey=False)	
			fig.set_size_inches(5,5)
			ax[0].scatter(xmatch_table[opt_mag_col+'_opt'], xmatch_table['LRs'], s=1, c='#87340D')
			ax[0].set_ylabel('LRs', fontproperties=prop)
			ax[0].set_yscale('log')
			ax[1].scatter(xmatch_table[opt_mag_col+'_opt'], xmatch_table['Rels'], s=1, c='#87340D')
			ax[1].set_ylabel('Rels', fontproperties=prop)
			ax[1].set_xlabel(opt_mag_col+'-band magnitude', fontproperties=prop)
			plt.suptitle(catalogName, fontproperties=prop)
			plt.subplots_adjust(hspace=0.0,wspace=0.0)
			plt.savefig(LRRelplotOutName , dpi=300, bbox_inches = 'tight')
			plt.close()
		
	return xmatch_table

def xmatch_berk(radio_cat, radio_band, xmatchTabOutName, opt_survey='DECaLS', opt_survey_dr='DR10', opt_mag_col = 'r', search_radius_asec = 4.0, makePlots=False, radRACol='RA', radDecCol='DEC', eRadRACol='E_RA', eRadDecCol='E_DEC', outSubscript=''):

	xmatchDirPath = startup.config['productsDir']+os.path.sep+'xmatches'
	
	catalogName = radio_cat.split(os.path.sep)[-1]
	captureBlockId = catalogName.split('_')[3]
	targetName = (catalogName.split('_1024ch_')[1]).split('_srl_')[0]
	
	# checking if this catalog is listed as having no optical sources in the area
	no_optical_sources_filename = startup.config['productsDir']+os.path.sep+'no_%s_sources.txt' %opt_survey
	if(os.path.exists(no_optical_sources_filename)):
		with open(no_optical_sources_filename, 'r') as file:
			no_optical_sources_catnames = [line.strip() for line in file]
		if catalogName in no_optical_sources_catnames:
			return
			
	# checking if this catalog is listed as having no optical counterparts in the optical
	no_xmatches_filename = startup.config['productsDir']+os.path.sep+'no_matches_%s.txt' %('_'.join(outSubscript.split('_')[1:]))
	if(os.path.exists(no_xmatches_filename)):
		with open(no_xmatches_filename, 'r') as file:
			no_xmatches_catnames = [line.strip() for line in file]
		if catalogName in no_xmatches_catnames:
			return

	if(opt_survey == 'DECaLS'):
	
		try:
			xmatch_tab = xmatchOpt(catalog=radio_cat, opt_survey=opt_survey, opt_survey_dr=opt_survey_dr, opt_mag_col=opt_mag_col, search_radius_asec=search_radius_asec, outPath=xmatchDirPath, makePlots=makePlots, radRACol=radRACol, radDecCol=radDecCol, eRadRACol=eRadRACol, eRadDecCol=eRadDecCol, outSubscript=outSubscript)
			
			if 'captureBlockId' not in xmatch_tab.columns:
				xmatch_tab.add_column(captureBlockId, name='captureBlockId', index=0) 
			if 'object' not in xmatch_tab.columns:
				xmatch_tab.add_column(targetName, name='object', index=1)
			if 'band' not in xmatch_tab.columns:
				xmatch_tab.add_column(radio_band, name='band')
				
		except RuntimeError as e:
			print(e)
			if "No optical sources found" in str(e):
				_write_no_optical_catFile(catalogName, no_optical_sources_filename)
			if "No cross-matched objects" in str(e):
				_write_no_optical_catFile(catalogName, no_xmatches_filename)
			return
			
		xmatch_tab.write(xmatchTabOutName, format='fits', overwrite=True)
		print("\nWrote cross-matched table %s." %xmatchTabOutName)
	
	else:
		print("\nERROR: Optical survey not mentioned.")
		return 
		
	
	
	return xmatch_tab
