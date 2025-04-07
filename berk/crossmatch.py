"""

Routines for cross matching, using likelihood ratio method.

"""

import os
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table, hstack
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
from zCluster import retrievers
# from matplotlib import font_manager as fm
# this_platform = os.environ['THIS_PLATFORM']
# if(this_platform == 'hp455'): #krishna's laptop
#     fpath = "/home/krishna/Dropbox/fonts/cmunss.ttf"
#     prop = fm.FontProperties(fname=fpath,size=12,math_fontfamily='stixsans')
# else:
# TODO: Fix this - this global variable is not what we want, but left in for now
from matplotlib import font_manager as fm
prop = fm.FontProperties(size=12,math_fontfamily='stixsans')

from . import startup, tasks


def _write_no_optical_cat_file(catalog_name, outfilename):
    """Insert docstring info

    """
    if os.path.exists(outfilename):
        with open(outfilename, 'a', encoding='utf8') as outfile:
            outfile.write("%s\n" %(catalog_name))
    else:
        with open(outfilename, 'w', encoding='utf8') as outfile:
            outfile.write("%s\n" %(catalog_name))
    return 0


def get_mag_hist(opt_mag, nbins):
    """Insert docstring info

    """
    opt_mag_hist, opt_mag_bins = np.histogram(opt_mag, bins=nbins)
    return opt_mag_hist, opt_mag_bins


def get_n_m(opt_mag, nbins, area_asec):
    """Insert docstring info

    """
    opt_mag_hist = get_mag_hist(opt_mag, nbins)[0]
    # TODO : check if it is whole area of survey or area of the overlapping regions
    n_m = opt_mag_hist/area_asec.value
    return n_m


def get_matches(radio_cat_table, optical_cat_table, search_radius_asec, rad_ra_col, rad_dec_col,
                opt_ra_col, opt_dec_col):
    """Insert docstring info

    """

    radio_coords = SkyCoord(ra=radio_cat_table[rad_ra_col].value * u.deg,
                            dec=radio_cat_table[rad_dec_col].value * u.deg)
    optical_coords = SkyCoord(ra=optical_cat_table[opt_ra_col].value * u.deg,
                              dec=optical_cat_table[opt_dec_col].value * u.deg)

    radio_cartesian = np.vstack(radio_coords.cartesian.xyz).T
    optical_cartesian = np.vstack(optical_coords.cartesian.xyz).T

    optical_tree = cKDTree(optical_cartesian)
    search_radius_rad = search_radius_asec.to(u.rad).value

    matches = optical_tree.query_ball_point(radio_cartesian, r=search_radius_rad)

    return matches


def get_separation(radio_cat_table, optical_cat_table, search_radius_asec, rad_ra_col, rad_dec_col,\
                   opt_ra_col, opt_dec_col):
    """Insert docstring info

    """
    radio_coords = SkyCoord(ra=radio_cat_table[rad_ra_col].value * u.deg,
                            dec=radio_cat_table[rad_dec_col].value * u.deg)
    optical_coords = SkyCoord(ra=optical_cat_table[opt_ra_col].value * u.deg,
                              dec=optical_cat_table[opt_dec_col].value * u.deg)
    matches = get_matches(radio_cat_table, optical_cat_table, search_radius_asec,
                          rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col)

    separation_list = []
    no_match_radio_idx = [] # List of radio sources with no match
    match_radio_idx = [] # List of radio sources with atleast one match
    match_optical_idx = [] # List of ALL optical sources that are matched to >= one radio source

    for radio_idx, optical_indices in enumerate(matches):
        if len(optical_indices) == 0: # radio sources with no optical matches in the search radius
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
    """Insert docstring info

    """
    return np.sqrt(sigma_rad**2 + sigma_opt**2)


def f_r(radio_coord, sigma_rad, optical_coord, sigma_opt):
    """Insert docstring info

    """
    r = radio_coord.separation(optical_coord).value
    sigma_pos = combined_error(sigma_rad, sigma_opt)

    return (1 / (2 * np.pi * sigma_pos**2)) * np.exp(-0.5 * (r**2 / sigma_pos**2))


def random_points_in_circle(center_ra, center_dec, radius_deg, n_points):
    """Insert docstring info

    """
    center_coord = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg, frame='icrs')

    rand_r = np.sqrt(np.random.uniform(0, radius_deg**2, n_points))  # Radius
    rand_theta = np.random.uniform(0, 2 * np.pi, n_points)  # Angle

    rand_ra = center_ra + rand_r * np.cos(rand_theta) / np.cos(center_dec * u.deg.to(u.rad))
    rand_dec = center_dec + rand_r * np.sin(rand_theta)

    coords = SkyCoord(ra=rand_ra * u.deg, dec=rand_dec * u.deg, frame='icrs')
    distances = center_coord.separation(coords)
    within_circle = distances < radius_deg * u.deg
    return rand_ra[within_circle], rand_dec[within_circle]


def make_random_cat(center_ra, center_dec, radius_deg, n_random_points, rad_ra_col, rad_dec_col):
    """Insert docstring info

    """
    rand_ra, rand_dec = random_points_in_circle(center_ra, center_dec, radius_deg, n_random_points)
    random_cat = Table([rand_ra, rand_dec], names=(rad_ra_col, rad_dec_col))
    return random_cat


def get_uobs_urand_ratio(radio_sources, rand_radio_sources, optical_sources, search_radius_asec,
                         rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col):
    """Insert docstring info

    """
    uobs_r = len(get_separation(radio_sources, optical_sources, search_radius_asec,
                                rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col)[3])
    urand_r = len(get_separation(rand_radio_sources, optical_sources, search_radius_asec,
                                 rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col)[3])
    return uobs_r/urand_r


def uobs_urand_ratio_model(r , Q0, sigma):
    """Insert docstring info

    """
    return 1. - (Q0*(1. - np.exp(-1*np.power(r,2.)/2.*np.power(sigma,2.))))


def get_Q0(radio_sources, rand_radio_sources, optical_sources, beam_size, rad_ra_col, rad_dec_col,
           opt_ra_col, opt_dec_col):
    """Insert docstring info

    """
    radii = np.arange( 1., beam_size, 0.5 )
    uobs_urand_ratios = []
    for radius in radii:
        search_radius_asec = float(radius)*u.arcsec
        uobs_urand_ratios.append(get_uobs_urand_ratio(radio_sources, rand_radio_sources,
                                                      optical_sources, search_radius_asec,
                                                      rad_ra_col, rad_dec_col, opt_ra_col,
                                                      opt_dec_col))
    params_fit, params_cov = curve_fit(uobs_urand_ratio_model, radii, uobs_urand_ratios)
    Q0, sig = params_fit
    Q0_err, sig_err = np.sqrt(np.diag(params_cov))
    return Q0, Q0_err


def get_q_m(radio_sources, optical_sources, n_m, search_radius_asec,
            opt_mag_col, Q0, n_mag_bins, rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col):
    """Insert docstring info

    """
    n_radio = len(radio_sources)
    ####
    # New
    result=get_separation(radio_sources, optical_sources, search_radius_asec, rad_ra_col,
                          rad_dec_col, opt_ra_col, opt_dec_col)

    match_optical_idx=result[2]
    ####
    # Old
    # separation_list, match_radio_idx, match_optical_idx, no_match_radio_idx = \
        # get_separation(radio_sources, optical_sources, search_radius_asec, rad_ra_col,
                         # rad_dec_col, opt_ra_col, opt_dec_col)
    ####
    opt_match_mags = optical_sources[opt_mag_col][match_optical_idx]
    mag_bins = get_mag_hist(optical_sources[opt_mag_col], n_mag_bins)[1]
    total_m = np.histogram(opt_match_mags, bins=mag_bins)[0]
    #real_m = total_m - (n_m * n_radio * np.pi * search_radius_asec.value**2)
    real_m = np.maximum(total_m - (n_m * n_radio * np.pi * search_radius_asec.value**2), 0) #TODO
    q_m = (real_m*Q0)/np.sum(real_m)
    return q_m


def get_qm_nm(mag, mag_bins, qm_nm_values):
    """Insert docstring info

    """
    bin_index = np.digitize(mag, mag_bins, right=True) - 1
    if 0 <= bin_index < len(qm_nm_values):
        return qm_nm_values[bin_index]
    #print("Issue with getting qm_nm")
    return np.nan


def make_mag_bins(magnitudes, nbins):
    """Insert docstring info

    """
    binwidth = (np.max(magnitudes)-np.min(magnitudes))/nbins
    mag_bins = np.arange( np.floor(np.min(magnitudes)), np.ceil(np.max(magnitudes)), binwidth )
    if np.ceil(np.max(magnitudes)) > np.max( mag_bins ):
        mag_bins = np.append( mag_bins, np.max(mag_bins)+0.4 )
    if np.floor(np.min(magnitudes)) < np.min( mag_bins ):
        mag_bins = np.insert( mag_bins, 0, np.min( mag_bins )-0.4 )
    # NOTE: below was (mag_bins) for some reason - not sure this was intentional?
    return mag_bins


def add_postfix_to_columns(table, postfix):
    """Insert docstring info

    """
    new_names = {colname: colname + postfix for colname in table.colnames}
    table.rename_columns(list(new_names.keys()), list(new_names.values()))
    return table
    
def get_centre_radius(radio_cat, rad_ra_col, rad_dec_col):
    """Insert docstring info

    """
    ra_list = list(radio_cat[rad_ra_col])
    dec_list = list(radio_cat[rad_dec_col])
    coords = SkyCoord(ra=ra_list * u.deg, dec=dec_list * u.deg, frame='icrs')

    # Centre from bounding box
    center_ra = (np.min(ra_list) + np.max(ra_list)) / 2
    center_dec = (np.min(dec_list) + np.max(dec_list)) / 2
    center_coord = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg, frame='icrs')

    # True maximum angular distance to any source
    separations = center_coord.separation(coords)
    radius_deg_max = np.max(separations).deg

    return center_ra, center_dec, radius_deg_max
    


def retrieve_decals(center_ra, center_dec, radius_deg, DR='DR10'):
    """Insert docstring info

    """
    print("\nRetrieving DECaLS %s sources with RAcen=%f, Deccen=%f, and radius=%f deg" \
          % (DR, center_ra, center_dec, radius_deg))
    if DR == 'DR10':
        decals_cat = Table(retrievers.DL_DECaLSDR10Retriever(center_ra, center_dec,
                                                             halfBoxSizeDeg = radius_deg,
                                                             DR = None))
    elif DR == 'DR8':
        # NOTE: There is no such routine in zCluster currently
        decals_cat = Table(retrievers.DL_DECaLSDR8Retriever(center_ra, center_dec,
                                                            halfBoxSizeDeg = radius_deg,
                                                            DR = None))
    else:
        raise Exception("Need to give DR10 or DR8")
    return decals_cat


def ompute_lr_rel(radio_sources, optical_sources, search_radius_asec, opt_mag_col, mag_bins, Q0,
                   qm_nm_values, rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col,
                   e_rad_ra_col, e_rad_dec_col, opt_pos_err_col):
    """Insert docstring info

    """
    matches = get_matches(radio_sources, optical_sources, search_radius_asec,
                          rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col)

    xmatch_optical_idx_list = []
    xmatch_radio_idx_list = []
    xmatch_lr_list = []
    xmatch_rel_list = []
    xmatch_phax_xm_list = []
    xmatch_pthix_xm_list = []

    for radio_idx, optical_indices in enumerate(matches):
        n_cand = len(optical_indices)

        if n_cand == 0: # TODO
            continue

        radio_coord = SkyCoord(ra=radio_sources[rad_ra_col][radio_idx] * u.deg,
                               dec=radio_sources[rad_dec_col][radio_idx] * u.deg)
        sigma_radio = np.sqrt(radio_sources[e_rad_ra_col][radio_idx]**2 + \
                              radio_sources[e_rad_dec_col][radio_idx]**2)

        lrs_radio_source = []  # List of LRs of all matches for a given radio source
        for optical_idx in optical_indices:
            optical_coord = SkyCoord(ra=optical_sources[opt_ra_col][optical_idx] * u.deg,
                                     dec=optical_sources[opt_dec_col][optical_idx] * u.deg)
            sigma_optical = optical_sources[opt_pos_err_col][optical_idx]

            F_R = f_r(radio_coord=radio_coord, sigma_rad=sigma_radio, optical_coord=optical_coord,
                      sigma_opt=sigma_optical)
            QM_NM = get_qm_nm(optical_sources[opt_mag_col][optical_idx], mag_bins, qm_nm_values)

            lr_given_match = QM_NM * F_R

            lrs_radio_source.append(lr_given_match)

        # Compute reliability
        sum_lrs_radio_source = np.sum(lrs_radio_source)
        denominator = sum_lrs_radio_source + (1 - Q0)

        # Handle edge cases
        if denominator <= 0:
            raise ValueError(f"Invalid denominator for radio source {radio_idx}: {denominator}")

        # List of Rels of all matches for a given radio source
        rels_radio_source = [LR / denominator for LR in lrs_radio_source]
        p_hasxm_radio_source = sum_lrs_radio_source / denominator

        # Collect data for cross-matched table
        for i, optical_idx in enumerate(optical_indices):
            xmatch_optical_idx_list.append(optical_idx)
            xmatch_radio_idx_list.append(radio_idx)
            xmatch_lr_list.append(lrs_radio_source[i])
            xmatch_rel_list.append(rels_radio_source[i])
            xmatch_phax_xm_list.append(p_hasxm_radio_source)
            xmatch_pthix_xm_list.append(rels_radio_source[i])

    radio_match_rows = radio_sources[xmatch_radio_idx_list]
    optical_match_rows = optical_sources[xmatch_optical_idx_list]

    radio_match_rows = add_postfix_to_columns(radio_match_rows, '_rad')
    optical_match_rows = add_postfix_to_columns(optical_match_rows, '_opt')

    # Combine data
    xmatch_table = hstack([radio_match_rows, optical_match_rows])
    xmatch_table['LRs'] = xmatch_lr_list
    xmatch_table['Rels'] = xmatch_rel_list
    xmatch_table['prob_has_match'] = xmatch_phax_xm_list
    xmatch_table['prob_this_match'] = xmatch_pthix_xm_list

    return xmatch_table


def xmatch_opt(catalog, opt_survey, opt_survey_dr, opt_mag_col, search_radius_asec, out_path,
              make_plots, rad_ra_col, rad_dec_col, e_rad_ra_col, e_rad_dec_col, out_subscript,
              save_files = True):
    """Insert docstring info

    """
    catalog_name = catalog.split(os.path.sep)[-1].replace(".fits","")

    search_radius_asec = search_radius_asec * u.arcsec
    n_mag_bins = 15
    beam_size = 7.0
    opt_pos_err = (0.2*u.arcsec).to(u.deg).value

    radio_sources = Table.read(catalog, format='fits', hdu=1)

    rad_pos_columns = [rad_ra_col, rad_dec_col, e_rad_ra_col, e_rad_dec_col]
    if not all(col in radio_sources.colnames for col in rad_pos_columns):
        print("RA and Dec columns of radio catalogue is not well set. Exiting!")
        return None

    if any(radio_sources[rad_ra_col] < 0.0):
        print("\nCrossmatching - Fixing RA for %s" %catalog_name)
        radio_sources = catalogs.fixRA(radio_sources, racol=rad_ra_col, wrap_angle=360)

    center_ra, center_dec, rad_ra, rad_dec, radius_deg = get_centre_radius(radio_sources,
                                                                           rad_ra_col, rad_dec_col)
    n_radio = len(radio_sources)
    sky_area_sqdeg = np.pi*rad_ra*rad_dec*u.deg**2
    sky_area_sqasec = sky_area_sqdeg.to(u.arcsec**2)

    # Collecting optical sources
    if opt_survey == 'DECaLS':
        opt_ra_col, opt_dec_col = 'RADeg', 'decDeg'
        opt_pos_err_col = 'pos_err'
        optical_sources = retrieve_decals(center_ra, center_dec, radius_deg, DR=opt_survey_dr)
        if len(optical_sources) == 0:
            raise RuntimeError("\n%s: No optical sources found in %s database...!" \
                                % (catalog_name, opt_survey))
        print("\n%d %s sources found in this region..." % (len(optical_sources), opt_survey))
    else:
        raise RuntimeError("\n%s: Optical survey not mentioned....!" %catalog_name)

    if opt_pos_err_col not in optical_sources.colnames:
        optical_sources[opt_pos_err_col] = opt_pos_err

    rand_radio_sources=make_random_cat(center_ra, center_dec, radius_deg, n_radio,
                                       rad_ra_col, rad_dec_col)

    print("\nFinding Q0...")
    #Q0 = 0.7024573416919023 #TODO
    Q0, Q0_err = get_Q0(radio_sources, rand_radio_sources, optical_sources, beam_size,
                        rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col)
    print("Q0=", Q0)
    n_m = get_n_m(opt_mag=list(optical_sources[opt_mag_col]), nbins=n_mag_bins,
                  area_asec=sky_area_sqasec)
    q_m = get_q_m(radio_sources, optical_sources, n_m, search_radius_asec,
                  opt_mag_col, Q0, n_mag_bins, rad_ra_col, rad_dec_col, opt_ra_col, opt_dec_col)
    qm_nm_values = q_m/n_m
    mag_bins = make_mag_bins(optical_sources[opt_mag_col], n_mag_bins)

    if save_files:
        optical_sources.write(out_path+os.path.sep+'decals_sources_%s.fits' %out_subscript,
                              format='fits', overwrite=True)
        rand_radio_sources.write(out_path+os.path.sep+'Randoms_%s.fits' %out_subscript,
                                 format='fits', overwrite=True)
        np.savetxt(out_path+os.path.sep+'Q0_%s.txt' %out_subscript, [Q0], fmt='%f')

    print("\nComputing LR and Reliabilities...")
    xmatch_table = ompute_lr_rel(radio_sources, optical_sources, search_radius_asec, opt_mag_col,
                                  mag_bins, Q0, qm_nm_values, rad_ra_col, rad_dec_col,
                                  opt_ra_col, opt_dec_col, e_rad_ra_col, e_rad_dec_col,
                                  opt_pos_err_col)

    # computing cross-match separation
    coords_rad = SkyCoord(ra=xmatch_table[rad_ra_col+'_rad'].value*u.deg,
                          dec=xmatch_table[rad_dec_col+'_rad'].value*u.deg)
    coords_opt = SkyCoord(ra=xmatch_table[opt_ra_col+'_opt'].value*u.deg,
                          dec=xmatch_table[opt_dec_col+'_opt'].value*u.deg)
    separation_rad_opt = coords_rad.separation(coords_opt).deg

    xmatch_table['rad_opt_sep_deg'] = separation_rad_opt

    if len(xmatch_table) == 0:
        raise RuntimeError("\n%s: No cross-matched objects...!" %catalog_name)

    if make_plots:
        # Plotting optical and radio sources
        plot_out_name = "%s/RadOptSkyPlot_%s.png" %(out_path, out_subscript)
        if not os.path.exists(plot_out_name):
            plt.figure(figsize=(6, 6))
            plt.scatter(optical_sources[opt_ra_col], optical_sources[opt_dec_col], s=1,
                        c='#8AD5F1', label='%s (N=%d)'%(opt_survey, len(optical_sources)))
            plt.scatter(radio_sources[rad_ra_col], radio_sources[rad_dec_col], s=2, c='#87340D',
                        label='MeerKAT (N=%d)' %len(radio_sources))
            plt.scatter(xmatch_table[rad_ra_col+'_rad'], xmatch_table[rad_dec_col+'_rad'],
                        marker='o', facecolor='None', linewidth=0.5, s=15,
                        edgecolor='#06471D',
                        label='MeerKATx%s (N=%d)' %(opt_survey, len(xmatch_table)))
            plt.title(catalog_name)
            plt.title("%s\nSearch radius = %0.1f asec, %s band, Q0=%0.2f" \
                      % (catalog_name, search_radius_asec.value, opt_mag_col, Q0),
                      fontproperties=prop)
            plt.xlabel("RA (deg; J2000)", fontproperties=prop)
            plt.ylabel("Dec (deg; J2000)", fontproperties=prop)
            lgnd = plt.legend(loc="lower left", scatterpoints=1, fontsize=10, prop=prop)
            lgnd.legend_handles[0]._sizes = [30]
            lgnd.legend_handles[1]._sizes = [30]
            lgnd.legend_handles[2]._sizes = [30]
            plt.savefig(plot_out_name , dpi=300, bbox_inches = 'tight')
            plt.close()

        # Plotting LR and Rel
        lrrelplot_out_name = "%s/LRRelMagPlot_%s.png" %(out_path, out_subscript)
        if not os.path.exists(lrrelplot_out_name):
            fig,ax=plt.subplots(nrows=2,ncols=1,sharex=True, sharey=False)
            fig.set_size_inches(5,5)
            ax[0].scatter(xmatch_table[opt_mag_col+'_opt'], xmatch_table['LRs'], s=1, c='#87340D')
            ax[0].set_ylabel('LRs', fontproperties=prop)
            ax[0].set_yscale('log')
            ax[1].scatter(xmatch_table[opt_mag_col+'_opt'], xmatch_table['Rels'], s=1, c='#87340D')
            ax[1].set_ylabel('Rels', fontproperties=prop)
            ax[1].set_xlabel(opt_mag_col+'-band magnitude', fontproperties=prop)
            plt.suptitle(catalog_name, fontproperties=prop)
            plt.subplots_adjust(hspace=0.0,wspace=0.0)
            plt.savefig(LRRelplot_out_name , dpi=300, bbox_inches = 'tight')
            plt.close()

    return xmatch_table


def xmatch_berk(radio_cat, radio_band, xmatch_tab_out_name, opt_survey='DECaLS',
                opt_survey_dr='DR10', opt_mag_col = 'r', search_radius_asec = 4.0,
                make_plots=False, rad_ra_col='RA', rad_dec_col='DEC',
                e_rad_ra_col='E_RA', e_rad_dec_col='E_DEC', out_subscript=''):
    """Insert docstring info

    """

    xmatch_dir_path = startup.config['productsDir']+os.path.sep+'xmatches'

    catalog_name = radio_cat.split(os.path.sep)[-1]
    capture_block_id = catalog_name.split('_')[3]
    target_name = (catalog_name.split('_1024ch_')[1]).split('_srl_')[0]

    # checking if this catalog is listed as having no optical sources in the area
    no_optical_sources_filename = startup.config['productsDir']+os.path.sep+'no_%s_sources.txt'\
                                  %opt_survey
    if os.path.exists(no_optical_sources_filename):
        with open(no_optical_sources_filename, 'r', encoding="utf-8") as infile:
            no_optical_sources_catnames = [line.strip() for line in infile]
        if catalog_name in no_optical_sources_catnames:
            return None

    # checking if this catalog is listed as having no optical counterparts in the optical
    no_xmatches_filename = startup.config['productsDir']+os.path.sep+'no_matches_%s.txt' \
                           %('_'.join(out_subscript.split('_')[1:]))
    if os.path.exists(no_xmatches_filename):
        with open(no_xmatches_filename, 'r', encoding="utf-8") as infile:
            no_xmatches_catnames = [line.strip() for line in infile]
        if catalog_name in no_xmatches_catnames:
            return None

    if opt_survey == 'DECaLS':
        try:
            xmatch_tab = xmatch_opt(catalog=radio_cat, opt_survey=opt_survey,
                                    opt_survey_dr=opt_survey_dr, opt_mag_col=opt_mag_col,
                                    search_radius_asec=search_radius_asec,
                                    out_path=xmatch_dir_path, make_plots=make_plots,
                                    rad_ra_col=rad_ra_col, rad_dec_col=rad_dec_col,
                                    e_rad_ra_col=e_rad_ra_col, e_rad_dec_col=e_rad_dec_col,
                                    out_subscript=out_subscript, save_files=False)

            if 'captureBlockId' not in xmatch_tab.columns:
                xmatch_tab.add_column(capture_block_id, name='captureBlockId', index=0)
            if 'object' not in xmatch_tab.columns:
                xmatch_tab.add_column(target_name, name='object', index=1)
            if 'band' not in xmatch_tab.columns:
                xmatch_tab.add_column(radio_band, name='band')

        except RuntimeError as e:
            print(e)
            if "No optical sources found" in str(e):
                _write_no_optical_cat_file(catalog_name, no_optical_sources_filename)
            if "No cross-matched objects" in str(e):
                _write_no_optical_cat_file(catalog_name, no_xmatches_filename)
            return None

        xmatch_tab.write(xmatch_tab_out_name, format='fits', overwrite=True)
        print("\nWrote cross-matched table %s." %xmatch_tab_out_name)
    else:
        print("\nERROR: Optical survey not mentioned.")
        return None

    return xmatch_tab
