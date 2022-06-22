# config.py
import socket
import numpy as np
import csky as cy
import getpass
import utils
import os
import healpy as hp

hostname = socket.gethostname()
username = getpass.getuser()
print('Running as User: {} on Hostname: {}'.format(username, hostname))
job_base = 'baseline_analysis'
#job_base = 'systematics_full'
if 'condor00' in hostname or 'cobol' in hostname or 'gpu' in hostname:
    print("WRONG SUBMIT!!!")
    quit()
    repo = cy.selections.Repository(
        local_root='/data/i3store/users/ssclafani/data/analyses'.format(username))
    template_repo = repo
    ana_dir = cy.utils.ensure_dir(
        '/data/i3store/users/{}/data/analyses'.format(username))
    base_dir = cy.utils.ensure_dir(
        '/data/i3store/users/{}/data/analyses/{}'.format(username, job_base))
    job_basedir = '/data/i3home/{}/submitter_logs'.format(username)
else:
    repo = cy.selections.Repository(
        local_root='/home/shiqiyu/analysis/wg-nu-sources/ESTES_seyfert/data/', username=username)
    template_repo = cy.selections.mrichman_repo
    base_dir = cy.utils.ensure_dir('/home/shiqiyu/analysis/wg-nu-sources/ESTES_seyfert/data/analyses/{}'.format(job_base))
    ana_dir = '{}/ana'.format (base_dir)
    job_basedir = '/scratch/{}/'.format(username) 

# Path to submit config file. This needs to be a relative path to $HOME
# Example content of this file:
#    eval `/cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/setup.sh`
#    source  ~/path/to/venv/bin/activate
submit_cfg_file = 'analysis/wg-nu-sources/ESTES_seyfert/submitter_config'

# ---------------------------------------------
# Define csky config settings for trial runners
# ---------------------------------------------

def set_min_max_dec(gp_map, min_dec, max_dec):
    npix  = len(gp_map)
    nside = hp.npix2nside(npix)
    max_theta = np.pi/2. - np.radians(min_dec)
    min_theta = np.pi/2. - np.radians(max_dec)

    (theta, phi)  = hp.pix2ang(nside, np.arange(npix))
    mask = np.logical_or((theta < min_theta),(theta > max_theta))

    gp_map[mask] = 0
    return gp_map

def get_ps_conf(src, gamma, cutoff_GeV=np.inf, sigsub=True):
    """Get csky trial runner config for Point Source Likelihood

    Parameters
    ----------
    src : csky.utils.sources
        The sources.
    gamma : float
        The spectral index gamma to use for the powerlaw flux.
    cutoff_GeV : float, optional
        The cutoff value for the powerlaw flux.

    Returns
    -------
    dict
        The config, which may be passed to csky.get_trial_runner
    """
    if sigsub is False:
        print(utils.bcolors.YELLOW)
        print('=========================================================')
        print('=== Warning: trial runner is using no sigsub!         ===')
        print('=========================================================')
        print(utils.bcolors.ENDC)

    conf = {
        'src': src,
        'flux': cy.hyp.PowerLawFlux(gamma, energy_cutoff=cutoff_GeV),
        'update_bg': True,
        'sigsub':  sigsub,
        'randomize': ['ra', cy.inj.DecRandomizer],
        'sindec_bandwidth': np.radians(5),
        'dec_rand_method': 'gaussian_fixed',
        'dec_rand_kwargs': dict(randomization_width=np.radians(3)),
        'dec_rand_pole_exlusion': np.radians(8)
    }

    return conf

def get_seyfert_ps_conf(src, src_dist, src_log_lumin, gamma=0, weighted_src=None, seyfert_flux_table_path='/home/shiqiyu/northern_sky_seyferts/generate_splinetables/splinetables/', norm=1., sigsub=False):
#/data/user/shiqiyu/northern_sky_seyferts/generate_splinetables/splinetables/', norm=1., sigsub=True):
    """Get csky trial runner config for Point Source Likelihood

    Parameters
    ----------
    src : csky.utils.sources
        The sources.

    Returns
    -------
    dict
        The config, which may be passed to csky.get_trial_runner
    """
    if sigsub is False:
        print(utils.bcolors.YELLOW)
        print('=========================================================')
        print('=== Warning: trial runner is using no sigsub!         ===')
        print('=========================================================')
        print(utils.bcolors.ENDC)

    from lib_bass_seyferts import SeyfertFlux
    fluxs=[]
    src_dist=np.atleast_1d(src_dist)
    src_log_lumin=np.atleast_1d(src_log_lumin)
    for dist, log_lumin in zip(src_dist, src_log_lumin):
        sflux = SeyfertFlux(seyfert_flux_table_path, log_lumin)
        psp_table = sflux.get_splinetable()
        rel_flux_scale = sflux._relative_flux_scale

        flux_model=cy.hyp.SeyfertCoreCoronaFlux(psp_table, log_lumin, dist, norm=1.0,  lumin_scale=rel_flux_scale)
        fluxs.append(flux_model)
    if weighted_src is None:
        weighted_src = src
    conf = {
        'src': weighted_src,
        'corona_flux': True,
        'energy': 'customflux', 
        'flux': fluxs, ##flux for fit
        'update_bg': False,
        'sigsub':  sigsub,
        'inj_conf':{
            'src':src,
            'corona_flux': True,
            'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
            'randomize'      :['ra', 'dec'],
            'update_bg':False,
            'flux': fluxs,
            }
    }

    if gamma != 0:
        print(" Fitting in power law. pop out energy from config so that energy can be set as fit later by csky by default")
        conf.pop('energy')

    return conf

def get_gp_bg_ps_conf(template_str, sigmas, 
        src, src_dist, src_log_lumin, gp_norm = 0.89, #(DNN's paper)
        min_dec=-90, max_dec=90, gamma=None, cutoff_GeV=np.inf,
        base_dir=base_dir, repo=template_repo,
        weighted_src=None, seyfert_flux_table_path='/home/shiqiyu/northern_sky_seyferts/generate_splinetables/splinetables/', norm=1., sigsub=False):
    """Get csky trial runner config for Galactic Plane Template
    Parameters
    ----------
    template_str : str
        The name of the template to use. Must be one of:
        ['pi0', 'fermibubbles', 'kra5', 'kra50']
    gamma : float, optional
        The spectral index to use. This may only be set for Pi0 or
        Fermi Bubbles. Defaults to 2.7 for Pi0 and 2.0 for Fermi Bubbles.
    cutoff_GeV : float, optional
        The cutoff value for the powerlaw spectrum used for Fermi Bubble flux.
        This is only relevant for the Fermi Bubble template.
    base_dir : str, optional
        The path to the base directory. Will be used to cache the templates.
    repo : csky.selections.Repository, optional
        Csky data repository to use.
    Returns
    -------
    dict
        The config, which may be passed to csky.get_trial_runner
    Raises
    ------
    ValueError
        Description
    """

    # Print warning: GP templates have fixed gamma in our analysis.
    # Setting it to custom values should only be done for debugging/testing
    # purposes and the user should be aware of this.
    if gamma is not None:
        print(utils.bcolors.YELLOW)
        print('=========================================================')
        print('=== Warning: trial runner is using non-default gamma! ===')
        print('=========================================================')
        print(utils.bcolors.ENDC)
    smears = np.logspace(np.log10(min(sigmas)), np.log10(max(sigmas)), 50+1)


    if 'kra' in template_str:

        # check that gamma isn't set
        #if gamma is not None:
            #raise ValueError(
                #'Gamma must not be specified for KRA, but is:', gamma)

        if template_str == 'kra5':
            template, energy_bins = repo.get_template(
                'KRA-gamma_5PeV_maps_energies', per_pixel_flux=True)
            krag5_map = set_min_max_dec(template, min_dec, max_dec)

            kra_flux = cy.hyp.BinnedFlux(
                bins_energy=energy_bins,
                flux=template.sum(axis=0))
            template_dir = cy.utils.ensure_dir(
                '{}/templates/kra5'.format(base_dir))
        elif template_str == 'kra50':
            template, energy_bins = repo.get_template(
                      'KRA-gamma_maps_energies', per_pixel_flux=True)
            kra_flux = cy.hyp.BinnedFlux(
                bins_energy=energy_bins,
                flux=template.sum(axis=0))
            template_dir = cy.utils.ensure_dir(
                '{}/templates/kra50'.format(base_dir))
        temp_model = cy.pdf.CustomFluxEnergyPDFRatioModel
        gp_conf = {
            'template': template,
            'bins_energy': energy_bins,
            'fitter_args': dict(gamma=2.5),
            'randomize': ['ra'],
            'update_bg': False,
            'sigmas': smears,
            'sigsub': False,
            'dir': template_dir,
            'bg_temp':True,
            #'flux':kra_flux,
            'space': 'template',
            'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
            'inj_conf':{
#                'nsig': 68,
                'gp_norm': gp_norm,
#.0053756309758,
                'space': 'template',
                'flux':kra_flux,
                'sig': 'template',
                'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
                'randomize'      :['ra'],
                'update_bg': False,
                },
            cy.pdf.CustomFluxEnergyPDFRatioModel: dict(
                hkw=dict(bins=(
                       np.linspace(-1, 1, 20),
                       np.linspace(np.log10(500), 8.001, 20)
                       )),
                flux=kra_flux,
                features=['sindec', 'log10energy'],
                normalize_axes=([1])),
            'energy': False,
        }
    else:
        raise ValueError('Unknown template name: {}'.format(template_str))

    from lib_bass_seyferts import SeyfertFlux
    fluxs=[]
    src_dist=np.atleast_1d(src_dist)
    src_log_lumin=np.atleast_1d(src_log_lumin)
    for dist, log_lumin in zip(src_dist, src_log_lumin):
        sflux = SeyfertFlux(seyfert_flux_table_path, log_lumin)
        psp_table = sflux.get_splinetable()
        rel_flux_scale = sflux._relative_flux_scale

        flux_model=cy.hyp.SeyfertCoreCoronaFlux(psp_table, log_lumin, dist, norm=1.0,  lumin_scale=rel_flux_scale)
        fluxs.append(flux_model)
    if weighted_src is None:
        weighted_src = src
    ps_conf = {
        'src': weighted_src,
        'corona_flux': True,
        'energy': 'customflux',
        'flux': fluxs, ##flux for fit
        'update_bg': False,
        'sigsub':  sigsub,
        'inj_conf':{
            'src':src,
            'corona_flux': True,
            'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
            'randomize'      :['ra', 'dec'],
            'update_bg':False,
            'flux': fluxs,
            }
    }
    if gamma is not None:
        #raise ValueError('gamma has to be none, do not support pl fit atm')
        print(" Fitting in power law. pop out energy from config so that energy can be set as fit later by csky by default")
        ps_conf.pop('energy')
        ps_conf['inj_conf']['corona_flux']=False
        ps_conf['inj_conf']['flux']= cy.hyp.PowerLawFlux(gamma, energy_cutoff=cutoff_GeV)
    conf=ps_conf
    conf['gpbg_conf']=gp_conf

    return conf
