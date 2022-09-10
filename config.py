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
        local_root='/data/user/shiqiyu/ESTES/data/', username=username)
#/home/shiqiyu/analysis/wg-nu-sources/ESTES_seyfert/data/', username=username)
    template_repo = cy.selections.mrichman_repo
    #base_dir = cy.utils.ensure_dir('/home/shiqiyu/analysis/wg-nu-sources/ESTES_seyfert/data/analyses/{}'.format(job_base))
    base_dir = cy.utils.ensure_dir('/data/user/shiqiyu/ESTES/data/analyses/{}'.format(job_base))

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
    if gamma!=0:               
        print("config use gamma = ",  gamma)
        gamma=3
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
    pl_fluxs = []
    src_dist=np.atleast_1d(src_dist)
    src_log_lumin=np.atleast_1d(src_log_lumin)
    for dist, log_lumin in zip(src_dist, src_log_lumin):
        sflux = SeyfertFlux(seyfert_flux_table_path, log_lumin)
        psp_table = sflux.get_splinetable()
        rel_flux_scale = sflux._relative_flux_scale
        flux_model=cy.hyp.SeyfertCoreCoronaFlux(psp_table, log_lumin, dist, norm=1.0,  lumin_scale=rel_flux_scale)
        pl_flux = cy.hyp.PowerLawFlux(gamma, energy_cutoff=np.inf)
        pl_fluxs.append(pl_flux)
        fluxs.append(flux_model)
    if weighted_src is None:
        weighted_src = src
    conf = {
        'src': weighted_src,
        'corona_flux': True,
        'energy': 'customflux', 
        'flux': fluxs, ##flux for fit
        'update_bg': True,
        'sigsub':  sigsub,
        'randomize': ['ra', cy.inj.DecRandomizer],
        'sindec_bandwidth': np.radians(5),
        'dec_rand_method': 'gaussian_fixed',
        'dec_rand_kwargs': dict(randomization_width=np.radians(3)),
        'dec_rand_pole_exlusion': np.radians(8),
        'inj_conf':{
            'src':src,
            'corona_flux': True,
            'update_bg':True,
            'flux': fluxs,
            'randomize': ['ra', cy.inj.DecRandomizer],
            'sindec_bandwidth': np.radians(5),
            'dec_rand_method': 'gaussian_fixed',
            'dec_rand_kwargs': dict(randomization_width=np.radians(3)),
            'dec_rand_pole_exlusion': np.radians(8)
            }
    }

    if False: #not sigsub:
        mc_conf = {
          'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
          'randomize'      :['ra', 'dec'],
          'update_bg': True,
          'sigsub': False,
          'inj_conf':{
            'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
            'randomize'      :['ra', 'dec'],
            'update_bg':False,
            }
        }
        conf.update(mc_conf)

    if gamma != 0:
        print(" Fitting in power law. pop out energy from config so that energy can be set as fit later by csky by default")
        conf.pop('energy')
        conf['flux'] = pl_fluxs

    return conf

def get_gp_bg_ps_conf(template_str, sigmas, 
        src, src_dist, src_log_lumin, gp_norm = 21.8e-12, #(DNN's paper) kra=0.89
        min_dec=-90, max_dec=90, gamma=None, cutoff_GeV=np.inf,
        base_dir=base_dir, repo=template_repo,
        weighted_src=None, seyfert_flux_table_path='/home/shiqiyu/northern_sky_seyferts/generate_splinetables/splinetables/', norm=1., sigsub=False, mcbg=False):
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


    if 'kra' in template_str or 'pi0' in template_str:

        # check that gamma isn't set
        #if gamma is not None:
            #raise ValueError(
                #'Gamma must not be specified for KRA, but is:', gamma)

        if template_str == 'kra5':
            template, energy_bins = repo.get_template(
                'KRA-gamma_5PeV_maps_energies', per_pixel_flux=True)
            krag5_map = set_min_max_dec(template, min_dec, max_dec)
            tmp_flux = cy.hyp.BinnedFlux(
                bins_energy=energy_bins,
                flux=template.sum(axis=0))
            template_dir = cy.utils.ensure_dir(
                '{}/templates/kra5'.format(base_dir))
            gamma_tmp = 2.7
        elif template_str == 'pi0':
            gamma_tmp = 2.7
            template = repo.get_template('Fermi-LAT_pi0_map')
            template = set_min_max_dec(template, min_dec, max_dec)
            template_dir = cy.utils.ensure_dir(
                '{}/templates/pi0/gamma/{:.3f}'.format(base_dir, gamma_tmp))
            tmp_flux = cy.hyp.PowerLawFlux(gamma_tmp, energy_cutoff=cutoff_GeV)
        #temp_model = cy.pdf.CustomFluxEnergyPDFRatioModel
        gp_conf = {
            'mcbg': mcbg,
            'template': template,
            'fitter_args': dict(gamma=gamma_tmp),
            'randomize': ['ra'],
            'update_bg': True,
            'sigmas': smears,
            'sigsub': False,
            'dir': template_dir,
            'bg_temp':True,
            'space': 'template',
            'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
            'flux':tmp_flux,
            'inj_conf':{
                'template': template,
#                'nsig': 68,
                'gp_norm': gp_norm,
                'space': 'template',
                'fitter_args': dict(gamma=gamma_tmp),
                'flux':tmp_flux,
                'sig': 'template',
#                'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
                'randomize'      :['ra'],
                'update_bg': True,
                },
#            'energy': False,
        }
    else:
        raise ValueError('Unknown template name: {}'.format(template_str))
    if template_str == 'kra5':
        update = {
            cy.pdf.CustomFluxEnergyPDFRatioModel: dict(
                hkw=dict(bins=(
                       np.linspace(-1, 1, 20),
                       np.linspace(np.log10(500), 8.001, 20)
                       )),
                flux=tmp_flux,
                features=['sindec', 'log10energy'],
                normalize_axes=([1])),
            'bins_energy': energy_bins,
        }
        gp_conf.pop('flux')
        gp_conf.update(update)

    from lib_bass_seyferts import SeyfertFlux
    fluxs=[]
    pl_fluxs=[]
    src_dist=np.atleast_1d(src_dist)
    src_log_lumin=np.atleast_1d(src_log_lumin)
    for dist, log_lumin in zip(src_dist, src_log_lumin):
        sflux = SeyfertFlux(seyfert_flux_table_path, log_lumin)
        psp_table = sflux.get_splinetable()
        rel_flux_scale = sflux._relative_flux_scale
        flux_model=cy.hyp.SeyfertCoreCoronaFlux(psp_table, log_lumin, dist, norm=1.0,  lumin_scale=rel_flux_scale)
        fluxs.append(flux_model)

        pl_flux = cy.hyp.PowerLawFlux(gamma, energy_cutoff=np.inf)
        pl_fluxs.append(pl_flux)

    if weighted_src is None:
        weighted_src = src

    ps_conf = {
        'src': weighted_src,
        'corona_flux': True,
        'energy': 'customflux',
        'flux': fluxs, ##flux for fit
        'update_bg': True,
        'sigsub':  sigsub,
        'inj_conf':{
            'src':src,
            'corona_flux': True,
            'bg_weight_names':['astro_weight', 'atmo_weight', 'muon_weight'],
            'randomize'      :['ra', 'dec'],
            'update_bg':True,
            'flux': fluxs,
            }
    }
    if gamma is not None:
        print(" Fitting in power law. pop out energy from config so that energy can be set as fit later by csky by default")
        ps_conf.pop('energy')
        ps_conf['inj_conf']['corona_flux']=False
        ps_conf['inj_conf']['flux']= pl_fluxs
    conf=ps_conf
    conf['gpbg_conf']=gp_conf

    return conf
