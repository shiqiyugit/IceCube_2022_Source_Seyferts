# config.py
import socket
import numpy as np
import csky as cy
import getpass
import utils
import os

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
        '/home/shiqiyu/analysis/wg-nu-sources/ESTES_seyfert/data/', '/data/user/smancina/ESTES_DataSet',  username=username)
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

def get_seyfert_ps_conf(src, src_dist, src_log_lumin, gamma=0, weighted_src=None, seyfert_flux_table_path='/data/user/shiqiyu/northern_sky_seyferts/generate_splinetables/splinetables/', norm=1., sigsub=False):
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

