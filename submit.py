#!/usr/bin/env python
try:
    import matplotlib
except ModuleNotFoundError:
    os.environ['MPLBACKEND'] = 'AGG'
    import matplotlib

import csky as cy
import numpy as np
import datetime, socket
from submitter import Submitter
now = datetime.datetime.now
import config as cg
import click, sys, os, time
flush = sys.stdout.flush

repo, ana_dir, base_dir, job_basedir = cg.repo, cg.ana_dir, cg.base_dir, cg.job_basedir
hostname = cg.hostname
username = cg.username
submit_cfg_file = cg.submit_cfg_file


class State (object):
    def __init__ (self, ana_name, ana_dir, save,  base_dir,  job_basedir, mask_deg=0, source_r=0, nsrc_tomask=None, mask_self=False):
        self.ana_name, self.ana_dir, self.save, self.job_basedir = ana_name, ana_dir, save, job_basedir
        self.base_dir = base_dir
        self.mask_deg=mask_deg
        self.nsrc_tomask = nsrc_tomask
        self.source_r = source_r
        self.mask_self=mask_self

        if self.mask_deg !=0:
            self.ana_name += '_masking{:08.3f}'.format(self.mask_deg)
        if self.source_r !=0:
            self.ana_name += '_top{}_{}deg'.format(nsrc_tomask, source_r)
        if self.mask_self:
            self.ana_name += '_selfmask'
        self._ana = None    
    @property
    def state_args (self):
        return '--ana {} --ana-dir {} --base-dir {} --mask_deg {}'.format (
            self.ana_name, self.ana_dir, self.base_dir, self.mask_deg)

pass_state = click.make_pass_decorator (State)

@click.group (invoke_without_command=True, chain=True)
@click.option ('-a', '--ana', 'ana_name', default='ESTES', help='Dataset title')
@click.option ('--ana-dir', default=ana_dir, type=click.Path ())
@click.option ('--job_basedir', default=job_basedir, type=click.Path ())
@click.option ('--save/--nosave', default=False)
@click.option ('--base-dir', default=base_dir,
               type=click.Path (file_okay=False, writable=True))
@click.option ('--mask_deg', default=0, type=float, help='+- degree along GP')
@click.option ('--source_r', default=0, type=float, help='radius degree around sources')
@click.option ('--nsrc_tomask', default=None, type=int)
@click.option ('--mask_self', default=False, type=bool)
@click.pass_context
def cli (ctx, ana_name, ana_dir, save, base_dir, job_basedir, mask_deg, source_r, nsrc_tomask, mask_self):
    ctx.obj = State.state = State (ana_name, ana_dir, save, base_dir, job_basedir, mask_deg, source_r, nsrc_tomask, mask_self)


@cli.resultcallback ()
def report_timing (result, **kw):
    exe_t1 = now ()
    print ('c7: end at {} .'.format (exe_t1))
    print ('c7: {} elapsed.'.format (exe_t1 - exe_t0))

@cli.command ()
@pass_state
def setup_ana (state):
    state.ana

@cli.command ()
@click.option ('--n-trials', default=10000, type=int)
@click.option ('--n-jobs', default=10, type=int)
@click.option ('-n', '--n-sig', 'n_sigs', multiple=True, default=[0], type=float)
@click.option ('--gamma', default=3.0, type=float)
@click.option ('-c', '--cutoff', default=np.inf, type=float, help='exponential cutoff energy (TeV)')      
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--sigsub/--nosigsub', default=True)
@click.option ('--dec', 'dec_degs', multiple=True, type=float, default=())
@click.option ('--dry/--nodry', default=False)
@click.option ('--seed', default=0)
@pass_state
def submit_do_ps_trials (
        state, n_trials, n_jobs, n_sigs, gamma, 
        cutoff,  poisson, sigsub, dec_degs, dry, 
        seed):
    ana_name = state.ana_name
    T = time.time ()
    poisson_str = 'poisson' if poisson else 'nopoisson'
    sigsub_str = 'sigsub' if sigsub else 'nosigsub'
    job_basedir = state.job_basedir 
    job_dir = '{}/{}/ps_trials/T_E{}_{:17.6f}'.format (
        job_basedir, ana_name, int(gamma * 100),  T)
    sub = Submitter (job_dir=job_dir, memory=5, 
        max_jobs=1000, config = submit_cfg_file)
    commands, labels = [], []
    #reqs = '(Machine != "cobol97.private.pa.umd.edu") & (Machine != "cobol94.private.pa.umd.edu")'
    trial_script = os.path.abspath('trials.py')
    dec_degs = dec_degs or np.r_[-81:+81:2]
    for dec_deg in dec_degs:
        for n_sig in n_sigs:
            for i in range (n_jobs):
                s = i + seed
                fmt = ' {} do-ps-trials --dec_deg={:+08.3f} --n-trials={}' \
                        ' --n-sig={} --gamma={:.3f} --cutoff={}' \
                        ' --{} --seed={} --{}'
  
                command = fmt.format (trial_script,  dec_deg, n_trials,
                                      n_sig, gamma, cutoff, poisson_str, s,sigsub_str)
                fmt = 'csky__dec_{:+08.3f}__trials_{:07d}__n_sig_{:08.3f}__' \
                        'gamma_{:.3f}_cutoff_{}_{}__seed_{:04d}'
                label = fmt.format (dec_deg, n_trials, n_sig, gamma,
                                    cutoff, poisson_str, s )
                commands.append (command)
                labels.append (label)
    sub.dry = dry
    print(hostname)
    if 'condor00' in hostname:
        sub.submit_condor00 (commands, labels) #, reqs=reqs)
    else:
        sub.submit_npx4 (commands, labels)

@cli.command ()
@click.option ('--n-trials', default=1000, type=int)
@click.option ('--n-jobs', default=1, type=int)
@click.option ('-n', '--n-sig', 'n_sigs', multiple=True, default=[0], type=float)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--sigsub/--nosigsub', default=True)
@click.option ('-nsrc', '--n-sources', 'n_srcs', default=14, type=int, help="the top n src to run")
@click.option ('--dry/--nodry', default=False)
@click.option ('--seed', default=0)
@click.option ('--gamma', default=0, type=float, help='0 is fit to model, otherwise call pl')
@pass_state
def submit_do_seyfert_ps_trials (
        state, n_trials, n_jobs, n_sigs,
        poisson, sigsub, n_srcs,
        dry, seed, gamma):
    ana_name = state.ana_name
    ana_command=''
    if state.mask_deg !=0 :
        ana_command +=' --mask_deg={}'.format(state.mask_deg)
    if state.source_r !=0:
        if state.nsrc_tomask>0:
            ana_command +=' --nsrc_tomask={} --source_r={}'.format(state.nsrc_tomask, state.source_r)

    T = time.time ()
    poisson_str = 'poisson' if poisson else 'nopoisson'
    sigsub_str = 'sigsub' if sigsub else 'nosigsub'
    job_basedir = state.job_basedir
    job_dir = '{}/{}/ps_trials/T_{:17.6f}'.format (
        job_basedir, ana_name, T)
        
    sub = Submitter (job_dir=job_dir, memory=2,
        max_jobs=1000, config = submit_cfg_file)
    commands, labels = [], []
    trial_script = os.path.abspath('trials.py')
    for ind in range(n_srcs):
        for n_sig in n_sigs:
            for i in range (n_jobs):
                s = i + seed
                fmt = ' {} {} do-seyfert-ps-trials --n-src={} --n-trials={}' \
                        ' --n-sig={} ' \
                        ' --{} --seed={} --{} --gamma={}'
                if state.mask_self:
                    mask_command =ana_command + ' --nth_tomask={}'.format(ind)
                    command = fmt.format (trial_script, mask_command, ind, n_trials,
                                      n_sig, poisson_str, s,sigsub_str, gamma)
                else:
                    command = fmt.format (trial_script, ana_command, ind, n_trials,
                                      n_sig, poisson_str, s,sigsub_str, gamma)
                fmt = 'csky__src_{:02d}__trials_{:07d}__n_sig_{:08.3f}__' \
                        '{}__seed_{:04d}'
                label = fmt.format (ind, n_trials, n_sig,
                                    poisson_str, s )
                commands.append (command)
                labels.append (label)
    sub.dry = dry
    if 'condor00' in hostname:
        sub.submit_condor00 (commands, labels) #, reqs=reqs)
    else:
        sub.submit_npx4 (commands, labels)

@cli.command ()
@click.option ('--n-trials', default=10000, type=int)
@click.option ('--gamma', default=2, type=float)
#@click.option ('--dec_deg',   default=0, type=float, help='Declination in deg')
@click.option ('--dec', 'dec_deg', multiple=True, type=float, default=[])
@click.option ('--dry/--nodry', default=False)
@click.option ('--seed', default=0)
@pass_state                                                                                                               
def submit_do_ps_sens (
        state, n_trials,  gamma,dec_deg,  dry, seed):
    ana_name = state.ana_name
    T = time.time ()
    job_basedir = state.job_basedir 
    job_dir = '{}/{}/ECAS_11yr/T_{}'.format (
        job_basedir, ana_name,  T)
    sub = Submitter (job_dir=job_dir, memory=5, 
        max_jobs=1000, config = submit_cfg_file)
    #env_shell = os.getenv ('I3_BUILD') + '/env-shell.sh'
    commands, labels = [], []
    this_script = os.path.abspath (__file__)
    trial_script = os.path.abspath('trials.py')
    
    sindecs = np.arange(-1,1.01,.1)
    sindecs[0] = -.99
    sindecs[-1] = .99
    dec_degs = dec_deg or np.degrees(np.arcsin(sindecs))
    for dec_deg in dec_degs:
        s =  seed
        fmt = '{} do-ps-sens  --n-trials {}' \
                            ' --gamma={:.3f} --dec_deg {}' \
                            ' --seed={}'
        command = fmt.format ( trial_script,  n_trials,
                              gamma, dec_deg, s)
        fmt = 'csky_sens_{:07d}_' \
                'gamma_{:.3f}_decdeg_{:04f}_seed_{:04d}'
        label = fmt.format (
                n_trials, 
                gamma, dec_deg, s)
        commands.append (command)
        labels.append (label)
    sub.dry = dry
    if 'condor00' in hostname:
        sub.submit_condor00 (commands, labels)
    else:
        sub.submit_npx4 (commands, labels)


@cli.command ()
@click.option ('--n-trials', default=10000, type=int)
@click.option ('--n-jobs', default=10, type=int)
@click.option ('-n', '--n-sig', 'n_sigs', multiple=True, default=[0], type=float)
@click.option ('--gamma', default=0, type=float)
@click.option ('-c', '--cutoff', default=np.inf, type=float, help='exponential cutoff energy (TeV)')      
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--sigsub/--nosigsub', default=True)
@click.option ('--catalog', type=str, default='seyfert_southernsky')
@click.option ('--seed', default=0)
@click.option ('--weightedfit', default=False, type=bool)
@click.option ('--corona', default=True, type=bool)
@click.option ('--debug', default=False, type=bool)
@click.option ('--nu_max', default=1000, type=float, help="cenA is 6.5")
@click.option ('--nu_min', default =0.1557, type=float)
@pass_state
def submit_do_seyfert_stacking_trials (
        state, n_trials, n_jobs, n_sigs, gamma, cutoff,  poisson,  sigsub,
        catalog,  seed, weightedfit, corona, debug, nu_min, nu_max):
    ana_name = state.ana_name
    ana_command=''
    if state.mask_deg !=0:
        ana_command +=' --mask_deg={}'.format(state.mask_deg)
    if state.source_r !=0:
        assert state.nsrc_tomask is not 0, "only mask top few sources!"
        ana_command +=' --nsrc_tomask={} --source_r={}'.format(state.nsrc_tomask, state.source_r)
        
    T = time.time ()
    poisson_str = 'poisson' if poisson else 'nopoisson'
    sigsub_str = 'sigsub' if sigsub else 'nosigsub'

    job_basedir = state.job_basedir 
    job_dir = '{}/{}/stacking_trials/T_E{}_{:17.6f}'.format (
        job_basedir, ana_name, int(gamma * 100),  T)
    sub = Submitter (job_dir=job_dir, memory=2, 
        max_jobs=1000, config = submit_cfg_file)
    commands, labels = [], []
    trial_script = os.path.abspath('trials.py')
    if catalog:
        catalogs = [catalog]
    else:
        catalogs = ['snr', 'unid', 'pwn']
    for cat in catalogs:
        for n_sig in n_sigs:
            for i in range (n_jobs):
                s = i + seed
                fmt = ' {} {} do-seyfert-stacking-trials --catalog={} --n-trials={}' \
                        ' --n-sig={} --gamma={:.3f} --cutoff={}' \
                        ' --{} --seed={} --weightedfit={} --corona={} --debug={} --nu_max={} --nu_min={} --{}'
  
                command = fmt.format (trial_script, ana_command,  cat, n_trials,
                                      n_sig, gamma, cutoff, poisson_str, s, weightedfit, corona, debug, nu_max, nu_min, sigsub_str)
                fmt = 'csky__cat_{}__trials_{:07d}__n_sig_{:08.3f}__' \
                        'gamma_{:.3f}_cutoff_{}_{}__seed_{:04d}'
                label = fmt.format (cat, n_trials, n_sig, gamma,
                                    cutoff, poisson_str, s)

                commands.append (command)
                labels.append (label)
    print(hostname)
    if 'condor00' in hostname:
        sub.submit_condor00 (commands, labels)
    else:
        sub.submit_npx4 (commands, labels)

@cli.command ()
@click.option ('--n-trials', default=1000, type=int)
@click.option ('--n-jobs', default=10, type=int)
@click.option ('--gamma', default=0, type=float)
@click.option ('-c', '--cutoff', default=np.inf, type=float, help='exponential cutoff energy (TeV)')
@click.option ('--dry/--nodry', default=False)
@click.option ('--catalog', type=str, default='seyfert_southernsky')
@click.option ('--seed', default=0)
@click.option ('--nsigma', default=0, type=float)
@pass_state
def submit_do_stacking_sens (
        state, n_trials, n_jobs,  gamma, cutoff,  dry,
        catalog,  seed, nsigma):
    ana_name = state.ana_name
    T = time.time ()
    job_basedir = state.job_basedir
    job_dir = '{}/{}/stacking_trials/T_E{}_{:17.6f}'.format (
        job_basedir, ana_name, int(gamma * 100),  T)
    sub = Submitter (job_dir=job_dir, memory=10,
        max_jobs=1000, config = submit_cfg_file)
    commands, labels = [], []
    trial_script = os.path.abspath('trials.py')
    if catalog:
        catalogs = [catalog]
    else:
        catalogs = ['snr', 'unid', 'pwn']
    for cat in catalogs:
            for i in range (n_jobs):
                s = i + seed
                fmt = ' {} do-stacking-sens --catalog={} --n-trials={}' \
                        ' --gamma={:.3f} --cutoff={}' \
                        ' --seed={}  '\
                        ' --nsigma={} '
                command = fmt.format (trial_script,  cat, n_trials,
                                      gamma, cutoff, s, nsigma)
                fmt = 'csky__cat_{}__trials_{:07d}__' \
                        'gamma_{:.3f}_cutoff_{}__seed_{:04d}_{}'
                if nsigma !=0:
                    l = 'dp'+str(nsigma)
                else:
                    l = 'sens'
                label = fmt.format (cat, n_trials, gamma,
                                    cutoff,  s, l)
                commands.append (command)
                labels.append (label)
    sub.dry = dry
    print(hostname)
    if 'condor00' in hostname:
        sub.submit_condor00 (commands, labels)
    else:
        sub.submit_npx4 (commands, labels)

@cli.command ()
@click.argument ('temp')
@click.option ('--n-trials', default=1000, type=int)
@click.option ('--n-jobs', default=10, type=int)
@click.option ('-n', '--n-sig', 'n_sigs', multiple=True, default=[0], type=float)
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--dry/--nodry', default=False)
@click.option ('-c', '--cutoff', default=np.inf, type=float, help='exponential cutoff energy (TeV)')
@click.option ('--seed', default=0, type=int)
@pass_state
def submit_do_gp_trials (
        state, temp, n_trials, n_jobs, n_sigs,
        poisson, dry, cutoff, seed):
    #example command using click python submit.py submit-do-gp-trials --n-sig=0 --n-jobs=1 --n-trials=1000 pi0
    ana_name = state.ana_name
    T = time.time ()
    job_basedir = state.job_basedir
    poisson_str = 'poisson' if poisson else 'nopoisson'
    job_dir = '{}/{}/gp_trials/{}/T_{:17.6f}'.format (
        job_basedir, ana_name, temp, T)
    sub = Submitter (job_dir=job_dir, memory=5,
        max_jobs=1000, config = submit_cfg_file)
    commands, labels = [], []
    reqs = '(Machine != "cobol93.private.pa.umd.edu")'
    trial_script = os.path.abspath('trials.py')
    print(n_sigs)
    for n_sig in n_sigs:
        for i in range (n_jobs):
            s = i + seed
            fmt = '{} do-gp-trials --n-trials={}' \
                    ' --n-sig={} ' \
                    ' --{} --seed={} --cutoff {} {}'
            command = fmt.format (trial_script,  n_trials,
                                  n_sig,  poisson_str,  s, cutoff, temp)
            fmt = 'csky__trials_{:07d}__n_sig_{:08.3f}__' \
                    '{}__{}__seed_{:04d}__cutoff_{}'
            label = fmt.format (
                    n_trials,  n_sig, temp, poisson_str,
                    s,   cutoff)
            commands.append (command)
            labels.append (label)
    if 'condor00' in hostname:
        print('submitting from condor00')
        sub.submit_condor00 (commands, labels, reqs=reqs)
    else:
        sub.submit_npx4 (commands, labels)

@cli.command ()
@click.argument('temp', default="pi0")
@click.option ('--n-trials', default=1000, type=int)
@click.option ('--n-jobs', default=10, type=int)
@click.option ('-n', '--n-sig', 'n_sigs', multiple=True, default=[0], type=float)
@click.option ('-nsrc', '--n-sources', 'n_srcs', default=14, type=int, help="the top n src to run")
@click.option ('--poisson/--nopoisson', default=True)
@click.option ('--dry/--nodry', default=False)
@click.option ('-c', '--cutoff', default=np.inf, type=float, help='exponential cutoff energy (TeV)')
@click.option ('--seed', default=0, type=int)
@click.option ('--gamma', default=0, type=float)
@pass_state
def submit_do_gp_bg_ps_trials (
        state, temp, n_trials, n_jobs, n_sigs, n_srcs,
        poisson, dry, cutoff, seed, gamma):
    #example command using click python submit.py submit-do-gp-trials --n-sig=0 --n-jobs=1 --n-trials=1000 pi0
    ana_name = state.ana_name
    T = time.time ()
    job_basedir = state.job_basedir
    poisson_str = 'poisson' if poisson else 'nopoisson'
    job_dir = '{}/{}/gp_trials/{}/T_E{}_{:17.6f}'.format (
        job_basedir, ana_name, temp, int(gamma * 100), T)
    sub = Submitter (job_dir=job_dir, memory=2,
        max_jobs=1000, config = submit_cfg_file)
    commands, labels = [], []
    trial_script = os.path.abspath('trials.py')

    for ind in range(n_srcs):
        for n_sig in n_sigs:
            for i in range (n_jobs):
                s = i + seed
                fmt = '{} do-gp-bg-ps-trials --n-trials={}' \
                        ' --n-sig={} ' \
                        ' --nsrc={} ' \
                        ' --{} --seed={} --cutoff {} --gamma={}'
                command = fmt.format (trial_script,  n_trials,
                                  n_sig, ind, poisson_str,  s, cutoff, gamma)
                fmt = 'csky__trials_{:07d}__n_sig_{:08.3f}__' \
                        'src_{}th__{}__{}__seed_{:04d}__cutoff_{}_gamma_{:08.3f}'
                label = fmt.format (
                    n_trials,  n_sig, ind, temp, poisson_str,
                    s,   cutoff, gamma)
                commands.append (command)
                labels.append (label)
    if 'condor00' in hostname:
        print('submitting from condor00')
        sub.submit_condor00 (commands, labels, reqs=reqs)
    else:
        sub.submit_npx4 (commands, labels)


if __name__ == '__main__':
    exe_t0 = now ()
    print ('start at {} .'.format (exe_t0))
    cli ()
