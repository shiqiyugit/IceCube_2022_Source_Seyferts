import h5py
import argparse
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.colors as colors
#import Plots

parser = argparse.ArgumentParser()
parser.add_argument("-f", type=str, nargs='+',
    dest="input_files",
    help="paths to input files (absolute or relative)")

args = parser.parse_args()
files = args.input_files

energies_all = []
zenith_all = []
azimuth_all = []
pulse_charge_all = []

for file in files:
    file_in = h5py.File(file, "r")

    features = file_in['features']
    labels = file_in['labels']

    dom_index = np.asarray(features['dom_index'])
    pulse_time = np.asarray(features['pulse_time'])
    pulse_charge = np.asarray(features['pulse_charge'])

    try:
        reco = file_in['reco']
    except:
        pass

    energies = np.asarray(labels['energy'])
    zenith = np.asarray(labels['zenith'])
    azimuth = np.asarray(labels['azimuth'])

    # this is NOT the most efficient way to do this, but anyway...
    for e in energies:
        energies_all.append(e)
    for z in zenith:
        zenith_all.append(z)
    for a in azimuth:
        azimuth_all.append(a)
    for p in pulse_charge:
        pulse_charge_all.append(p)

    file_in.close()


bins_e = np.logspace(3,8,10)
bins_z = np.arange(-1.0,1.2,.2)
bins_a = np.arange(0,2*np.pi+.3,.75)
bins_q = np.arange(0,500,100)

no_hits_all = np.zeros(len(pulse_charge_all))

if len(no_hits_all) != len(energies_all):
    raise RuntimeError("Length of energies (%i) and hits per DOM (%i) do not match"%(len(energies), len(no_hits_all)))

for i in range(len(no_hits_all)):
    no_hits_all[i] = len(pulse_charge_all[i])

fig, ax = plt.subplots()
fig.set_figheight(4.5)
fig.set_figwidth(7)

bins_hits = np.logspace(.01,6,7)

fig, ax = plt.subplots()
fig.set_figheight(4.5)
fig.set_figwidth(7)

bins_energy = np.logspace(2,8,20)
a, b, c, im = ax.hist2d(energies_all,no_hits_all, bins = [bins_energy,bins_hits],
    cmap=plt.cm.plasma, norm=colors.LogNorm())
ax.set_ylabel("Number of Hits")
ax.set_xlabel("Energy [GeV]")
cbar = plt.colorbar(im, ax=ax, label='Events')
ax.set_xscale('log')
ax.set_yscale('log')
fig.tight_layout()
fig.savefig('test_2d.png')
