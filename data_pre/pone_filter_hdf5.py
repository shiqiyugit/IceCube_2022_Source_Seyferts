#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/icetray-start
#METAPROJECT: simulation/V06-01-00-RC4

import numpy as np
import h5py
import glob
import sys
import math
import argparse
np.set_printoptions(threshold=sys.maxsize)

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", type=str, default=None,
                    dest="input_file", help="name of input hdf5 file")
parser.add_argument("--base_name", type=str, default=None,
                    dest="base_name", help="base of output file name")
parser.add_argument("--min_DOM", type=int, default=10,
                    dest="min_DOM", help="minimum number of DOMs required for event to pass")
parser.add_argument("--energy", type=float, default=100000,
                    dest="energy", help="minimum energy")
parser.add_argument("--only_reco", type=bool, default=False,
                    dest="only_reco", help="include only reco events in sample")
args = parser.parse_args()

min_DOM = args.min_DOM # makes optional cut of events that only hit a certain number of DOMs
only_reco = args.only_reco
energy = args.energy

print("Printing cuts...")
if min_DOM > 0: print("Minimum number of DOMs hit: %i"%min_DOM)
else: print("No cut on minimum number of DOMs hit")
if only_reco == True: print("Only including reco events")
else: print("Keeping events with and without PegLeg reconstructions")
if energy>0: print("Minimum number of energy: ", energy)

fin = h5py.File(args.input_file, 'r')

label_keys = list(fin["labels"].keys())
print(label_keys)
feature_keys = list(fin["features"].keys())
if "reco" in fin:
    reco_keys = list(fin["reco"].keys())

features = dict()
labels = dict()
if "reco" in fin:
    reco = dict()

for k in label_keys:
    labels[k] = np.array(fin["labels"][k])
for k in feature_keys:
    features[k] = np.array(fin["features"][k])
if "reco" in fin:
    for k in reco_keys:
        reco[k] = np.array(fin["reco"][k])
weights = np.array(fin["weights"])

total_entries = len(weights)

#output_file = args.base_name
import os
thisFile = os.path.split(args.input_file)[1]
output_file = os.path.splitext(thisFile)[0]

if min_DOM:
    output_file += "_%iminDOM"%min_DOM
if energy>0:
    output_file += "_min%iGeV"%energy

fout = h5py.File("/mnt/scratch/agarw132/"+output_file+".hdf5", 'w')

#shuffle entries
order = np.arange(total_entries)
np.random.shuffle(order)
for k in label_keys:
    try:
        labels[k] = labels[k][order]
    except:
        continue
for k in feature_keys:
    features[k] = features[k][order]
if "reco" in fin:
    for k in reco_keys:
        reco[k] = reco[k][order]
weights = weights[order]

print("Finished reading file")

#------------------------------------------------------------

if min_DOM != 0:
    print("Masking for events that have at least %i hit DOMs"%min_DOM)
    dom_mask = []
    for i in range(len(features["dom_index"])):
        dom_ids = features["dom_index"][i]
        unique_doms = np.unique(dom_ids)
        if len(unique_doms) >= min_DOM:
            dom_mask.append(True)
        else:
            dom_mask.append(False)

    for k in label_keys:
        try:
            labels[k] = labels[k][dom_mask]
        except:
            continue
    for k in feature_keys:
        if k != "dom_index": features[k] = features[k][dom_mask]
    if "reco" in fin:
        for k in reco_keys:
            reco[k] = reco[k][dom_mask]
    weights = weights[dom_mask]
    features["dom_index"] = features["dom_index"][dom_mask]

if only_reco == True and "reco" in fin:
    print("Only keeping events with reconstructions")
    boolarray = np.array(reco["zenith"]) != 0
    for k in label_keys:
        labels[k] = labels[k][boolarray]
    for k in feature_keys:
        features[k] = features[k][boolarray]
    for k in reco_keys:
        if k != "energy":
            reco[k] = reco[k][boolarray]
    weights = weights[boolarray]
    reco["zenith"] = reco["zenith"][boolarray]

if energy>0:
    print("Only keeping high energy events")
    boolarray = np.array(labels["energy"])>energy
    for k in label_keys:
        labels[k] = labels[k][boolarray]
    for k in feature_keys:
        features[k] = features[k][boolarray]
    for k in reco_keys:
        if k != "energy":
            reco[k] = reco[k][boolarray]
    weights = weights[boolarray]
    reco["energy"] = reco["energy"][boolarray]

#------------------------------------------------------------

#shuffle entries
remaining_entries = len(weights)
order = np.arange(remaining_entries)
np.random.shuffle(order)
for k in label_keys:
    try:
        labels[k] = labels[k][order]
    except:
        continue
for k in feature_keys:
    features[k] = features[k][order]
if "reco" in fin:
    for k in reco_keys:
        reco[k] = reco[k][order]
weights = weights[order]

print("Finished making cuts")

print("Total events before cuts:", total_entries)
print("Total events after cuts:", remaining_entries)

grp_features = fout.create_group("features")
grp_labels   = fout.create_group("labels")
if "reco" in fin:
    grp_reco = fout.create_group("reco")
grp_weights = fout.create_dataset("weights", data=weights)

for k in label_keys:
    grp_labels.create_dataset(k, data=labels[k])
for k in feature_keys:
    grp_features.create_dataset(k, data=features[k])
if "reco" in fin:
    for k in reco_keys:
        grp_reco.create_dataset(k, data=reco[k])

fin.close()
fout.close()
