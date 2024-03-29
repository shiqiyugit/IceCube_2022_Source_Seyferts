#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3.0.1/icetray-start
#METAPROJECT: simulation/V06-01-00-RC4

import numpy as np
import h5py
import glob
import argparse
import sys

def get_file_entries(filename):
    f = h5py.File(filename, 'r')
    return len(f["weights"])
    f.close()
    del f

def check_reco(filename):
    f = h5py.File(filename, 'r')
    return ("reco" in f.keys())
    f.close()
    del f

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_files", type=str, default=None,
                    dest="input_files", help="name for input files")
parser.add_argument("-o", "--output_file", type=str, default=None,
                    dest="output_file", help="name for output file")
args = parser.parse_args()

filenames = sorted(glob.glob("/mnt/scratch/agarw132/pone_50TeV_to_1PeV/*.hdf5"))
outfilename = "/mnt/scratch/agarw132/pone_cleaned_pulses_linefit.hdf5"
#
#filenames = sorted(glob.glob(args.input_files))
#outfilename = args.output_file

np.set_printoptions(threshold=sys.maxsize)

reco = False
entries_per_file = []
for filename in filenames:
    print("checking {}".format(filename))
    entries_per_file.append(get_file_entries(filename))
    if not reco:
        reco = check_reco(filename)
total_entries = len(entries_per_file)
print("total entries: {}".format(total_entries))

f = h5py.File(outfilename, 'w')
grp_features = f.create_group("features")
grp_labels   = f.create_group("labels")
if reco:
    grp_reco = f.create_group("reco")

# create arrays from template
f_template = h5py.File(filenames[0], 'r')

label_keys = list(f_template["labels"].keys())
feature_keys = list(f_template["features"].keys())
if "reco" in f_template:
    reco_keys = list(f_template["reco"].keys())

out_features = dict()
out_labels = dict()
if "reco" in f_template:
    out_reco = dict()

for k in label_keys:
    out_labels[k] = np.array([])
for k in feature_keys:
    out_features[k] = np.array([])
if "reco" in f_template:
    for k in reco_keys:
        out_reco[k] = np.array([])
out_weights = np.array([])

f_template.close()
del f_template

# get minimum number of entries in files
min_entries = min(entries_per_file)
print("Minimum number of entries in all files:", min_entries)

# read data from input files
current_entry = 0
file_zeniths = []
file_time_lists = []
for filename in filenames:
    f_input = h5py.File(filename, 'r')
    entries = len(f_input["weights"])

    if entries == 0:
        continue

    # generate random numbers for checking
    save_length = np.random.randint(5,10)
    save_index = np.random.randint(0,entries-(save_length+1))

    print("reading input file {} with {} entries".format(filename, entries))

    out_weights = np.concatenate((out_weights, f_input["weights"][:]))
    for k in label_keys:
        out_labels[k] = np.concatenate((out_labels[k], f_input["labels"][k][:]))
    for k in feature_keys:
        out_features[k] = np.concatenate((out_features[k], f_input["features"][k][:]))
    if "reco" in f_input:
        for k in reco_keys:
            out_reco[k] = np.concatenate((out_reco[k], f_input["reco"][k][:]))

    test_zeniths = np.array(f_input["labels"]["zenith"][save_index:save_index+save_length])
    file_zeniths.append(test_zeniths)
    test_times = np.array(f_input["features"]["pulse_time"][save_index])
    file_time_lists.append(test_times)
    
    f_input.close()
    del f_input

    current_entry += entries

if current_entry != len(out_weights):
    raise RuntimeError("Unexpected number of events -- got %i, expected %i"%(current_entry, len(out_weights)))

for k in label_keys:
    grp_labels.create_dataset(k, data=out_labels[k])
print("Finished creating labels")
for k in feature_keys:
    dt = h5py.special_dtype(vlen=out_features[k][0].dtype)
    dset = grp_features.create_dataset(k, (len(out_features[k]), ), dtype=dt)
    for i in range(len(out_features[k])):
        dset[i] = out_features[k][i]
print("Finished creating features")
if reco:
    for k in reco_keys:
        grp_reco.create_dataset(k, data=out_reco[k])
print("Finished creating reco")
grp_weights = f.create_dataset("weights", data=out_weights)
print("Finished creating weights")

f.close()
print(' ')

# check for concatenation errors
def isSubArray(long_array, short_array):
    i = 0
    j = 0
    m = len(long_array)
    n = len(short_array)

    max_match = 0
    match = 0

    if type(long_array[0]) != np.ndarray:
        while i < m and j < n:
            if long_array[i] == short_array[j]:
                i += 1
                j += 1
                match += 1
                if match > max_match:
                    max_match = int(match)

                if j == n:
                    return True

            else:
                i = i - j + 1
                j = 0
                match = 0

            if i == m:
                print("Longest match:", max_match)
                return False

    else:
        done_return = False
        for k, sub_array in enumerate(long_array):
            p = len(sub_array)
            while i < p and j < n:
                if sub_array[i] == short_array[j]:
                    i += 1
                    j += 1
                    match += 1
                    if match > max_match:
                        max_match = int(match)

                    if j == n:
                        done_return = True
                        return True

                else:
                    if match > 0:
                        print(match)
                    i = 0
                    j = 0
                    match = 0
                    break 
            
        if not done_return:
            print("Longest match:", max_match)
            return False

bool_array = []
for i, filename in enumerate(filenames):
    f_input = h5py.File(filename, 'r')
    zeniths = f_input["labels"]["zenith"][:]
    bool_array.append(isSubArray(zeniths, file_zeniths[i]))
    f_input.close()
    del f_input

for i in range(len(bool_array)):
    print("Kept regression labels from %s: %s"%(filenames[i], bool_array[i]))
if False in bool_array:
    raise RuntimeError("Regression label information not kept from all files")

for i in range(len(bool_array)):
    print("Kept classification labels from %s: %s"%(filenames[i], bool_array[i]))
if False in bool_array:
    raise RuntimeError("Classification label information not kept from all files")

bool_array = []
for i, filename in enumerate(filenames):
    f_input = h5py.File(filename, 'r')
    time_lists = f_input["features"]["pulse_time"][:]
    bool_array.append(isSubArray(time_lists, file_time_lists[i]))
    f_input.close()
    del f_input

for i in range(len(bool_array)):
    print("Kept features from %s: %s"%(filenames[i], bool_array[i]))
if False in bool_array:
    raise RuntimeError("Feature information not kept from all files")

f = h5py.File(outfilename, "r")
print(f.keys())
out_zeniths = f["labels"]["zenith"][:]
out_time_lists = f["features"]["pulse_time"][:]

bool_array = []
for i, filename in enumerate(filenames):
    bool_array.append(isSubArray(out_zeniths, file_zeniths[i]))
for i in range(len(bool_array)):
    print("Found outfile regression labels from %s: %s"%(filenames[i], bool_array[i]))
if False in bool_array:
    raise RuntimeError("Regression label information not found in outfile")

bool_array = []
for i in range(len(bool_array)):
    print("Found outfile classification labels from %s: %s"%(filenames[i], bool_array[i]))
if False in bool_array:
    raise RuntimeError("Classification label information not found in outfile")

bool_array = []
for i, filename in enumerate(filenames):
    bool_array.append(isSubArray(out_time_lists, file_time_lists[i]))
for i in range(len(bool_array)):
    print("Found outfile features from %s: %s"%(filenames[i], bool_array[i]))
if False in bool_array:
    raise RuntimeError("Feature information not found in outfile")

f.close()
del f
