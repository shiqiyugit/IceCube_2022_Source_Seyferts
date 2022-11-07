import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import h5py
import scipy.stats
import glob
import os
import sys
import math
import argparse
import time

from Generators import DataGenerator, SplitGenerator
from Attention import AttentionWithContext
from Plots import plot_uncertainty, plot_2dhist, plot_1dhist, plot_error, plot_loss, plot_error_vs_reco, plot_inputs, plot_outputs, plot_outputs_classify, plot_hit_info

import keras
import tensorflow as tf
from keras.utils.generic_utils import Progbar
from keras import backend as K, initializers, regularizers, constraints
from keras.engine.topology import Layer
from keras.models import Model
from keras.layers import Input, Dense, Concatenate, Embedding, BatchNormalization
from keras.optimizers import Adam

from keras.layers import Lambda, Flatten, Reshape, CuDNNLSTM, LSTM, Bidirectional, Activation, Dropout
from keras.layers import Conv1D, SpatialDropout1D, MaxPooling1D

np.set_printoptions(threshold=sys.maxsize)
from pone_RNN import *

def main(config=1):
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-h", "--hits",type=int,default=250, dest="hits", help="number of dom hits used for training")
    parser.add_argument("-e", "--epochs",type=int,default=8, dest="epochs", help="number of training epochs")
    parser.add_argument("-b", "--beta_1",type=float,default=0.8, dest="beta_1", help="learning rate decay parameter")
    parser.add_argument("-r", "--lr", type=float,default=0.0001, dest="lr", help="learning rate")
    parser.add_argument("-o", "--dropout", type=float,default=0., dest="dropout", help="change network dropout for each layer")
    parser.add_argument("-l", "--log_energy", type=int,default=0, dest="log_energy", help="use log energy rather than absolute for training")
    parser.add_argument("-f", "--file", type=str, default="pone_cleaned_pulses_linefit_10minDOM_min100000GeV.hdf5", dest="file_name", help="file to use for training")
    parser.add_argument("-p", "--path", type=str, default="/mnt/scratch/agarw132/", dest="path", help="path to input file")
    parser.add_argument("-u", "--output", type=str, default="./plots_new/", dest="output", help="output folder destination")
    parser.add_argument("-s", "--standardize", type=int,default=0, dest="standardize", help="perform data standardization")
    parser.add_argument("-c", "--checkpoints", type=int,default=0, dest="checkpoints", help="use training checkpoints from previous run")
    parser.add_argument("-w", "--weights", type=int,default=1, dest="weights", help="use sample weights for training")
    parser.add_argument("-t", "--data_type", type=str, default=None, help="name")
    parser.add_argument("-n", "--num_use", type=int,default=None, dest="num_use", help="number of samples to use for plotting")
    args = parser.parse_args()

    no_hits = args.hits
    no_epochs = args.epochs
    beta_1 = args.beta_1
    learning_rate = args.lr
    dropout = args.dropout
    use_log_energy = bool(args.log_energy)
    ff_name = args.path + args.file_name
    use_standardization = bool(args.standardize)
    use_checkpoints = bool(args.checkpoints)
    use_weights = bool(args.weights)
    num_use = args.num_use

    ff = h5py.File(ff_name, 'r')
    global gen_filename
    global save_folder_name
    gen_filename = "run_"+str(no_epochs)+"_test"
    save_folder_name = args.output + gen_filename + '/'
    if os.path.isdir(save_folder_name) != True:
        os.mkdir(save_folder_name)
    print("Saving to:", save_folder_name)

    reco = False
    if "reco" in ff.keys():
        reco = True

    network_labels = ["dir_x", "dir_y", "dir_z"]
    if use_standardization:
        normalization = normalize(ff, network_labels, use_log_energy)
    else:
        normalization = []
    gen = DataGenerator(ff, labels=network_labels, maxlen=no_hits, use_log_energy=use_log_energy,use_weights=use_weights,normal=normalization)
    gen_train = SplitGenerator(gen, fraction=0.70, offset=0.00)
    gen_val = SplitGenerator(gen, fraction=0.20, offset=0.70)
    gen_test = SplitGenerator(gen, fraction=0.10, offset=0.90)

    vocab_size = VOCAB_SIZE
    time_samples = no_hits

    print("Plotting input distributions")
    t_inputs_start = time.time()
    #pulse_time_data, pulse_charge_data
    plot_inputs(ff["features/pulse_time"][:], ff["features/pulse_charge"][:], num_use=num_use, log_charge=False, gen_filename=save_folder_name)
    t_inputs_end = time.time()
    print((t_inputs_end-t_inputs_start)/60., "minutes to plot inputs")

    labels_raw = None
    labels_predicted_raw = None
    if reco: labels_reco = None
    weights_raw = None
   
    print("Testing model")
    for i in range(len(gen_train)-1):
        batch_features, batch_labels, batch_weights = gen_train[i]

        if labels_raw is None:
            labels_raw = batch_labels
            weights_raw = batch_weights
        else:
            labels_raw           = np.append(labels_raw,           batch_labels,           axis=0)
            weights_raw          = np.append(weights_raw,          batch_weights,          axis=0)
        del batch_features
        del batch_labels
        del batch_weights

    train_labels = labels_raw#gen.untransform_labels(labels_raw)
    train_weights = weights_raw

    for i in range(len(gen_val)-1):
        batch_features, batch_labels, batch_weights = gen_val[i]

        if labels_raw is None:
            labels_raw = batch_labels
            weights_raw = batch_weights
        else:
            labels_raw           = np.append(labels_raw,           batch_labels,           axis=0)
            weights_raw          = np.append(weights_raw,          batch_weights,          axis=0)
        del batch_features
        del batch_labels
        del batch_weights

    val_labels = gen.untransform_labels(labels_raw)
    val_weights = weights_raw

    for i in range(len(gen_test)-1):
        batch_features, batch_labels, batch_weights = gen_test[i]

        if labels_raw is None:
            labels_raw = batch_labels
            weights_raw = batch_weights
        else:
            labels_raw           = np.append(labels_raw,           batch_labels,           axis=0)
            weights_raw          = np.append(weights_raw,          batch_weights,          axis=0)
        del batch_features
        del batch_labels
        del batch_weights
 
    test_labels = gen.untransform_labels(labels_raw)
    test_weights = weights_raw

    labels = np.concatenate((train_labels, val_labels, test_labels))
    weights = np.concatenate((train_weights, val_weights, test_weights))

    dx_true = labels[:, 0]
    dy_true = labels[:,1]
    dz_true = labels[:,2]
    total_entries = len(weights)
    #shuffle entries
    order = np.arange(total_entries)
    np.random.seed(86)
    np.random.shuffle(order)
    weights = weights[order]

    dx_true = dx_true[order]
    dy_true = dy_true[order]
    dz_true = dz_true[order]

    from scipy.stats import norm

    #isTrack_predicted = labels_predicted[:,4]
    #isCascade_predicted = labels_predicted[:,5]
    #isTrack_true = labels[:,4]
    #isCascade_true = labels[:,5]
    
    #isTrack_predicted = [isTrack_predicted > isCascade_predicted]
    #isCascade_predicted = [isCascade_predicted > isTrack_predicted]
    
    #trueTracks = np.sum(np.logical_and(isTrack_true, isTrack_predicted))
    #falseTracks = np.sum(np.logical_and(np.logical_not(isTrack_true), isTrack_predicted))
    #trueCascades = np.sum(np.logical_and(isCascade_true, isCascade_predicted))
    #falseCascades = np.sum(np.logical_and(np.logical_not(isCascade_true), isCascade_predicted))
    
    #fig, ax = plt.subplots()
    #bars1 = ax.bar(np.arange(2), [trueTracks, trueCascades], 0.25, color="SkyBlue")
    #bars2 = ax.bar(np.arange(2)+0.5*np.ones(2), [falseTracks, falseCascades], 0.25, color="IndianRed")
    #ax.set_title("Track vs. Cascade classification results")
    #ax.set_xticks(np.arange(4)/2)
    #ax.set_xticklabels(("True Tracks", "False Tracks", "True Cascades", "False Cascades"))
    #imgname = save_folder_name+"class.png"
    #plt.savefig(imgname)
    
    zenith_true, azimuth_true = np.degrees(to_zenazi(dx_true, dy_true, dz_true))

    #Make plots
    print("Plotting regression output distributions")
    t_regression_start = time.time()

    plot_outputs(dx_true, -1.0, 1.0, "dx [m]", weights, num_use=num_use, logscale=False, gen_filename=save_folder_name)
    plot_outputs(dy_true, -1.0, 1.0, "dy [m]", weights, num_use=num_use, logscale=False, gen_filename=save_folder_name)
    plot_outputs(dz_true, -1.0, 1.0, "dz [m]", weights, num_use=num_use, logscale=False, gen_filename=save_folder_name)
    plot_outputs(azimuth_true, 0, 360, "Azimuth [degrees]", weights, num_use=num_use, logscale=False, gen_filename=save_folder_name)
    plot_outputs(zenith_true, 0, 180, "Zenith [degrees]", weights, num_use=num_use, logscale=False, gen_filename=save_folder_name)
    t_regression_end = time.time()
    print((t_regression_end-t_regression_start)/60., "minutes to plot regression outputs")

    print("Plotting classification output distributions")
    t_classification_start = time.time()
    t_classification_end = time.time()
    print((t_classification_end-t_classification_start)/60., "minutes to plot classification outputs")

#    print("Plotting hit information")
#    t_hit_start = time.time()
#    plot_hit_info(ff, energy_true, order, num_use=num_use, logscale=False, gen_filename=save_folder_name)
#    t_hit_end = time.time()
#    print((t_hit_end-t_hit_start)/60., "minutes to plot hit information")

    return 0#network_history.history['val_loss']

if __name__ == "__main__":
    main()
