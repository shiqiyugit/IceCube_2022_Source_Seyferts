#!/bin/bash 

#python trials.py find-stacking-n-sig --inputname 30_corona_powerlaw_weightedfit
python trials.py find-stacking-n-sig --inputname 30_corona_powerlaw_weightedfit --nsigma 3
python trials.py find-stacking-n-sig --inputname 30_corona_powerlaw_weightedfit --nsigma 5

#python trials.py find-stacking-n-sig --inputname 30_corona_flux 
python trials.py find-stacking-n-sig --inputname 30_corona_flux --nsigma 3
python trials.py find-stacking-n-sig --inputname 30_corona_flux --nsigma 5

python trials.py find-ps-n-sig --inputname 30_corona_powerlaw --nsigma 5
python trials.py find-ps-n-sig --inputname 30_corona_powerlaw --nsigma 3 
#python trials.py find-ps-n-sig --inputname 30_corona_powerlaw 

python trials.py find-ps-n-sig --inputname 30_corona_flux  --nsigma 5
python trials.py find-ps-n-sig --inputname 30_corona_flux  --nsigma 3
#python trials.py find-ps-n-sig --inputname 30_corona_flux 

# ------------test before submission-------------
#python trials.py do-seyfert-ps-trials --n-trials 20 
