#!/bin/bash 
#-------------------------------------------------------------------
#do catelog trials
#-------------------------------------------------------------------

#model inj model fit
#submit_do_gp_bg_ps_trials
#signal trials
python submit.py --mask_deg 10 --nsrc_tomask 4 --source_r 7 --mask_self True submit-do-gp-bg-ps-trials --n-trials 10000 --n-jobs 1 --n-sig=1 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10  --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20 -nsrc 14

#bkg trails
python submit.py --mask_deg 10 --nsrc_tomask 4 --source_r 7 --mask_self True submit-do-gp-bg-ps-trials --n-trials 2000000 --n-jobs 50 -nsrc 14

#model inj pl fit
#signal trials
#python submit.py --mask_deg 10 --nsrc_tomask 4 --source_r 7 --mask_self True submit-do-gp-bg-ps-trials --n-trials 10000 --n-jobs 1 --n-sig=1 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10  --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20 -nsrc 14 --gamma 3
#python submit.py --mask_deg 10 --nsrc_tomask 4 --source_r 7 --mask_self True submit-do-seyfert-ps-trials --osigsub --gamma 3 --n-trials 10000 --n-jobs 1 --n-sig=1 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10  --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20 -nsrc 14

#bg trials
#python submit.py --mask_deg 10 --nsrc_tomask 4 --source_r 7 --mask_self True submit-do-gp-bg-ps-trials --n-trials 2000000 --n-jobs 50 -nsrc 14 --gamma 3
#python submit.py --mask_deg 10 --nsrc_tomask 4 --source_r 7 --mask_self True submit-do-seyfert-ps-trials --n-trials 5000000 --n-jobs 20 --nosigsub --gamma 3 -nsrc 14

#-------------------------------------------------------------------
#STACKING
#-------------------------------------------------------------------

#RUN  STACKING: MODEL INJ PL FIT
#bg
# with  cenA, top 30
#python submit.py --mask_deg 1.5 --nsrc_tomask 4 --source_r 7 submit-do-seyfert-stacking-trials --corona True --gamma 3 --weightedfit True --n-trials 500000 --n-jobs 200 --nosigsub --nu_min 0.3
# without cenA, top 30-1
#python submit.py --mask_deg 1.5 --nsrc_tomask 4 --source_r 7 submit-do-seyfert-stacking-trials --corona True --gamma 3 --weightedfit True --n-trials 500000 --n-jobs 200 --nosigsub --nu_max 6 --nu_min 0.3

#signal
#with  cenA top 30
#python submit.py --mask_deg 1.5 --nsrc_tomask 4 --source_r 7 submit-do-seyfert-stacking-trials --n-trials 10000 --n-jobs 1 --n-sig=1 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10 --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20 --n-sig=21 --n-sig=22 --n-sig=23 --gamma 3 --weightedfit True --corona True --nosigsub --nu_min 0.3

# wo cenA 30-1
#python submit.py --mask_deg 1.5 --nsrc_tomask 4 --source_r 7 submit-do-seyfert-stacking-trials --n-trials 10000 --n-jobs 1 --n-sig=1 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10 --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20 --n-sig=21 --n-sig=22 --n-sig=23 --gamma 3 --weightedfit True --corona True --nosigsub --nu_max 6 --nu_min 0.3

#RUN STACKING TRIALS: MODEL INJ MODEL FIT
#bg
#with cen
#python submit.py --mask_deg 1.5 --nsrc_tomask 4 --source_r 7 submit-do-seyfert-stacking-trials --n-trials 50000 --n-jobs 200 --corona True --nosigsub --nu_min 0.3
#wo cenA
#python submit.py --mask_deg 1.5 --nsrc_tomask 4 --source_r 7 submit-do-seyfert-stacking-trials --n-trials 50000 --n-jobs 200 --corona True --nosigsub --nu_max 6 --nu_min 0.3

#signal 
#with cenA
#python submit.py --mask_deg 1.5 --nsrc_tomask 4 --source_r 7 submit-do-seyfert-stacking-trials --n-trials 10000 --n-jobs 1 --n-sig=1 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10 --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20 --n-sig=21 --n-sig=22 --n-sig=23 --corona True --nosigsub --nu_min 0.3
#wo cenA
#python submit.py --mask_deg 1.5 --nsrc_tomask 4 --source_r 7 submit-do-seyfert-stacking-trials --n-trials 10000 --n-jobs 1 --n-sig=1 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10 --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20 --n-sig=21 --n-sig=22 --n-sig=23 --corona True --nosigsub --nu_max 6 --nu_min 0.3

# ----------------- test before submit ------------------
#python trials.py --mask_degree 1.5 do-seyfert-ps-trials --sigsub --n-trials 10 -nsrc 1
#python trials.py do-seyfert-stacking-trials --n-trials 1000 --debug True
#python trials.py --mask_deg 1.5 do-stacking-sens --nsigma 0 --n-trials 3000
