#!/bin/bash 

#run catelog trials with model inj model fit

#signal trials
#python submit.py submit-do-seyfert-ps-trials --nosigsub --n-trials 2000 --n-jobs 1 --dec=-43.0191 --dist=3.7 --logl=42.39 --dec=-65.339 --dist=6.0 --logl=42.63 --dec=-42.3706 --dist=20.9 --logl=43.48 --dec=-59.2348 --dist=39.2 --logl=44.09 --dec=-49.4682 --dist=8.1 --logl=42.07 --dec=-38.083 --dist=51.0 --logl=43.77 --dec=-44.1746 --dist=16.9 --logl=42.43 --dec=-39.9093 --dist=51.0 --logl=43.51 --dec=-37.7386 --dist=47.8 --logl=43.43 --dec=-30.9489 --dist=36.6 --logl=43.2 --dec=-30.3094 --dist=69.4 --logl=43.83 --dec=-17.2532 --dist=30.6 --logl=42.86 --dec=-22.8263 --dist=26.5 --logl=42.72 --dec=-5.3443 --dist=33.0 --logl=42.81 --dec=-26.0503 --dist=18.6 --logl=42.33 --dec=-36.1404 --dist=18.2 --logl=42.32 --dec=-65.4276 --dist=57.5 --logl=43.38 --dec=-29.134 --dist=93.2 --logl=43.87 --dec=-34.8537 --dist=46.2 --logl=43.12 --dec=-34.2956 --dist=33.2 --logl=42.74 --dec=-50.3026 --dist=158.0 --logl=44.44 --dec=-31.8698 --dist=33.9 --logl=42.76 --dec=-21.9278 --dist=70.7 --logl=43.5 --dec=-10.3235 --dist=22.8 --logl=42.31 --dec=-7.4562 --dist=35.6 --logl=42.69 --dec=-27.3078 --dist=109.1 --logl=43.92 --dec=-60.6178 --dist=55.8 --logl=43.25 --dec=-16.7286 --dist=71.6 --logl=43.41 --dec=-57.0688 --dist=49.3 --logl=43.08 --dec=-58.08 --dist=34.7 --logl=42.65  --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10  --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20
#bkg trails
#python submit.py submit-do-seyfert-ps-trials --n-trials 5000 --n-jobs 20 --dec=-43.0191 --dist=3.7 --logl=42.39 --dec=-65.339 --dist=6.0 --logl=42.63 --dec=-42.3706 --dist=20.9 --logl=43.48 --dec=-59.2348 --dist=39.2 --logl=44.09 --dec=-49.4682 --dist=8.1 --logl=42.07 --dec=-38.083 --dist=51.0 --logl=43.77 --dec=-44.1746 --dist=16.9 --logl=42.43 --dec=-39.9093 --dist=51.0 --logl=43.51 --dec=-37.7386 --dist=47.8 --logl=43.43 --dec=-30.9489 --dist=36.6 --logl=43.2 --dec=-30.3094 --dist=69.4 --logl=43.83 --dec=-17.2532 --dist=30.6 --logl=42.86 --dec=-22.8263 --dist=26.5 --logl=42.72 --dec=-5.3443 --dist=33.0 --logl=42.81 --dec=-26.0503 --dist=18.6 --logl=42.33 --dec=-36.1404 --dist=18.2 --logl=42.32 --dec=-65.4276 --dist=57.5 --logl=43.38 --dec=-29.134 --dist=93.2 --logl=43.87 --dec=-34.8537 --dist=46.2 --logl=43.12 --dec=-34.2956 --dist=33.2 --logl=42.74 --dec=-50.3026 --dist=158.0 --logl=44.44 --dec=-31.8698 --dist=33.9 --logl=42.76 --dec=-21.9278 --dist=70.7 --logl=43.5 --dec=-10.3235 --dist=22.8 --logl=42.31 --dec=-7.4562 --dist=35.6 --logl=42.69 --dec=-27.3078 --dist=109.1 --logl=43.92 --dec=-60.6178 --dist=55.8 --logl=43.25 --dec=-16.7286 --dist=71.6 --logl=43.41 --dec=-57.0688 --dist=49.3 --logl=43.08 --dec=-58.08 --dist=34.7 --logl=42.6 --nosigsub 

#run catelog trials with model inj pl fit
#signal trials
#python submit.py submit-do-seyfert-ps-trials --nosigsub --gamma 3 --n-trials 2000 --n-jobs 1 --dec=-43.0191 --dist=3.7 --logl=42.39 --dec=-65.339 --dist=6.0 --logl=42.63 --dec=-42.3706 --dist=20.9 --logl=43.48 --dec=-59.2348 --dist=39.2 --logl=44.09 --dec=-49.4682 --dist=8.1 --logl=42.07 --dec=-38.083 --dist=51.0 --logl=43.77 --dec=-44.1746 --dist=16.9 --logl=42.43 --dec=-39.9093 --dist=51.0 --logl=43.51 --dec=-37.7386 --dist=47.8 --logl=43.43 --dec=-30.9489 --dist=36.6 --logl=43.2 --dec=-30.3094 --dist=69.4 --logl=43.83 --dec=-17.2532 --dist=30.6 --logl=42.86 --dec=-22.8263 --dist=26.5 --logl=42.72 --dec=-5.3443 --dist=33.0 --logl=42.81 --dec=-26.0503 --dist=18.6 --logl=42.33 --dec=-36.1404 --dist=18.2 --logl=42.32 --dec=-65.4276 --dist=57.5 --logl=43.38 --dec=-29.134 --dist=93.2 --logl=43.87 --dec=-34.8537 --dist=46.2 --logl=43.12 --dec=-34.2956 --dist=33.2 --logl=42.74 --dec=-50.3026 --dist=158.0 --logl=44.44 --dec=-31.8698 --dist=33.9 --logl=42.76 --dec=-21.9278 --dist=70.7 --logl=43.5 --dec=-10.3235 --dist=22.8 --logl=42.31 --dec=-7.4562 --dist=35.6 --logl=42.69 --dec=-27.3078 --dist=109.1 --logl=43.92 --dec=-60.6178 --dist=55.8 --logl=43.25 --dec=-16.7286 --dist=71.6 --logl=43.41 --dec=-57.0688 --dist=49.3 --logl=43.08 --dec=-58.08 --dist=34.7 --logl=42.65  --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10 --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 

#bg trials


#STACKING

#run stacking with model inj. pl fit
#bg
#python submit.py submit-do-seyfert-stacking-trials --corona True --gamma 3.0 --weightedfit True --n-trials 2000 --n-jobs 50 --nosigsub 

#signal
#python submit.py submit-do-seyfert-stacking-trials --n-trials 1000 --n-jobs 10 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10 --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --gamma 3.0 --weightedfit True --corona True --nosigsub

#run stacking with model inj. model fit
#bg
#python submit.py submit-do-seyfert-stacking-trials --n-trials 2000 --n-jobs 50 --corona True --nosigsub
#signal 

#python submit.py submit-do-seyfert-stacking-trials --n-trials 1000 --n-jobs 10 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10 --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --corona True --nosigsub


#stacking wo cenA

#run stacking with model inj. pl fit top 10 wo cenA
#top 10 wo cenA
#python submit.py submit-do-seyfert-stacking-trials --n-trials 10000 --n-jobs 100 --corona True --nosigsub --nu_min 0.35 --nu_max 6.5
#top 5 wo cenA
#python submit.py submit-do-seyfert-stacking-trials --n-trials 10000 --n-jobs 100 --corona True --nosigsub --nu_min 0.61 --nu_max 6.5
#top 30 stacking bg without cenA
#python submit.py submit-do-seyfert-stacking-trials --n-trials 10000 --n-jobs 100 --corona True --nosigsub --nu_max 6.5

# catelog trials including GP in bg 
# top 30 sources by default
# inj model fit model 

# signal tials
#python submit.py submit-do-gp-bg-ps-trials --n-trials 2000 --n-jobs 1 --n-sig=1 --n-sig=2 --n-sig=3 --n-sig=4 --n-sig=5 --n-sig=6 --n-sig=7 --n-sig=8 --n-sig=9 --n-sig=10 --n-sig=11 --n-sig=12 --n-sig=13 --n-sig=14 --n-sig=15 --n-sig=16 --n-sig=17 --n-sig=18 --n-sig=19 --n-sig=20 --n-sig=21 --n-sig=22 --n-sig=23 --n-sig=24 --n-sig=25 --nosigsub
#bg
python submit.py submit-do-gp-bg-ps-trials --n-trials 10000 --n-jobs 1000 --n-sig=0 -nsrc 5

# ------------test before submission-------------
#python trials.py do-seyfert-ps-trials --n-trials 20 




# ----------------- test before submit ------------------
#python trials.py do-seyfert-ps-trials --nosigsub --n-trials 10 -dec=-43.0191 -dist=3.7 -logl=42.39  --n-sig=10

