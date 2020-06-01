## mERG data processing scripts for data recorded on MED64 system

* **analyse_mERG_MED64.R** is used to process each recording from a sample, the input is a csv file containing the average of 3 mERG recordings. The data is passed through a 1-50Hz band pass filter and peaks are detected using **analyse_merg_peaks.R**. It outputs, summary graphs as well as individual channel traces with detected peaks indicated. 

* **analyse_merg_peaks.R** contain functions to detect peaks from mERG recordings. It applays a smoothing function and detects all the maxima or minima and return their positions and values. 

* **summarise_peaks.R** aggregates all of the results processed with **analyse_mERG_MED64.R**. Experimental conditions (solutions, stimuli, graft type, time after transplantation) is obtained from the file path. 

* **bernoulli.stan** is the stan script used to analyse the probability of b-wave. 

