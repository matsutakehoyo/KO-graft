# Retinal Organoid Fluorescence Analysis

Analysis of fluorescence signal increase in stem cell derived retinal organoids. The fluorescence signal increase across DD (diferentiation day) is modeled using the modified Gompertz curve, which has three main parameters. 

- **A**: maximum fluorescence intensity
- **lambda**: The delay of onset, e.i., the begining of fluorescence increa. 
- **mu**: the maximum fluorescence incerase rate

- `get_exif.sh` ia a bash script used to extract image metadata and save it to a text file. 
- The following R srcitps were used to prerpocess and model the fluorescence images of retinal organoids (auxillary files not included).
  - `fluorescence analysis.R`: extract the channel data (RGB values) from images and stores it in a dataframe together with the above mentioned metadata. 
  - `fluorescence analyisis modeling.R`: process the channel information to subtract background and model the signal.

- The Stan script (`growth_curve_nested.stan`) shows the code for the Bayesina inference using Rstan. 
