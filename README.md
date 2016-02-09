# hcr16_ploscb
Code for generating figures in Henriksen, Cumming, &amp; Read (2016) 
PLoS Computational Biology. All code is released under GPL v2, a version
of which you can find on https://gnu.org/licenses/gpl-2.0.html. All simulations
require BEMtoolbox, which is available on https://github.com/sidh0/BEMtoolbox.
Clone the repository, run SetupBEMtoolbox in Matlab and run the code for the
figures. Note: most of these simulations will save data to disk in order to 
speed up simulations. This may take up to several GB, and will be stored in 
your BEMtoolbox/.data directory.

The code has not been extensively tested on different platforms and was 
developed on Linux so there might very well be some compatibility issues. 
Please let me know if you run into any errors 
(email: sid(dot)henriksen(at)gmail(dot)com). 

The organization of this repository is pretty straightforward. The directories
alternation_data and dotsize_data hold the psychophysical data for the
correlation alternation and dotsize experiments, respectively. The functions 
alternating_psych_analysis and dotsize_psych_analysis analyse these two
datasets.