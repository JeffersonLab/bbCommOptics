This is the file directory for optimizing BB and uses the original BB sieve design. Note that yTar and zVertex are opposite in sign for BB on beam left. 

Add the run details in list_of_optics_run.dat. Here, I have a dummy run number "100", and the sim file is hardcoded in make_hist_bb_optics.C. The central angle of BB, number of files, and number of additional cuts is also an input. For now, the additional cut is wide and available, but not used. 

The following steps can be used to run and analyze the optics data:
a. make_hist_bb_optics.C(100,kFALSE,kFALSE,kFALSE,-1)
b. set_ytar_delta_cuts.C(100,-1)
c. make_hist_bb_optics.C(100,kTRUE,kFALSE,kFALSE,-1)
d. set_ypfp_yfp_cuts.C(100,-1)
e. plot_yfp_cuts.C(100,-1)
f. set_xpfp_xfp_cuts.C(100,-1)
g. plot_xfp_cuts.C(100,-1)
h. make_fit_ntuple.C(100,-1)
i. plot_yptar_diff.C(100,-1), plot_ytar_diff.C(100,-1), plot_xptar_diff.C(100,-1)
j. fit_opt_matrix_v2.C()

When running step (j), one should check the input and output matrix elements files. Also, do not include the 0000 order input ME, but these will be re-fit.

Make directories:
hist
cuts
plots

Still TODO:
-include mis-pointing
-include beam offsets

