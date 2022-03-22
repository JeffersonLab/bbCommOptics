Using straight-track replayed data, we can obtain the alignment of the GEMs relative to the spectrometer coordinate system. We will determine the GEM X, Y, and Z origin of the front layer, the Z position of the sieve, and the GEM theta and phi orientation.

Begin by modifying the config.dat to include any physics cuts and the relevant run files.

The code is run in ROOT as:
.x AlignZeroField.C("config.dat")

Currently the code allows the user to make cuts on the x and y sieve hole positions using the convention documented on page 2 here:
https://sbs.jlab.org/cgi-bin/DocDB/private/ShowDocument?docid=166

This could be improved by fitting the peaks, and using the cuts to assign the sieve events in a particular hole to that hole/position.

One must enter the minimized values in the second iteration of the code to see the results on the sieve distribution.
