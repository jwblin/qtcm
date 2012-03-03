README FILE
By Johnny Lin


This directory consists of output from runs of the pure-Fortran
version of QTCM.  This output can be considered the "correct"
versions of output to validate against.

/create         Directory of files needed to create the
		benchmarks.

/aquaplanet     Default aquaplanet run (landon set to 0) for 15
                days.

/landon         Default landon run (landon set to 1) for 15 days.

/aquaplanet_*   A variety of other aquaplanet runs, with parameter
                changes.

/landon_*       A variety of other landon runs, with parameter
                changes.

/timing         Script to run an aquaplanet benchmark run for 365
		days, for the purposes of timing comparison as in
		the introduction chapter of the manual.
