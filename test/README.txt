README FILE
By Johnny Lin


This directory contains example Python scripts and unit tests for the
package.  Unit tests mainly use the unittest module.

./benchmarks    The "correct" versions of output to validate against.

./bnddir        Directory of boundary data used by the test scripts.

./rundir        A directory with run space the test scripts will use
                to run in.

create_*.py     Create various qtcm runs that will later be analyzed
                by the test_*.py scripts.

run_long.py     Executes all the long create_*.py files.

run_short.py    Executes all the short create_*.py files.

timing_365.py   Executes and calculates full and parts Qtcm timings.

test_all.py     Runs all tests on the package, if run as main, i.e.
                at the command line you execute:

                $ python test_all.py

		This script just calls all the other test_*.py
		scripts in this directory, so it needs to be run
		in the same directory of those scripts.  This script
		also executes all the run_*.py scripts, so you don't
		have to execute create_*.py or run_*.py scripts
		independently.

test_*.py       Besides test_all.py, test scripts to test various
                aspects of the package.
