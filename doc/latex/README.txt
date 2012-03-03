README FILE
By Johnny Lin


This directory contains LaTeX documentation for the package.

/archive     Older versions of various files for the manual.

/figs        Directory of figures for the manual.

README.txt   This file.

code_to_latex.py  A Python program to create tables of various
                  model values using the package code itself.

index.html   Redirect to ../html/index.html.

manual.tex   Main LaTeX file for the manual, with the preamble.
	     Note that this file gives a call to my master BibTeX
	     file, which is not included in the distribution.

manual.*     Besides manual.tex, these are the auxilary files
	     leftover from the last compile of manual.tex.  The
	     .bbl file is here for users who may want to do their
	     own compile of the manual (you'll need it since you
	     won't have my master BibTeX database to run BibTeX
	     with).  The .aux file is used by latex2html.

*.tex        Besides manual.tex, all other LaTeX files that are
             called by manual.tex.

A copy of the license this package is released under is found in
the file ../GNU_LGPL.txt.
