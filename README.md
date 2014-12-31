exotimo
=======

Molecular orbital calculations including exotic particles such as positrons and muons. Emphasis will be on the following: 

1. computational approaches that will provide reasonable results for large molecular or solid state systems containing one or more exotic particles
2. methods for computing experimental observables such as positron annihilation lifetime spectroscopy and Doppler broadening of annihilation radiation

development approach
====================

*TODO*: describe use of Makefile.devel to generate tests from spec

ExotiMO is setup to use the `start_test` script that ships with Chapel. 
Create `alias testexotimo="$CHPL_HOME/util/start_test -no-chpl-home-warn -compopts '-M $EXOTIMO_HOME/src'"` and execute
`testexotimo` in a directory containing ExotiMO tests (i.e. `.chpl` and `.good` files).
