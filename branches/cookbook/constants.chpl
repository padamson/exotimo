/*
 constants.chpl: Constants for ExotiMO

 This program is part of the ExotiMO quantum chemistry program suite.

 Copyright (c) 2014, Paul E. Adamson. All Rights Reserved. 

 ExotiMO version XXX and later is covered by the XXX
 license. Please see the file LICENSE that is part of this
 distribution. 
*/

// Misc
const PI: real = 3.14159265358979323846;
const ZERO: real = 0.0e0;
const HALF: real = 0.5e0;
const ONE: real = 1.0e0;
const TWO: real = 2.0e0;
const crit:real = 1.0E-6;

// Misc units
const Clight: real=2.99792458e8;     // speed of light in m/s
const Kboltz: real=3.166830e-6;      // Boltzmann constant
const e2: real = 14.399;             // Coulomb's law coeff if R in \AA and resulting E in eV
const planck: real=6.6260755e-34;    // Planck's constant, in Js

// Distance units
const bohr2ang: real = 0.529177249;  // Conversion of length from bohr to angstrom
const ang2bohr: real = 1/bohr2ang;

// Energy units
const hartree2kcal: real = 627.5095; // Hartree to kcal/mol conversion
const kcal2hartree: real = 1/hartree2kcal;

const ev2kcal: real = 23.061;        // Conversion of energy in eV to energy in kcal/mol
const kcal2ev: real = 1/ev2kcal;

const hartree2joule: real = 4.3597482e-18;   // Hatree to Joule conversion factor
const joule2hartree: real = 1/hartree2joule;

const ev2hartree: real = hartree2kcal/ev2kcal;
const hartree2ev: real = 1/ev2hartree;

// Mass units
const amu2me: real = 1822.882;       // Conversion from mass in amu to mass in au (m_e)
const me2amu: real = 1/amu2me;       // Conversion from mass in au (m_e) to mass in amu 


// Time units
const tau2ps: real = 41341.447;      // Conversion from time in au to time in ps
const ps2tau: real = 1/tau2ps;       // inverse

// Derived quantities
const Rgas: real = Kboltz*hartree2kcal*1000.0; // gas constant R = 1.98722 cal/mole/K

//
const TOOSMALL: real = 1.0E-09;
const MAX_TRIES: int = 30;
const OK: int = 0;
const BASIS_SIZE: int = 9;
const MATRIX_SIZE: int = BASIS_SIZE * BASIS_SIZE;
