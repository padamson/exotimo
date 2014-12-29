/*
import sys
from NumWrap import matrixmultiply,transpose
from math import exp,log
from Constants import Kboltz
from LA2 import mkdens
import logging
*/

/*
def mkdens_fermi(nel,orbe,orbs,e_temp):
    efermi = get_efermi(nel,orbe,e_temp)
    occs = get_fermi_occs(efermi,orbe,e_temp)
    D = mkdens_occs(orbs,occs)
    entropy = get_entropy(occs,e_temp)
    return D,entropy
*/

/*
def mkdens_occs(c,norbs:int,occs:norbs*real,tol:real=1.0E-5,verbose:bool=False):
    //Density matrix from a set of occupations (e.g. from FD expression).
    //Determine how many orbs have occupations greater than 0
    var noccs:int = 0;
    for i in [1..norbs] {
        if occs(i) >= tol then noccs += 1;
    }
    if verbose then {
        writeln("mkdens_occs: %d occupied orbitals found", noccs);
    }
    //Determine how many doubly occupied orbitals we have
    var nclosed:int = 0
    for i in [1..norbs] {
        if abs(1.0-occs(i)) <= tol then nclosed += 1;
    if verbose then {
        writeln("mkdens_occs: %d closed-shell orbitals found", nclosed);
    }
    //made it to here.
    D = mkdens(c,0,nclosed);
    for i in range(nclosed,norb):
        D = D + occs[i]*matrixmultiply(c[:,i:i+1],transpose(c[:,i:i+1]))
    return D
    
def get_fermi_occ(efermi,en,temp):
    kT = Kboltz*temp
    x = (en-efermi)/kT
    if x < -50.: return 1.
    elif x > 50.: return 0
    return 1/(1+exp(x))

def get_entropy(occs,temp):
    kT = Kboltz*temp
    entropy = 0
    for fi in occs:
        if abs(fi) < 1e-10: break # stop summing when occs get small
        if fi > 1e-10:
            entropy += kT*fi*log(fi)
        if (1-fi) > 1e-10:
            entropy += kT*(1.-fi)*log(1.-fi)
    return entropy
    

def get_fermi_occs(efermi,orbe,temp):
    occs = []
    for en in orbe:
        occs.append(get_fermi_occ(efermi,en,temp))
    return occs

def get_t0_occs(nel,nbf):
    occs = [0]*nbf
    nc,no = divmod(nel,2)
    for i in range(nc): occs[i] = 1.
    for i in range(nc,nc+no): occs[i] = 0.5
    return occs

def get_efermi(nel,orbe,temp,**opts):
    "Bisection method to get Fermi energy from Fermi-Dirac dist"
    tol = opts.get('tol',1e-9)
    verbose = opts.get('verbose',True)

    elow,ehigh = orbe[0]-100.,orbe[-1]
    nlow = 2*sum(get_fermi_occs(elow,orbe,temp))
    nhigh = 2*sum(get_fermi_occs(ehigh,orbe,temp))

    if nlow > nel:
        logging.error("elow incorrect %f -> %f " % (elow,nlow))
        raise Exception("elow incorrect %f -> %f " % (elow,nlow))
    if nhigh < nel:
        logging.error("ehigh incorrect %f -> %f " % (ehigh,nhigh))
        raise Exception("ehigh incorrect %f -> %f " % (ehigh,nhigh))

    for i in range(100):
        efermi = (elow+ehigh)/2
        n = 2*sum(get_fermi_occs(efermi,orbe,temp))
        if abs(n-nel) < tol:
            break
        elif n < nel:
            elow = efermi
        else:
            ehigh = efermi
    else:
        print "get_fd_occs: Too many iterations"
    return efermi
*/

