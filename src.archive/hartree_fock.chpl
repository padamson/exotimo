use Time;
use Ints;
use LA2;

/*
from fermi_dirac import get_efermi, get_fermi_occs,mkdens_occs,get_entropy
from LA2 import geigh,mkdens,trace2
from Ints import get2JmK,getbasis,getints,getJ,getK
from Convergence import DIIS
import logging

from math import sqrt,pow
from PyQuante.cints import dist

from NumWrap import array2string, fabs, dot, transpose
*/


proc get_fock(D:[] real,Ints: [] real,h:[] real){
    //Return the HF Fock matrix
    return h + get2JmK(Ints,D);
}

proc get_energy(h:[] real,F:[] real,D:[] real,enuke:real=0.0):real{
    //Form the total energy of the closed-shell wave function.
    var eone:real = trace2(D,h);
    var etwo:real = trace2(D,F);
    writeln("eone,etwo,enuke",(2*eone),(etwo-eone),enuke);
    return eone+etwo+enuke;
}

proc get_guess(h:[] real,S:[] real){
  const D: domain(2) = [h.domain];
  //Form an initial guess from the one-electron Hamiltonian
  var evals:[D.dim(1)] real;
  var evecs:[h.domain] real;
  (evals,evecs) = geigh1(h,S);
  return evecs;
}

/*

proc get_nel(atoms,charge=0):
    print "Warning, hartree_fock.get_nel deprecated"
    print "Use the Molecular instance function"
    return atoms.get_nel(charge)

proc get_enuke(atoms):
    print "Warning, hartree_fock.get_enuke deprecated"
    print "Use the Molecular instance function"
    return atoms.get_enuke()

proc hf(atoms,**opts):
    """\
    Experimental inteface to rhf/uhf routines. See those
    routines for more information.
    """
    # Ultimately this will check spin_type to determine whether
    # to call uhf or rohf for open shell
    nclosed,nopen = atoms.get_closedopen()
    if nopen: return uhf(atoms,**opts)
    return rhf(atoms,**opts)

proc rhf(atoms,**opts):
    """\
    rhf(atoms,**opts) - Closed-shell HF driving routine
    
    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    ConvCriteria  1e-4    Convergence Criteria
    MaxIter       20      Maximum SCF iterations
    DoAveraging   True    Use DIIS for accelerated convergence (default)
                  False   No convergence acceleration
    ETemp         False   Use ETemp value for finite temperature DFT (default)
                  float   Use (float) for the electron temperature
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not none, the guess orbitals
    """
    ConvCriteria = opts.get('ConvCriteria',1e-4)
    MaxIter = opts.get('MaxIter',20)
    DoAveraging = opts.get('DoAveraging',False)
    ETemp = opts.get('ETemp',False)

    bfs = opts.get('bfs',None)
    if not bfs:
        basis_data = opts.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = opts.get('integrals', None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()

    orbs = opts.get('orbs',None)
    if orbs is None: orbe,orbs = geigh(h,S)

    nclosed,nopen = atoms.get_closedopen()
    nocc = nclosed
    assert(nopen == 0), "SCF currently only works for closed-shell systems"
    enuke = atoms.get_enuke()
    eold = 0.

    logging.info("RHF calculation on %s" % atoms.name)
    logging.info("Nbf = %d" % len(bfs))
    logging.info("Nclosed = %d" % nclosed)

    if DoAveraging:
        logging.info("Using DIIS averaging")
        avg = DIIS(S)
    logging.debug("Optimization of HF orbitals")
    for i in range(MaxIter):
        if ETemp:
            efermi = get_efermi(nel,orbe,ETemp)
            occs = get_fermi_occs(efermi,orbe,ETemp)
            D = mkdens_occs(orbs,occs)
            entropy = get_entropy(occs,ETemp)
        else:
            D = mkdens(orbs,0,nocc)
        G = get2JmK(Ints,D)
        F = h+G
        if DoAveraging: F = avg.getF(F,D)
        orbe,orbs = geigh(F,S)
        energy = get_energy(h,F,D,enuke)
        if ETemp:
            energy += entropy
        logging.debug("%d %f" % (i,energy))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final HF energy for system %s is %f" % (atoms.name,energy))
    #print 2*trace2(D,h),trace2(D,G)
    return energy,orbe,orbs

proc uhf(atoms,**opts):
    """\
    uhf(atoms,**opts) - Unrestriced Open Shell Hartree Fock
    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    ConvCriteria  1e-4    Convergence Criteria
    MaxIter       20      Maximum SCF iterations
    DoAveraging   True    Use DIIS averaging for convergence acceleration
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not None, the guess orbitals
    """
    ConvCriteria = opts.get('ConvCriteria',1e-5)
    MaxIter = opts.get('MaxIter',40)
    DoAveraging = opts.get('DoAveraging',True)
    averaging = opts.get('averaging',0.5)
    ETemp = opts.get('ETemp',False)

    bfs = opts.get('bfs',None)
    if not bfs:
        basis_data = opts.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = opts.get('integrals', None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()

    nalpha,nbeta = atoms.get_alphabeta() #pass in opts for multiplicity
    S,h,Ints = getints(bfs,atoms)

    orbs = opts.get('orbs',None)
    if orbs!=None:
        orbsa = orbsb = orbs
    else:
        orbe,orbs = geigh(h,S)
        orbea = orbeb = orbe
        orbsa = orbsb = orbs

    enuke = atoms.get_enuke()
    eold = 0.

    logging.info("UHF calculation on %s" % atoms.name)
    logging.info("Nbf = %d" % len(bfs))
    logging.info("Nalpha = %d" % nalpha)
    logging.info("Nbeta = %d" % nbeta)
    logging.info("Averaging = %s" % DoAveraging)
    logging.debug("Optimization of HF orbitals")
    for i in range(MaxIter):
        print "SCF Iteration:",i,"Starting Energy:",eold
        if ETemp:
            # We have to multiply nalpha and nbeta by 2
            #  to get the Fermi energies to come out correct:
            efermia = get_efermi(2.0*nalpha,orbea,ETemp)
            occsa = get_fermi_occs(efermia,orbea,ETemp)
            #print "occsa = ",occsa
            Da = mkdens_occs(orbsa,occsa)
            efermib = get_efermi(2.0*nbeta,orbeb,ETemp)
            occsb = get_fermi_occs(efermib,orbeb,ETemp)
            #print "occsb = ",occsb
            Db = mkdens_occs(orbsb,occsb)
            entropy = 0.5*(get_entropy(occsa,ETemp)+get_entropy(occsb,ETemp))
        else:
            Da = mkdens(orbsa,0,nalpha)
            Db = mkdens(orbsb,0,nbeta)
        if DoAveraging:
            if i: 
                Da = averaging*Da + (1-averaging)*Da0
                Db = averaging*Db + (1-averaging)*Db0
            Da0 = Da
            Db0 = Db

        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)
        Ka = getK(Ints,Da)
        Kb = getK(Ints,Db)
        Fa = h+Ja+Jb-Ka
        Fb = h+Ja+Jb-Kb
        orbea,orbsa = geigh(Fa,S)
        orbeb,orbsb = geigh(Fb,S)
        energya = get_energy(h,Fa,Da)
        energyb = get_energy(h,Fb,Db)
        energy = (energya+energyb)/2+enuke
        Dab = Da+Db
        Eone = trace2(Dab,h)
        Ej = 0.5*trace2(Dab,Ja+Jb)
        Ek = -0.5*(trace2(Da,Ka)+trace2(Db,Kb))
        if ETemp: energy += entropy
        logging.debug("%d %f %f %f %f" % (i,energy,Eone,Ej,Ek))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
    logging.info("Final UHF energy for system %s is %f" % (atoms.name,energy))
    return energy,(orbea,orbeb),(orbsa,orbsb)


#added by Hatem H Helal hhh23@cam.ac.uk
proc mk_auger_dens(c, occ):
    "Forms a density matrix from a coef matrix c and occupations in occ"
    #count how many states we were given
    nstates = occ.shape[0]
    D = 0.0
    for i in range(nstates):
        D += occ[i]*dot( c[:,i:i+1], transpose(c[:,i:i+1]))
    #pad_out(D)
    return D

proc uhf_fixed_occ(atoms,occa, occb,**opts):
    """\
    occa and occb represent the orbital occupation arrays for 
    the calculating spin orbitals with holes
    
    uhf(atoms,**opts) - Unrestriced Open Shell Hartree Fock
    atoms       A Molecule object containing the molecule

    Options:      Value   Description
    --------      -----   -----------
    ConvCriteria  1e-4    Convergence Criteria
    MaxIter       20      Maximum SCF iterations
    DoAveraging   True    Use DIIS averaging for convergence acceleration
    bfs           None    The basis functions to use. List of CGBF's
    basis_data    None    The basis data to use to construct bfs
    integrals     None    The one- and two-electron integrals to use
                          If not None, S,h,Ints
    orbs          None    If not None, the guess orbitals
    """

    from biorthogonal import biorthogonalize,pad_out

    ConvCriteria = opts.get('ConvCriteria',1e-5)
    MaxIter = opts.get('MaxIter',40)
    DoAveraging = opts.get('DoAveraging',True)
    averaging = opts.get('averaging',0.5)
    ETemp = opts.get('ETemp',False)
    
    bfs = opts.get('bfs',None)
    if not bfs:
        basis_data = opts.get('basis_data',None)
        bfs = getbasis(atoms,basis_data)

    integrals = opts.get('integrals', None)
    if integrals:
        S,h,Ints = integrals
    else:
        S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()

    nalpha,nbeta = atoms.get_alphabeta() #pass in opts for multiplicity
    S,h,Ints = getints(bfs,atoms)

    orbsa = opts.get('orbsa',None)
    orbsb = opts.get('orbsb',None)
    if (orbsa!=None and orbsb!=None):
        orbsa = orbsa
        orbsb = orbsb
    else:
        orbe,orbs = geigh(h,S)
        orbea = orbeb = orbe
        orbsa = orbsb = orbs
    
    #print "A Trial Orbital Energies:\n", orbea

    print "A Trial Orbitals:\n"
    pad_out(orbsa)

    #print "B Trial Orbital Energies:\n",orbeb

    print "B Trial Orbitals:\n"
    pad_out(orbsb)
    
    enuke = atoms.get_enuke()
    eold = 0.

    for i in range(MaxIter):
        print "SCF Iteration:",i,"Starting Energy:",eold
        #save the starting orbitals
        oldorbs_a=orbsa
        oldorbs_b=orbsb

        Da = mk_auger_dens(orbsa,occa)
        Db = mk_auger_dens(orbsb,occb)
        #Da_std = mkdens(orbsa,0,nalpha)
        #Db_std = mkdens(orbsb,0,nbeta)
        #pad_out(Da - Da_std ) #use to test mk_aug_dens with ground state occupations
        #pad_out(Db - Db_std )
        
        
        Ja = getJ(Ints,Da)
        Jb = getJ(Ints,Db)
        Ka = getK(Ints,Da)
        Kb = getK(Ints,Db)
        Fa = h+Ja+Jb-Ka
        Fb = h+Ja+Jb-Kb

        orbea,orbsa = geigh(Fa,S)
        orbeb,orbsb = geigh(Fb,S)
        
        #save the new orbitals
        neworbs_a=orbsa
        neworbs_b=orbsb
        
        #now we biorthogonalize the new orbitals to the old ones
        #to setup occupation arrays for the next scf cycle
        orbsa = biorthogonalize(neworbs_a,oldorbs_a,S,nalpha,occa)
        orbsb = biorthogonalize(neworbs_b,oldorbs_b,S,nbeta,occb)
        
        energya = get_energy(h,Fa,Da)
        energyb = get_energy(h,Fb,Db)
        energy = (energya+energyb)/2+enuke
        Dab = Da+Db
        Eone = trace2(Dab,h)
        Ej = 0.5*trace2(Dab,Ja+Jb)
        Ek = -0.5*(trace2(Da,Ka)+trace2(Db,Kb))
        
        #print "%d %f %f %f %f" % (i,energy,Eone,Ej,Ek)
        
        logging.debug("%d %f %f %f %f" % (i,energy,Eone,Ej,Ek))
        if abs(energy-eold) < ConvCriteria: break
        eold = energy
        if i==(MaxIter-1):
            print "Warning: Reached maximum number of SCF cycles may want to rerun calculation with more SCF cycles"
    logging.info("Final UHF energy for system %s is %f" % (atoms.name,energy))
    return energy,(orbea,orbeb),(orbsa,orbsb)


if __name__ == '__main__':
    from Molecule import Molecule
    h2 = Molecule('H2',atomlist=[(1,(1.,0,0)),(1,(-1.,0,0))])
    en = rhf(h2)
    print en

*/
