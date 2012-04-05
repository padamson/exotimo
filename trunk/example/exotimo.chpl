use CGBF;
use Time;
use basis;
use Ints;
use eigen;
use hartree_fock;
use io;
//use rho2gamma;

const r: real = 0.75;
const nuke1 = (0.0,0.0,r);
const nuke2 = (0.0,0.0,-r);
const sPowerList:[1..3] int = (0,0,0);
const pPowerList:[1..3,1..3] int = ((1,0,0),(0,1,0),(0,0,1));
const testMatrix:[1..3,1..3] real = ((1.0,2.0,3.0),(4.0,5.0,6.0),(7.0,8.0,9.0));
var transposeTestMatrix:[1..3,1..3] real;
const ZERO: real = 0.00000001;

proc main (){

  var Hnum_prims: int = 3;
  var Hpexps: [1..3] real = (3.425251,0.62391373,0.1688554);
  var Hpcoefs: [1..3] real = (0.15432897,0.53532814,0.44463454);
  //writeln("dot(Hpexps,Hpcoefs)=",dot(Hpexps,Hpcoefs));
  var i,j,k,l:int;

  var atoms=
    new classical(
        first_charge=1.0,
        first_pos=nuke1);
  atoms.add_nucleus(charge=3.0,pos=nuke2);
  writeln("atoms=",atoms);

  var electronBasis = 
    new basis(
        first_num_prims=Hnum_prims,
        first_pexps=Hpexps,
        first_pcoefs=Hpcoefs,
        particle="electron",
        first_origin=nuke1,
        first_powers=(0,0,0),
        first_atid=1);
  
  electronBasis.readBasis(basisFileName="../Basis/sto3g.dat",atom="LI",basis="STO3G",
    origin=nuke2,atid=2);

  writeln("testMatrix=",testMatrix);
  transposeTestMatrix = transpose(testMatrix);
  writeln("transposeTestMatrix=",transposeTestMatrix);
  writeln("transpose(testMatrix)=",transpose(testMatrix));
  writeln("matrixmultiply(testMatrix,transposeTestMatrix)=",matrixmultiply(testMatrix,transposeTestMatrix));


  //var testMatrixTranspose:[1..3,1..3] int = transpose(testMatrix);
  //writeln("testMatrixTranspose(1,1)=", testMatrixTranspose(1,1));

  //writeln("reshape(pPowerList)=",reshape(pPowerList,[1..9]));
  
  var nbf = electronBasis.num_cgbfs;
  var num_nuclei = atoms.num_nuclei;
  var S,h,T,V:[1..nbf,1..nbf] real = 0.0;
  (S,h) = getints(electronBasis.cgbfs,atoms.nuclei);
  writeln("S=");
  write_matrix(S);
  writeln("h=");
  write_matrix(h);
  //writeln("S=",S);
  //writeln("h=",h);

  //geepp(electronBasis.cgbfs(3),electronBasis.cgbfs(4));
  //geepp(electronBasis.cgbfs(3),electronBasis.cgbfs(5));
  //geepp(electronBasis.cgbfs(3),electronBasis.cgbfs(6));

  /*
  var vecs: [S.domain] real;
  vecs = get_guess(h,S);
  writeln("vecs = ", vecs);
  */

  /*
  const D: domain(2) = [S.domain];
  var d,e:[D.dim(1)] real;
  var Z:[S.domain] real;
  (d,e,Z) = tred2(S);
  writeln("d = ", d);
  writeln("e = ", e);
  writeln("Z = ", Z);
  Z = cholesky(S);
  writeln("Z = ", Z);
  */

  /*
  T = getT(electronBasis.cgbfs);
  writeln("T=",T);
  S = getS(electronBasis.cgbfs);
  writeln("S=",S);
  V = getV(electronBasis.cgbfs,atoms.nuclei);
  writeln("V=",V);
  var totlen:int = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8;
  var Ints:[1..totlen] real = 0.0;
  var time1: real = getCurrentTime();
  Ints = get2ints(electronBasis.cgbfs);
  //writeln("Ints=",Ints);
  writeln("Took ",getCurrentTime()-time1," seconds to compute very slow coulomb integral code ",totlen," times.");
  */

  /*
  var J:[1..nbf,1..nbf] real = 0.0;
  J = getJ(nbf,totlen,D);
  */

  /*
  forall (i) in [1..6] {
     writeln("cgbfs=",i);
     writeln(electronBasis.cgbfs(i));
  }
  */
  
  //forall (i,j,k,l) in [1..6,1..1,1..1,1..1] {
  
  /*
  for i in [1..6] {
    for j in [1..i] {
      for k in [1..6] {
        for l in [1..k] {
     ints(i,j,k,l) = CGBF_coulomb(
           electronBasis.cgbfs(i),
           electronBasis.cgbfs(j),
           electronBasis.cgbfs(k),
           electronBasis.cgbfs(l));
        }
      }
    }
  }
  */
  
  
  /*
  writeln("int 1,1,1,1");
     ints(1,1,1,1) = CGBF_coulomb(
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1));
  writeln("int 1,1,1,2");
     ints(1,1,1,2) = CGBF_coulomb(
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(2));
  writeln("int 1,1,1,4");
     ints(1,1,1,4) = CGBF_coulomb(
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(4));
  writeln("int 1,1,1,5");
     ints(1,1,1,5) = CGBF_coulomb(
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(5));
  writeln("int 1,1,1,6");
     ints(1,1,1,6) = CGBF_coulomb(
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(1),
           electronBasis.cgbfs(6));
     */

  //writeln("Took ",getCurrentTime()-time1," seconds to compute very slow coulomb integral code ",6*6*6*6," times.");

  //forall (i,j,k,l) in [1..6,1..1,1..1,1..1] {
  
  /*
  for i in [1..6] {
    for j in [1..i] {
      for k in [1..6] {
        for l in [1..k] {
          if ints(i,j,k,l) > ZERO {
            writeln(i," ",j," ",k," ",l," ",ints(i,j,k,l));
          }
        }
      }
    }
  }
  */
  
  
  /*
      writeln("ints(",1,",",1,",",1,",",1,")= ", 
        ints(1,1,1,1));
      writeln("ints(",1,",",1,",",1,",",2,")= ", 
        ints(1,1,1,2));
      writeln("ints(",1,",",1,",",1,",",4,")= ", 
        ints(1,1,1,4));
      writeln("ints(",1,",",1,",",1,",",5,")= ", 
        ints(1,1,1,5));
      writeln("ints(",1,",",1,",",1,",",6,")= ", 
        ints(1,1,1,6));
        */

  /*
  var nbf: int = electronBasis.num_cgbfs;
  writeln("nbf = ",nbf);
  var totlen: int = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8;
  writeln("totlen = ",totlen);
  var ints2: [1..totlen] real = 0.0;
  var ijkl,ij,kl:int;
  var m:int = 0;

  var time1: real = getCurrentTime();
  for i in [0..nbf-1] {
    for j in [0..i] {
      ij = i*(i+1)/2+j;
      for k in [0..nbf-1] {
        for l in [0..k] {
          kl = k*(k+1)/2+l;
          if ij >= kl {
            m+=1;
            ijkl = ijkl2intindex(i,j,k,l);
            writeln("computing ", i+1, " ",j+1, " ",k+1, " ",l+1);
            ints2(ijkl) = CGBF_coulomb(
                  electronBasis.cgbfs(i+1),
                  electronBasis.cgbfs(j+1),
                  electronBasis.cgbfs(k+1),
                  electronBasis.cgbfs(l+1));
          }
        }
      }
    }
  }
  writeln("Took ",getCurrentTime()-time1," seconds to compute very slow coulomb integral code ",totlen," times.");
  writeln("Took ",getCurrentTime()-time1," seconds to compute very slow coulomb integral code ",m," times.");

  for i in [0..nbf-1] {
    for j in [0..i] {
      ij = i*(i+1)/2+j;
      for k in [0..nbf-1] {
        for l in [0..k] {
          kl = k*(k+1)/2+l;
          if ij >= kl {
            ijkl = ijkl2intindex(i,j,k,l);
            writeln(i+1,"   ",j+1,"   ",k+1,"   ",l+1,"   ",ints2(ijkl));
          }
        }
      }
    }
  }

  var kints:[1..6,1..6] real;
  var cgbf1,cgbf2:CGBF;
  writeln("========");

  forall (i,j) in [1..6,1..6] {
    cgbf1=electronBasis.cgbfs(i);
    cgbf2=electronBasis.cgbfs(j);
    writeln("kints(",i,",",j,")=",cgbf1.contr_kinetic(cgbf2));
  }

  //writeln("electronBasis=",electronBasis);

  writeln("px px s s coulomb_repulsion=",coulomb_repulsion((0.0, 0.0, -0.75),0.810044,(1, 0, 0),0.63629,(0.0, 0.0, -0.75),0.810044,(1, 0, 0),0.63629,(0.0, 0.0, 0.75),1.79444,(0, 0, 0),3.42525,(0.0, 0.0, 0.75),1.79444,(0, 0, 0),3.42525));
  writeln("py py s s coulomb_repulsion=",coulomb_repulsion((0.0, 0.0, -0.75),0.810044,(0, 1, 0),0.63629,(0.0, 0.0, -0.75),0.810044,(0, 1, 0),0.63629,(0.0, 0.0, 0.75),1.79444,(0, 0, 0),3.42525,(0.0, 0.0, 0.75),1.79444,(0, 0, 0),3.42525));
  writeln("pz pz s s coulomb_repulsion=",coulomb_repulsion((0.0, 0.0, -0.75),0.810044,(0, 0, 1),0.63629,(0.0, 0.0, -0.75),0.810044,(0, 0, 1),0.63629,(0.0, 0.0, 0.75),1.79444,(0, 0, 0),3.42525,(0.0, 0.0, 0.75),1.79444,(0, 0, 0),3.42525));
  writeln("pz pz pz pz coulomb_repulsion=",coulomb_repulsion((0.0, 0.0, -0.75),0.810044,(0, 0, 1),0.63629,(0.0, 0.0, -0.75),0.810044,(0, 0, 1),0.63629,(0.0, 0.0, -0.75),1.79444,(0, 0, 1),3.42525,(0.0, 0.0, -0.75),1.79444,(0, 0, 1),3.42525));

//writeln(electronBasis.particle," = ",electronBasis);

for i in 1..3 {
electronBasis.add_CGBF(
    origin=(0.0,0.0,-r),powers=pPowerList(i,1..3),atid=2,
    num_prims=3,
    pexps=(0.63628969999999996,0.14786009999999999,0.048088699999999998),
    pcoefs=(-0.099967230000000004, 0.39951282999999999, 0.70011546999999996));

electronBasis.add_CGBF(
    origin=(0.0,0.0,-r),powers=pPowerList(i,1..3),atid=2,
    num_prims=3,
    pexps=(0.63628969999999996,0.14786009999999999, 0.048088699999999998),
    pcoefs= (0.15591627, 0.60768372000000004, 0.39195739000000002));
}
*/

}
