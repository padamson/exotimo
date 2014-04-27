use io;
use constants;
use getint;

var m, n:int;
var Dmat: domain(1) = {1..2};
var Dxy: domain(2) = {1..2,1..2}; 
var HF, R, Cbar, Rold: [Dmat] real;
var xy:[Dxy] real;

proc main () {
  var epsilon: [1..7] real;
  var E, Rsum, term: real;
  var val, pppone: real;
  var kount, ij, ji, icon, mm: int;

  readAtom("/Users/padamson/Research/exotimo/branches/cookbook/data/exotimo.in");
  
  mm = m*m;
  //Set initial R to ZERO, i.e. start from H not HF
  R = ZERO; 
  Rold = ZERO;
  //Assume that the first iteration will change R to be non-ZERO!
  kount = 0;
  do {
    kount += 1;
    E = ZERO; //Initialize E 
    icon = 0; //Initialize counter
    //Compute the one-electron Hamiltonian each time
    for i in 1..m {
      for j in 1..i {
        ij = m*(j-1) + i;
        ji = m*(i-1) + j;
        val = pppone(xy,m,i,j);
        HF(ij) = val;
        HF(ji) = val;
      }
    }

  } while (icon != 0 && kount < MAX_ITERATIONS);

}

proc readAtom(filename:string) {
  var infile = open(filename, iomode.r);
  var reader = infile.reader();

  m = reader.read(int);
  n = reader.read(int);

  Dxy = {1..m,1..2};

  read_matrix(xy, reader);

  writeln("m=",m);
  writeln("n=",n);
  write_matrix(xy);
}

