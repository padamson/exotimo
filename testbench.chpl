const MAX_ITERATIONS: int= 50;

const D: domain(1) = {1..400};
var HF, R, Cbar, Rold: [D] real;
var epsilon: [1..7] real;
var E, Rsum, term: real;
var pppone: real;
var m, n, i, icon, mm: int;

proc readAtom(filename:string) {
  var infile = open(filename, iomode.r);
  var reader = infile.reader();

  var m = reader.read(int),
      n = reader.read(int);

  const Dxy: domain(2) = {1..m,1..2};
  var xy: [Dxy] real;
  for i in 1..m {
    xy[i,1] = reader.read(real);
    xy[i,2] = reader.read(real);
  }
  writeln("m=",m);
  writeln("n=",n);
  for i in 1..m {
    writeln("xy[",i,",1]=",xy[i,1]);
    writeln("xy[",i,",2]=",xy[i,2]);
  }
}

readAtom("exotimo.in");
