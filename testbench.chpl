const MAX_ITERATIONS: int= 50;

var HF[400], R[400]: real;
var Cbar[400], epsilon[7], Rold[400]: real;
var E, Rsum, term: real;
var pppone: real;
var m, n, i, icon, mm: int;


proc readAtom(filename) {
  var infile = open(filename, iomode.r);
  var reader = infile.reader();

  var m = reader.read(int),
      n = reader.read(int);

  var 
