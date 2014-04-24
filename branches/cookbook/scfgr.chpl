use constants;

extern proc USE( i, j, k, l, val);

proc scfGR(R:[] real, G[] real, m: int, nfile: int){
  var val: real;
  var irs, itu, iru, its: int;
  var i, j, k, l, mu, ij, kl, getint: int;
  var pointer: int;

  
