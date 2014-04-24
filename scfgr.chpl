use constants;

extern proc USE( i, j, k, l, val);

proc scfGR(R:[] real, G[] real, m: int){
  var val: real;
  var irs, itu, iru, its: int;
  var i, j, k, l, mu, ij, kl, getint: int;
  var pointer: int;

  pointer = 0; // initialise 'file'
  while ( getint(i,j,k,l,mu,val,pointer) != 0 ) {
    USE(i,j,k,l,val); // 1
    if (i != j) then USE(j,i,k,l,val); // 2
    if (k != l) {
      USE(i,j,l,k,val); // 3
      if (i != j ) then USE(j,i,l,k,val); // 4
    }
    ij = (i*(i-1))/2 + j;
    kl = (k*(k-1))/2 + l;
    if (ij != kl ) {
      USE(k,l,i,j,val); // 5
      if (i != j) then USE(k,l,j,i,val); // 6
      if (k != l) {
        USE(l,k,i,j,val); // 7
        if (i != j) then USE(l,k,j,i,val); // 8
      }
    }
  } // end of while
  return;
}
