use constants;

//getint for the special case of PPP integrals
proc getint(i:int, j:int, k:int, l:int, mu:int, val:real, pointer: int): int {
  var status, ppptwo: int;
  var id, jd, kd, ld: int;
  if (pointer == 0) { //initialise labels
    id = 0;
    kd = 0;
  }
  pointer++;
  kd++; // increment labels
  if ( kd > id ) {
    kd = 1;
    id++;
  }
  jd = id;
  ld = kd; // ZDO approx for integrals

  status = ppptwo(id, kd, val); // PPP repulsion integral generator
  // regenerate labels
  i = id;
  j = jd;
  k = kd;
  l = ld;

  return status;
}

// PPP ZDO integrals, values from Tom Peacock's book p97
proc ppptwo(i:int, j:int, val:real, xy[]: real, m:int):int {
  var r: real;

  if ( i > m ) return 0;

  if ( i == j ) {
    val = 11.4e0; // one-center value in eV
    return 1;
  }

  r = sqrt((xy(i,1)-xy(j,1))**2 + (xy(i,2) - xy(j,2))**2);
  if ( abs(r-TWO*RCC) < 0.1e0) { // next nearest-neighbours
    val = 10.528e0 - 2.625e0*r + 0.2157*r*r;
  } else { // charged-sphere formula
    val = (1.0e0 + 2.0e0/r)**(-0.5e0);
    val = 7.2e0*(1.0e0 + val)/r;
  } // values in eV
  return 1;
}

proc pppone(xy[]:real, m:int, i:int, j:int):real{
  var junk:int;
  var r, val, gamma: real;
  
  val = ZERO;
  if (i == j) {//diagonal values
    val = -11.16e0;
    for k in 1..m{
      if (i == k) continue;
      junk = ppptwo(i,k,gamma);
      val = val - gamma; //Nuclear attraction modelled by e- repulsion
    }
  } else { // off-diagonal integrals
    r = sqrt((xy[i,1]-xy[j,1])**2 + (xy[i,2]-xy[j,2])**2);
    if (abs(r-1.414e0) < 0.1e0){ // C-C distance
      val = -2.395e0; //Beta in eV
    }
  }
  return val;
}


