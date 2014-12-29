proc SIGN(a,b){
  if (b >= 0.0) then {
    return abs(a);
  } else {
    return -abs(a);
  }
}

proc pythag(a:real,b:real):real{
  return sqrt(a**2 + b**2);
}

proc tred2(A:[] real){
  const D: domain(2) = [A.domain];
  const n: int = D.high(1);
  var d,e:[D.dim(1)] real;
  var a:[A.domain] real = A;
  var scale, hh, h, g, f: real;
  var l: int;

  for i in [n..2 by -1] {
    l=i-1;
    h=0.0;
    scale=0.0;
    if (l > 1) then {
      for k in [1..l] {
        scale += abs(a(i,k));
      }
      if (scale == 0.0) then {
        e(i) = a(i,l);
      } else {
        for k in [1..l] {
          a(i,k) /= scale;
          h += a(i,k) * a(i,k);
        }
        f = a(i,l);
        if (f >= 0.0) then {
          g = -sqrt(h);
        } else {
          g = sqrt(h);
        }
        e(i) = scale * g;
        h -= f*g;
        a(i,l) = f-g;
        f=0.0;
        for j in [1..l] {
          a(j,i) = a(i,j)/h;
          g=0.0;
          for k in [1..j] {
            g += a(j,k) * a(i,k);
          }
          for k in [j+1..l] {
            g += a(k,j)*a(i,k);
          }
          e(j) = g/h;
          f += e(j)*a(i,j);
        }
        hh = f/(h+h);
        for j in [1..l] {
          f = a(i,j);
          g = e(j)-hh*f;
          e(j) = g; 
          for k in [1..j] {
            a(j,k) -= (f*e(k)+g*a(i,k));
          }
        }
      }
    } else {
      e(i) = a(i,l);
    }
    d(i) = h;
  }
  d(1) = 0.0;
  e(1) = 0.0;
  for i in [1..n] {
    l=i-1;
    if (d(i) > 0.0) then {
      for j in [1..l] {
        g=0.0;
        for k in [1..l] {
          g += a(i,k)*a(k,j);
        }
        for k in [1..l] {
          a(k,j) -= g*a(k,i);
        }
      }
    }
    d(i) = a(i,i);
    a(i,i)=1.0;
    for j in [1..l] {
      a(j,i) = 0.0;
      a(i,j) = 0.0;
    }
  }
  return (d,e,a);
}

proc tqli(D:[] real, E:[] real, Z:[] real){
  const Dom: domain(1) = [D.domain];
  const n: int = Dom.high;
  var d: [D.domain] real = D;
  var e: [E.domain] real = E;
  var z: [Z.domain] real = Z;
  var ii, m, iiter: int;
  var s,r,p,g,f,dd,c,b: real;

  for i in [2..n] {
    e(i-1) = e(i);
  }
  e(n) = 0.0;
  for l in [1..n] {
    iiter = 0;
    do {
      for mm in [l..n-1 by 1] {
        dd = abs(d(mm)) + abs(d(mm+1));
        if (abs(e(mm)+dd == dd)) then {break;}
        m=mm;
        writeln("m=",m);
      }
      if (m != l) then {
        iiter += 1;
        if (iiter == 30) then {
          writeln("Too many iterations in tqli.");
          halt();
        }
        g = (d(l+1)-d(l))/(2.0*e(l));
        r = pythag(g,1.0);
        g=d(m)-d(l)+e(l)/(g+SIGN(r,g));
        c = 1.0;
        s = c;
        p = 0.0;
        for i in [m-1..l by -1] {
          f = s*e(i);
          b = c*e(i);
          r = pythag(f,g);
          e(i+1) = r;
          if (r == 0.0) then {
            d(i+1) -= p;
            e(m) = 0.0;
            break;
          }
          s = f/r;
          c = g/r;
          g = d(i+1) - p;
          r = (d(i)-g)*s+2.0*c*b;
          p = s*r;
          d(i+1) = g+p;
          g = c*r-b;
          for k in [1..n] {
            f=z(k,i+1);
            z(k,i+1) = s*z(k,i)+c*f;
            z(k,i) = c*z(k,i)-s*f;
          }
          ii=i;
        }
        if (r == 0.0 && ii >= 1) then continue;
        d(l) -= p;
        e(l) = g;
        e(m) = 0.0;
      }
    } while (m != l);
  }
  return (d,z);
}

proc eigh(S:[] real){
  const D: domain(2) = [S.domain];
  var d,e:[D.dim(1)] real;
  var Z:[S.domain] real;
  (d,e,Z) = tred2(S);
  return tqli(d,e,Z);
}
  



