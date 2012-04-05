/*******************************************************************************
/  miscellaneous constants
*******************************************************************************/

use constants;
const ITMAX: int = 100;
const EPS: real = 3.0e-7;
const FPMIN: real = 1.0e-30;
const SMALL: real = 1.0e-8;

/*******************************************************************************
/  miscellaneous functions
*******************************************************************************/

/*
proc lgamma2(z:real):real{
  if (z<=0) then return 0.0;

  var c: [0..6] real = (
    2.5066282746310005,
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5);

  var x:real = z;
  var y:real = x;
  var tmp: real = x + 5.5;
  tmp = (x+0.5)*log(tmp)-tmp;
  var ser: real = 1.000000000190015;

  for i in [1..6] {
    y += 1.0;
    ser += c[i]/y;
  }

  return tmp+log(c[0]*ser/x);

}
*/

proc fB(i: int, l1: int, l2: int, px: real, ax: real, bx: real, 
		 r: int, g: real): real {
  return binomial_prefactor(i,l1,l2,px-ax,px-bx)*Bfunc(i,r,g);
}

proc Bfunc(i: int, r: int, g: real): real{
  return fact_ratio2(i,r)*(4*g)**(r-i);
}

proc binomial_prefactor(s: int, ia: int, ib: int, xpa: real, xpb: real): real{
  var t: int;
  var sum: real=0.0;
  for t in [0..s]{
    if ((s-ia <= t) && (t <= ib)) then
      sum += binomial(ia,s-t)*binomial(ib,t)* xpa**(ia-s+t) * xpb**(ib-t);
  }
  return sum;
} 

proc binomial(a: int, b: int): int {return fact(a)/(fact(b)*fact(a-b));}

proc fact_ratio2(a: int, b: int): int { return fact(a)/fact(b)/fact(a-2*b); }

proc fact(n: int): int{
  if (n <= 1) then return 1;
  return n*fact(n-1);
}

/* double factorial function = 1*3*5*...*n */
proc fact2(n: int): int{ 
  if (n <= 1) then return 1;
  return n*fact2(n-2);
}

proc dist2(x1: 3*real, x2: 3*real): real{
  return (x1(1)-x2(1))*(x1(1)-x2(1))
    +(x1(2)-x2(2))*(x1(2)-x2(2))
    +(x1(3)-x2(3))*(x1(3)-x2(3));
}

proc dist(x1: 3*real, x2: 3*real): real{
  return sqrt(dist2(x1,x2));
}

proc product_center_1D(alphaa: real, xa: real, alphab: real, xb: real): real{
  return (alphaa*xa+alphab*xb)/(alphaa+alphab);
}

proc Fgamma(m: real, x0: real): real{
  var val: real;
  var x: real = x0;
  if ( abs(x) < SMALL ) then x = SMALL;
  val = gamm_inc(m+0.5,x);
  /* if (val < SMALL) return 0.; */ /* Gives a bug for D orbitals. */
  return 0.5* x**(-m-0.5) * val; 
}

proc gamm_inc(a: real, x: real): real{ /* Taken from NR routine gammap */
  var gamser,gammcf,gln: real;
  
  assert (x >= 0.0);
  assert (a > 0.0);
  if (x < (a+1.0)) then {
    (gamser,gln)=gser(a,x);
    return exp(gln)*gamser;
  } else {
    (gammcf,gln)=gcf(a,x);
    return exp(gln)*(1.0-gammcf);
  }
}
 
proc gser(a: real, x: real): [1..2] real {
  var n: int;
  var sum,del,ap: real;
  const gln: real=lgamma(a);
  var gamser: real;

  if (x <= 0.0) then {
    assert(x>=0.0);
    gamser=0.0;
    return (gamser,gln);
  } else {
    ap=a;
    del=1.0/a;
    sum=del;
    for n in [1..ITMAX] {
      ap+=1;
      del *= x/ap;
      sum += del;
      if (abs(del) < abs(sum)*EPS) then {
	gamser=sum*exp(-x+a*log(x)-(gln));
	return (gamser,gln);
      }
    }
    writeln("a too large, ITMAX too small in routine gser");
    return (0.0,0.0);
  }
}
 
proc gcf(a: real, x: real):[1..2] real {
  var i: int;
  var an,b,c,d,del,h: real;
  
  const gln: real = lgamma(a);
  var gammcf: real;

  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for i in [1..ITMAX] {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (abs(d) < FPMIN) then d=FPMIN;
    c=b+an/c;
    if (abs(c) < FPMIN) then c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (abs(del-1.0) < EPS) then break;
  }
  assert(i<=ITMAX);
  gammcf=exp(-x+a*log(x)-(gln))*h;
  return (gammcf,gln);
}

proc ijkl2intindex(i0: int, j0: int, k0: int, l0: int): int{
  var tmp,ij,kl: int;
  var i:int=i0;
  var j:int=j0;
  var k:int=k0;
  var l:int=l0;
  if (i<j) {
    tmp = i;
    i = j;
    j = tmp;
  }
  if (k<l) {
    tmp = k;
    k = l;
    l = tmp;
  }
  ij = i*(i+1)/2+j;
  kl = k*(k+1)/2+l;
  if (ij<kl) then {
    tmp = ij;
    ij = kl;
    kl = tmp;
  }
  return ij*(ij+1)/2+kl;
}

proc A_term(i:int, r:int, u:int, l1:int, l2:int, 
    PAx:real, PBx:real, CPx:real, gamma:real): real {
  /* THO eq. 2.18 */
  return -1**i * binomial_prefactor(i,l1,l2,PAx,PBx)*
    -1**u * fact(i)* CPx**(i-2*r-2*u) *
    (0.25/gamma)**(r+u)/fact(r)/fact(u)/fact(i-2*r-2*u);
}

proc A_array(l1:int, l2:int, PA:real, PB:real,
		CP:real, g:real): [0..l1+l2] real{
  /* THO eq. 2.18 and 3.1 */
  var i,r,u,I: int;
  const Imax: int = l1+l2;
  var A: [0..Imax] real;

  forall (i,r,u) in [0..Imax, 0..floor(i/2.0):int, 0..floor((i-2*r)/2.0):int] {
	I = i-2*r-u;
	A[I] += A_term(i,r,u,l1,l2,PA,PB,CP,g);
  }
  return A;
}

proc B_array(l1: int, l2: int, l3: int, l4: int, p: real, a: real,
		b: real, q: real, c: real, d: real,
		g1: real, g2: real, delta: real):[0..l1+l2+l3+l4] real{
  var I,i1,i2,r1,r2,u: int;
  var B: [0..l1+l2+l3+l4] real = 0.0;

  forall (i1,i2,r1,r2,u) in [0..l1+l2, 0..l3+l4, 0..i1/2, 0..i2/2, 0..(i1+i2)/2-r1-r2] {
    I = i1+i2-2*(r1+r2)-u;
    B[I] = B[I] + B_term(i1,i2,r1,r2,u,l1,l2,l3,l4,p,a,b,q,c,d,g1,g2,delta);
  }
  return B;
}

proc B_term(i1: int, i2: int, r1: int, r2: int, u: int, l1: int, l2: int,
	      l3: int, l4: int, Px: real, Ax: real, Bx: real,
	      Qx: real, Cx: real, Dx: real, gamma1: real,
	      gamma2: real, delta: real): real{
  /* THO eq. 2.22 */
  return fB(i1,l1,l2,Px,Ax,Bx,r1,gamma1)
    *((-1)**i2)*fB(i2,l3,l4,Qx,Cx,Dx,r2,gamma2)
    *((-1)**u)*fact_ratio2(i1+i2-2*(r1+r2),u)
    *(Qx-Px)**(i1+i2-2*(r1+r2)-2*u) / delta**(i1+i2-2*(r1+r2)-u);
}

proc coulomb_repulsion(
    xa: 3*real, norma: real, la: 3*int, alphaa: real,
    xb: 3*real, normb: real, lb: 3*int, alphab: real,
    xc: 3*real, normc: real, lc: 3*int, alphac: real,
    xd: 3*real, normd: real, ld: 3*int, alphad: real): real{

  var xp,xq: 3*real;
  var rab2,rcd2,rpq2,gamma1,gamma2,delta,sum: real;
  var Bx: [0..la(1)+lb(1)+lc(1)+ld(1)] real;
  var By: [0..la(2)+lb(2)+lc(2)+ld(2)] real;
  var Bz: [0..la(3)+lb(3)+lc(3)+ld(3)] real;
  var I,J,K: int;

  rab2 = dist2(xa,xb);
  rcd2 = dist2(xc,xd);
  xp(1) = product_center_1D(alphaa,xa(1),alphab,xb(1));
  xp(2) = product_center_1D(alphaa,xa(2),alphab,xb(2));
  xp(3) = product_center_1D(alphaa,xa(3),alphab,xb(3));
  xq(1) = product_center_1D(alphac,xc(1),alphad,xd(1));
  xq(2) = product_center_1D(alphac,xc(2),alphad,xd(2));
  xq(3) = product_center_1D(alphac,xc(3),alphad,xd(3));
  rpq2 = dist2(xp,xq);
  gamma1 = alphaa+alphab;
  gamma2 = alphac+alphad;
  delta = (1.0/gamma1+1.0/gamma2)/4.0;
  //writeln("rab2,rcd2,xp,xq,rpq2,gamma1,gamma2"," ",rab2," ",rcd2," ",xp," ",xq," ",rpq2," ",gamma1," ",gamma2);

  Bx = B_array(la(1),lb(1),lc(1),ld(1),xp(1),xa(1),xb(1),xq(1),xc(1),xd(1),gamma1,gamma2,delta);
  By = B_array(la(2),lb(2),lc(2),ld(2),xp(2),xa(2),xb(2),xq(2),xc(2),xd(2),gamma1,gamma2,delta);
  Bz = B_array(la(3),lb(3),lc(3),ld(3),xp(3),xa(3),xb(3),xq(3),xc(3),xd(3),gamma1,gamma2,delta);
  //writeln("Bx ",Bx);
  //writeln("By ",By);
  //writeln("Bz ",Bz);

  sum = 0.0;
  forall (I,J,K) in [0..la(1)+lb(1)+lc(1)+ld(1), 0..la(2)+lb(2)+lc(2)+ld(2), 0..la(3)+lb(3)+lc(3)+ld(3)] {
    sum += Bx[I]*By[J]*Bz[K]*Fgamma(I+J+K,0.25*rpq2/delta);
  }

  return 2.0*(PI**2.5)/(gamma1*gamma2*sqrt(gamma1+gamma2))
    *exp(-alphaa*alphab*rab2/gamma1) 
    *exp(-alphac*alphad*rcd2/gamma2)
    *sum*norma*normb*normc*normd;
}

proc contr_coulomb(
    lena: int, aexps: [1..lena] real, acoefs: [1..lena] real, anorms: [1..lena] real, xa: 3*real, la: 3*int, 
    lenb: int, bexps: [1..lenb] real, bcoefs: [1..lenb] real, bnorms: [1..lenb] real, xb: 3*real, lb: 3*int, 
    lenc: int, cexps: [1..lenc] real, ccoefs: [1..lenc] real, cnorms: [1..lenc] real, xc: 3*real, lc: 3*int, 
    lend: int, dexps: [1..lend] real, dcoefs: [1..lend] real, dnorms: [1..lend] real, xd: 3*real, ld: 3*int): real{

  var Jij: real = 0.0;

  forall (i,j,k,l) in [1..lena,1..lenb,1..lenc,1..lend]{
    Jij += acoefs(i)*bcoefs(j)*ccoefs(k)*dcoefs(l)
      * coulomb_repulsion(
          xa,anorms(i),la,aexps(i), xb,bnorms(j),lb,bexps(j), 
          xc,cnorms(k),lc,cexps(k), xd,dnorms(l),ld,dexps(l));
    //writeln("i j k l Jij ", i, " ",j, " ",k, " ",l, " ",Jij);
    //writeln("coefs ", " ",acoefs(i)," ",bcoefs(j)," ",ccoefs(k)," ",dcoefs(l));
    //writeln("xa anorms la aexps ",xa," ",anorms(i)," ",la," ",aexps(i));
    //writeln("xb bnorms lb bexps ",xb," ",bnorms(j)," ",lb," ",bexps(j));
    //writeln("xc cnorms lc cexps ",xc," ",cnorms(k)," ",lc," ",cexps(k));
    //writeln("xd dnorms ld dexps ",xd," ",dnorms(l)," ",ld," ",dexps(l));
  }
  return Jij;
}

proc nuclear_attraction(
    alpha1:real, l1:3*int, x1:3*real,
    alpha2:real, l2:3*int, x2:3*real, x3:3*real): real{
  var I,J,K:int;
  var xp:3*real;
  var gamma,sum,rab2,rcp2:real;
  var Ax:[0..l1(1)+l2(1)] real,Ay:[0..l1(2)+l2(2)] real,Az:[0..l1(3)+l2(3)] real;

  gamma = alpha1+alpha2;

  xp(1) = product_center_1D(alpha1,x1(1),alpha2,x2(1));
  xp(2) = product_center_1D(alpha1,x1(2),alpha2,x2(2));
  xp(3) = product_center_1D(alpha1,x1(3),alpha2,x2(3));

  rab2 = dist2(x1,x2);
  rcp2 = dist2(x3,xp);

  Ax = A_array(l1(1),l2(1),xp(1)-x1(1),xp(1)-x2(1),xp(1)-x3(1),gamma);
  Ay = A_array(l1(2),l2(2),xp(2)-x1(2),xp(2)-x2(2),xp(2)-x3(2),gamma);
  Az = A_array(l1(3),l2(3),xp(3)-x1(3),xp(3)-x2(3),xp(3)-x3(3),gamma);

  sum = 0.0;
  forall (I,J,K) in [0..l1(1)+l2(1), 0..l1(2)+l2(2), 0..l1(3)+l2(3)]{
    sum += Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K,rcp2*gamma);
  }

  return -2*PI/gamma*exp(-alpha1*alpha2*rab2/gamma)*sum;
}

proc three_center_1D(xi: real, ai: int, alphai: real,
      xj: real, aj: int, alphaj: real,
      xk: real, ak: int, alphak: real): real {

  var gamma, dx, px, xpi,xpj,xpk,intgl: real;
  var q,r,s,n: int;
  
  gamma = alphai+alphaj+alphak;
  dx = exp(-alphai*alphaj*(xi-xj)**2/gamma) *
    exp(-alphai*alphak*(xi-xk)**2/gamma) *
    exp(-alphaj*alphak*(xj-xk)**2/gamma);
  px = (alphai*xi+alphaj*xj+alphak*xk)/gamma;
    
  xpi = px-xi;
  xpj = px-xj;
  xpk = px-xk;
  intgl = 0.0;
  forall (q,r,s) in [0..ai, 0..aj, 0..ak] {
    if ((q+r+s)%2 == 0) then {
      n = (q+r+s)/2;
      intgl += binomial(ai,q)*binomial(aj,r)*binomial(ak,s)*
	    xpi**(ai-q) * xpj**(aj-r) * xpk**(ak-s) *
	    fact2(2*n-1)/(2*gamma)**n*sqrt(PI/gamma);
    }
  }
  return dx*intgl;
}

proc overlap(alpha1:real, l1:3*int, xa:3*real, 
            alpha2:real, l2:3*int, xb:3*real): real{
  /*Taken from THO eq. 2.12*/
  var xp,wx: 3*real;
  var rab2,gamma,pre: real;

  rab2 = dist2(xa,xb);
  gamma = alpha1+alpha2;

  for i in 1..3 do xp(i) = product_center_1D(alpha1,xa(i),alpha2,xb(i));

  pre = (PI/gamma)**1.5 * exp(-alpha1*alpha2*rab2/gamma);

  for i in 1..3 do wx(i) = overlap_1D(l1(i),l2(i),xp(i)-xa(i),xp(i)-xb(i),gamma);

  return pre*wx(1)*wx(2)*wx(3);
}

proc overlap_1D(l1:int, l2:int, PAx:real,PBx:real, gamma:real): real{
  /*Taken from THO eq. 2.12*/
  var i:int;
  var sum:real;
  for i in [0..floor(0.5*(l1+l2)):int] {
    sum += binomial_prefactor(2*i,l1,l2,PAx,PBx)* 
      fact2(2*i-1)/((2*gamma)**i);
  }
  return sum;
}

proc kinetic(
    alpha1:real, l1:3*int, xa:3*real,
    alpha2:real, l2:3*int, xb:3*real):real{

  var term0,term1,term2: real;
  term0 = alpha2*(2*(l2(1)+l2(2)+l2(3))+3)*overlap(alpha1,l1,xa,alpha2,l2,xb);
  term1 = -2.0 * alpha2**2 *
    (overlap(alpha1,l1,xa,alpha2,(l2(1)+2,l2(2),l2(3)),xb)
   + overlap(alpha1,l1,xa,alpha2,(l2(1),l2(2)+2,l2(3)),xb)
   + overlap(alpha1,l1,xa,alpha2,(l2(1),l2(2),l2(3)+2),xb));
  term2 = -0.5*(l2(1)*(l2(1)-1)*overlap(alpha1,l1,xa,alpha2,(l2(1)-2,l2(2),l2(3)),xb)
              + l2(2)*(l2(2)-1)*overlap(alpha1,l1,xa,alpha2,(l2(1),l2(2)-2,l2(3)),xb)
              + l2(3)*(l2(3)-1)*overlap(alpha1,l1,xa,alpha2,(l2(1),l2(2),l2(3)-2),xb));
  return term0+term1+term2;
}

proc normalization(alpha:real,powers:3*int):real{
    var (l,m,n):3*int = powers;
    return sqrt((2**(2*(l+m+n)+1.5))*(alpha**(l+m+n+1.5))/
          fact2(2*l-1)/fact2(2*m-1)/fact2(2*n-1)/(PI**1.5));
}


proc grad_nuc_att(
    alpha1:real, l1:3*int, x1:3*real, 
    alpha2:real, l2:3*int, x2:3*real, x3:3*real): 3*real{
                 
    var gamma:real = alpha1+alpha2;

    var xp:3*real = gaussian_product_center(alpha1,x1,alpha2,x2);
    var rab2:real = dist2(x1,x2);
    var rcp2:real = dist2(x3,xp);

    var Ax:[0..l1(1)+l2(1)] real = A_array(l1(1),l2(1),xp(1)-x1(1),xp(1)-x2(1),xp(1)-x3(1),gamma);
    var Ay:[0..l1(2)+l2(2)] real = A_array(l1(2),l2(2),xp(2)-x1(2),xp(2)-x2(2),xp(2)-x3(2),gamma);
    var Az:[0..l1(3)+l2(3)] real = A_array(l1(3),l2(3),xp(3)-x1(3),xp(3)-x2(3),xp(3)-x3(3),gamma);

    var grad_Ax:[0..l1(1)+l2(1)] real = grad_A_array(l1(1),l2(1),xp(1)-x1(1),xp(1)-x2(1),xp(1)-x3(1),gamma);
    var grad_Ay:[0..l1(2)+l2(2)] real = grad_A_array(l1(2),l2(2),xp(2)-x1(2),xp(2)-x2(2),xp(2)-x3(2),gamma);
    var grad_Az:[0..l1(3)+l2(3)] real = grad_A_array(l1(3),l2(3),xp(3)-x1(3),xp(3)-x2(3),xp(3)-x3(3),gamma);
    
    var sum_x:real = 0.0;
    var sum_y:real = 0.0;
    var sum_z:real = 0.0;

    forall (I,J,K) in [0..l1(1)+l2(1),0..l1(2)+l2(2),0..l1(3)+l2(3)]{
      sum_x += grad_Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K+1,rcp2*gamma);
      sum_y += Ax[I]*grad_Ay[J]*Az[K]*Fgamma(I+J+K+1,rcp2*gamma);
      sum_z += Ax[I]*Ay[J]*grad_Az[K]*Fgamma(I+J+K+1,rcp2*gamma);
    }

    var gradV_x:real = -4*PI*exp(-alpha1*alpha2*rab2/gamma)*sum_x;
    var gradV_y:real = -4*PI*exp(-alpha1*alpha2*rab2/gamma)*sum_y;
    var gradV_z:real = -4*PI*exp(-alpha1*alpha2*rab2/gamma)*sum_z;
        
    return (gradV_x,gradV_y,gradV_z);
    
}

proc gaussian_product_center(
    alpha1:real,A:3*real,alpha2:real,B:3*real):3*real{
    var gamma:real = alpha1+alpha2;
    return ((alpha1*A(1)+alpha2*B(1))/gamma,
           (alpha1*A(2)+alpha2*B(2))/gamma,
           (alpha1*A(3)+alpha2*B(3))/gamma);
}

proc grad_A_array(
    l1:int,l2:int,PA:real,PB:real,CP:real,g:real):[0..l1+l2] real{
    //"several pages of algebra were necessary to find this result..."
    var i,r,u,I: int;
    var Imax:int = l1+l2;
    var A: [0..Imax] real;

    forall (i,r,u) in [0..Imax,0..floor(i/2.0):int,0..floor((i-2*r)/2.0):int] {
      I = i-2*r-u;
      A[I] += (i-2.0*r+1.0)/(i-2.0*r-2.0*u+1.0)*CP*
        A_term(i,r,u,l1,l2,PA,PB,CP,g);
    }
    return A;
}
