use cints;

class PGBF{
  var exponent: real;
  var origin: 3*real;
  var momentum: string;
  var powers: 3*int = (0,0,0);
  var norm: real = normalization(exponent,powers);
  var coef: real = 1.0;

  proc prim_overlap(other: PGBF):real{
    return this.norm*other.norm*
      overlap(this.exponent, this.powers, this.origin, 
          other.exponent, other.powers, other.origin);
  }
    
  proc prim_kinetic(other: PGBF):real{
    return this.norm*other.norm*
      kinetic(this.exponent, this.powers, this.origin,
          other.exponent, other.powers, other.origin);
  }

  proc prim_nuclear(other: PGBF, C: 3*real):real{
    return this.norm*other.norm*
      nuclear_attraction(this.exponent, this.powers, this.origin,
          other.exponent, other.powers, other.origin, C);
  }

  proc prim_nuclear_gradient(other: PGBF, C: 3*real):3*real{
    var (grad1,grad2,grad3):3*real = grad_nuc_att(
          this.exponent, this.powers, this.origin,
          other.exponent, other.powers, other.origin, C);
    return (
        this.norm*other.norm*grad1,
        this.norm*other.norm*grad2,
        this.norm*other.norm*grad3);
  }

  proc prim_amp(x:3*real):real{
        //"Compute the amplitude of the PGBF at point x,y,z"
        return this.norm * this.coef
               * (x(1)-this.origin(1))**this.powers(1) 
               * (x(2)-this.origin(2))**this.powers(2) 
               * (x(3)-this.origin(3))**this.powers(3)
               * exp(-1.0 * this.exponent * dist2(x,this.origin));
  } 

  proc prim_move_center(dx,dy,dz){
    this.origin = (
        this.origin(1)+dx,
        this.origin(2)+dy,
        this.origin(3)+dz);
    return;
  }

  proc prim_laplacian(pos:3*real):real{
        var amp:real = this.prim_amp(pos);
        //writeln("In prim_laplacian.\n");
        //writeln("amp = ",amp,"\n");
        var alpha:real = this.exponent;
        var x:real = pos(1)-this.origin(1);
        if (x == 0.0) then
          halt("x == 0.0 in prim_laplacian results in divide-by-zero error");
        var y:real = pos(2)-this.origin(2);
        if (y == 0.0) then
          halt("y == 0.0 in prim_laplacian results in divide-by-zero error");
        var z:real = pos(3)-this.origin(3);
        if (z == 0.0) then
          halt("z == 0.0 in prim_laplacian results in divide-by-zero error");
        var x2:real = x*x;
        var y2:real = y*y;
        var z2:real = z*z;
        var r2:real = x2+y2+z2;
        var (L,M,N):3*int = this.powers;
        //writeln("(L,M,N)=",(L,M,N),"\n");
        var term:real = (L*(L-1)/x2 + M*(M-1)/y2 + N*(N-1)/z2) +
                4*alpha*alpha*r2 - 2*alpha*(2*(L+M+N)+3);
        //writeln("term=",term,"\n");
        return this.norm*this.coef*amp*term;
  }

  proc prim_grad(x:real,y:real,z:real):3*real{
        var alpha:real = this.exponent;
        var (I,J,K):3*int = this.powers;
        var C:real = this.norm*this.coef;
        var (x0,y0,z0):3*real = this.origin;
        var fx:real = (x-x0)**I * exp(-1.0*alpha*(x-x0)**2);
        var fy:real = (y-y0)**J * exp(-1.0*alpha*(y-y0)**2);
        var fz:real = (z-z0)**K * exp(-1.0*alpha*(z-z0)**2);
        var gx:real = -2.0*alpha*(x-x0)*fx;
        var gy:real = -2.0*alpha*(y-y0)*fy;
        var gz:real = -2.0*alpha*(z-z0)*fz;
        if (I > 0) then gx += (x-x0)**(I-1) * exp(-1.0*alpha*(x-x0)**2);
        if (J > 0) then gy += (y-y0)**(J-1) * exp(-1.0*alpha*(y-y0)**2);
        if (K > 0) then gz += (z-z0)**(K-1) * exp(-1.0*alpha*(z-z0)**2);
        return (C*gx*fy*fz,C*fx*gy*fz,C*fx*fy*gz);
  }
}

proc PGBF_coulomb(gA:PGBF,gB:PGBF,gC:PGBF,gD:PGBF):real{
    //Coulomb interaction between four cartesian Gaussians; THO eq. 2.22
    return coulomb_repulsion(
        gA.origin(1),gA.origin(2),gA.origin(3),gA.norm,
        gA.powers(1),gA.powers(2),gA.powers(3),gA.exponent,
        gB.origin(1),gB.origin(2),gB.origin(3),gB.norm,
        gB.powers(1),gB.powers(2),gB.powers(3),gB.exponent,
        gC.origin(1),gC.origin(2),gC.origin(3),gC.norm,
        gC.powers(1),gC.powers(2),gC.powers(3),gC.exponent,
        gD.origin(1),gD.origin(2),gD.origin(3),gD.norm,
        gD.powers(1),gD.powers(2),gD.powers(3),gD.exponent);
}

proc PGBF_three_center(gA:PGBF,gB:PGBF,gC:PGBF):real{
    //Three-center integral between Gaussians
  var na,nb,nc,ix,iy,iz:real;
    na = gA.norm;
    nb = gB.norm;
    nc = gC.norm;
    ix = three_center_1D(gA.origin(1),gA.powers(1),gA.exponent,
                         gB.origin(1),gB.powers(1),gB.exponent,
                         gC.origin(1),gC.powers(1),gC.exponent);
    iy = three_center_1D(gA.origin(2),gA.powers(2),gA.exponent,
                        gB.origin(2),gB.powers(2),gB.exponent,
                        gC.origin(2),gC.powers(2),gC.exponent);
    iz = three_center_1D(gA.origin(3),gA.powers[2],gA.exponent,
                        gB.origin(3),gB.powers[2],gB.exponent,
                        gC.origin(3),gC.powers[2],gC.exponent);
    return na*nb*nc*ix*iy*iz;
}
