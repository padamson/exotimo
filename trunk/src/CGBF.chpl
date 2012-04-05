use PGBF;
use cints;
use crys;

class CGBF{
  var origin: 3*real;
  var powers: 3*int = (0,0,0);
  var norm: real = 1.0;
  var atid:int = 1;
  var atno:int = 1;
  var num_prims:int = 1;
  var prims_dom = [1..num_prims];
  var pnorms: [prims_dom] real; 
  var first_pexp: real;
  var first_pcoef: real;
  var pexps: [prims_dom] real = (first_pexp);
  var pcoefs: [prims_dom] real = (first_pcoef);
  var prims: [prims_dom] PGBF = initialize_prims(first_pexp,first_pcoef);

  proc initialize_prims(exponent:real,coefficient:real):[this.prims_dom] PGBF{
        //Add a primitive BF to this contracted set
        var pbf = new PGBF(
              exponent=exponent,origin=this.origin,
              powers=this.powers,coef=coefficient);
        this.pnorms(1) = pbf.norm;
        return (pbf);
  }

  proc contr_center(other:CGBF):3*real{
    var (xa,ya,za):3*real = this.origin;
    var (xb,yb,zb):3*real = other.origin;
    return (0.5*(xa+xb),0.5*(ya+yb),0.5*(za+zb));
  }

  proc add_primitive(exponent:real,coefficient:real){
        //Add a primitive BF to this contracted set
        var pbf = new PGBF(
              exponent=exponent,origin=this.origin,
              powers=this.powers,coef=coefficient);
        this.num_prims += 1;
        this.prims_dom = [1..this.num_prims];
        this.prims(num_prims) = pbf;
        this.pnorms(num_prims) = pbf.norm;
        this.pexps(num_prims) = exponent;
        this.pcoefs(num_prims) = coefficient;
        return;
  }

  proc normalize{
        //Normalize the current CGBF
        var olap:real = this.contr_overlap(this);
        this.norm = 1.0/sqrt(olap);
        writeln("Normalized CGBF, norm =", this.norm);
  }

  proc contr_overlap(other:CGBF):real{
        //Overlap matrix element with another CGBF
        var Sij:real = 0.0;
        forall (i,j) in [this.prims_dom, other.prims_dom] do
          Sij += this.pcoefs(i)*other.pcoefs(j)*this.prims(i).prim_overlap(other.prims(j));
        //return this.norm*other.norm*Sij;
        return Sij;
  }

  proc contr_kinetic(other:CGBF):real{
        //KE matrix element with another CGBF"
        var Tij:real = 0.0;
        forall (i,j) in [this.prims_dom, other.prims_dom] do
          Tij += this.pcoefs(i)*other.pcoefs(j)*this.prims(i).prim_kinetic(other.prims(j)); 
        return this.norm*other.norm*Tij;
  }

  proc contr_nuclear(other:CGBF,C: 3*real):real{
        //Nuclear matrix element with another CGBF and a center C
        var Vij: real = 0.0;
        forall (i,j) in [this.prims_dom, other.prims_dom] do
          Vij += this.pcoefs(i)*other.pcoefs(j)*this.prims(i).prim_nuclear(other.prims(j),C); 
        return this.norm*other.norm*Vij;
  }

  proc contr_amp(x:real,y:real,z:real):real{
    //Compute the amplitude of the CGBF at point x,y,z
    var val:real = 0.0;
    for i in [this.prims_dom] do {
      val += this.prims(i).prim_amp(x,y,z);
    }
    return this.norm*val;
  }

  proc contr_move_center(dx,dy,dz){
    //Move the basis function to another center
    this.origin = (this.origin(1)+dx,this.origin(2)+dy,this.origin(3)+dz);
    //for prim in self._prims: prim.move_center(dx,dy,dz)
    return;
  }

}

proc CGBF_coulomb(a0:CGBF,b0:CGBF,c0:CGBF,d0:CGBF):real{
  var a:CGBF = a0;
  var b:CGBF = b0;
  var c:CGBF = c0;
  var d:CGBF = d0;
  var suma: int = a.powers(1)+a.powers(2)+a.powers(3);
  var sumb: int = b.powers(1)+b.powers(2)+b.powers(3);
  var sumc: int = c.powers(1)+c.powers(2)+c.powers(3);
  var sumd: int = d.powers(1)+d.powers(2)+d.powers(3);
  if suma < sumb then (a,b) = (b0,a0);
  if sumc < sumd then (c,d) = (d0,c0);

  return a.norm*b.norm*c.norm*d.norm*contr_coulomb_crys(
        a.num_prims,a.pexps,a.pcoefs,a.pnorms,a.origin,a.powers,
        b.num_prims,b.pexps,b.pcoefs,b.pnorms,b.origin,b.powers,
        c.num_prims,c.pexps,c.pcoefs,c.pnorms,c.origin,c.powers,
        d.num_prims,d.pexps,d.pcoefs,d.pnorms,d.origin,d.powers);
}

proc CGBF_three_center(a:CGBF,b:CGBF,c:CGBF):real{
  var sum:real = 0.0;
  forall (i,j,k) in [1..a.num_prims, 1..b.num_prims, 1..c.num_prims] {
    sum += a.pcoefs(i)*b.pcoefs(j)*c.pcoefs(k)*
      PGBF_three_center(a.prims(i),b.prims(j),c.prims(k));
  }
  return a.norm*b.norm*c.norm*sum;
}


