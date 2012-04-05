use Norm;
use constants;
use eigen;

//import logging
//from math import sqrt
//from NumWrap import matrixmultiply,transpose,diagonal,identity,zeros,eigh,cholesky,inv

/*
 * Note: to be really smart in a quantum chemistry program, we would
 * want to only symmetrically orthogonalize the S matrix once, since
 * the matrix doesn't change during the SCF procedure. Thus, we would
 * want to call X = SymOrth(S), and then geighD(H,X), rather
 * than calling geigh(H,S) every time, since the latter
 * recomputes the symmetric orthogonalization every SCF cycle.
*/
proc identity(D:domain):[D] real{
  var A:[D] real = 0.0;
  forall i in D.dim(1){
    A(i,i) = 1.0;
  }
  return A;
}

proc trace2(A:[] real,B:[] real):real{
  if A.domain != B.domain then {
    writeln("Matrices have unequal domains in function trace2.");
    halt();
  }
  const D: domain(2) = [A.domain];
  var sum1,sum2:real = 0.0;
  forall i in D.dim(1) {
    forall j in D.dim(2) {
      sum1 += A(i,j)*B(j,i);
    }
    sum2 += sum1;
  }
  return sum2;
}

proc dot(A:[] real,B:[] real):real{
  if A.numElements != B.numElements then {
    writeln("Arrays are of unequal length in function dot.");
    halt();
  }
  var sum:real = 0.0;
  forall i in [A.domain] {
    sum += A(i)*B(i);
  }
  return sum;
}

/*
proc norm(a:[] real):real{
    //val = norm(vec) : Return the 2-norm of a vector
    return sqrt(dot(a,a));
}
*/

proc transpose(A:[] real) {
  const D: domain(2) = [A.domain];
  var B: [D.dim(2), D.dim(1)] real; 
  forall i in D.dim(1) {
    forall j in D.dim(2) {
      writeln("setting B(",i,",",j,") to ",A(j,i));
      B(i,j) = A(j,i);
    }
  }
  return B;
}

proc sym(A:[]):[A.domain] real {
    //B = sym(A) : Symmetrize a matrix
    return 0.5*(A+transpose(A));
}

proc matrixmultiply(A:[],B:[]):[A.domain] real {
  //if A.domain != B.domain then {
    //writeln("Matrices have unequal domains in function matrixmultiply.");
    //halt();
  //}
  var AB: [A.domain] real = 0.0;
  const D: domain(2) = [A.domain];
  forall i in D.dim(1) {
    forall j in D.dim(2) {
      forall k in D.dim(1) {
        AB(i,j) += A(i,k)*B(k,j);
      }
    }
  }
  return AB;
}

proc simx(A:[] real,B:[] real,trans: string="N"){
  /*
   * C = simx(A,B,trans)
   * Perform the similarity transformation C = B'*A*B (trans='N') or
   * C = B*A*B' (trans='T').
    */
    if trans=="T" then {
      return matrixmultiply(B,matrixmultiply(A,transpose(B)));
    } else {
      return matrixmultiply(transpose(B),matrixmultiply(A,B));
    }
}

proc outprod(A:[] real):[A.domain] real {
    //D = outprod(A) : Return the outer product A*A'
    return matrixmultiply(A,transpose(A));
}

proc geigh1(H:[] real,A:[] real,orthog:string="Sym"){
    //Generalized eigenproblem using a symmetric matrix H.
    //This one forms the canonical transformation from S.
    //
    //Option:
    //orthog     'Sym'   Use Symmetric Orthogonalization (default)
               //'Can'   Use Canonical Orthogonalization
               //'Chol'  Use a Cholesky decomposition (not currently implemented)
                //
    //
    var X:[H.domain] real = 0.0;
    const D: domain(2) = [H.domain];
    var val: [D.dim(1)] real;
    var vec: [H.domain] real;
    if orthog == 'Can' then {
      X = CanOrth(A);
    } else if orthog == 'Chol' then {
      //   X = CholOrth(A);
      writeln("Cholesky decomposition not currently implemented in 'geigh'.");
      halt();
    }
    else {
      X = SymOrth(A);
    }
    (val,vec) = geigh2(H,X);
    return (val,vec);
}

proc geigh2(H:[] real,A:[] real){
    //Generalized eigenproblem using a symmetric matrix H.
    //This one takes the canonical transformation matrix as A.
    const D: domain(2) = [H.domain];
    var val: [D.dim(1)] real;
    var vec: [H.domain] real;
    (val,vec) = eigh(simx(H,A));
    vec = matrixmultiply(A,vec);
    return (val,vec);
}

proc SymOrth(S){
    //Symmetric orthogonalization of the real symmetric matrix S.
    //This is given by Ut(1/sqrt(lambda))U, where lambda,U are the
    //eigenvalues/vectors.
    const D: domain(2) = [S.domain];
    var val: [D.dim(1)] real;
    var vec: [S.domain] real;
    (val,vec) = eigh(S);
    var shalf: [S.domain] real = identity(D);
    forall i in D.dim(1){ 
        shalf(i,i) /= sqrt(val(i));
    }
    return simx(shalf,vec,'T');
}

proc CanOrth(S){
    //Canonical orthogonalization of matrix S. This is given by
    //U(1/sqrt(lambda)), where lambda,U are the eigenvalues/vectors.
    const D: domain(2) = [S.domain];
    var val: [D.dim(1)] real;
    var vec: [S.domain] real;
    (val,vec) = eigh(S);
    forall i in D.dim(1){ 
        vec[..,i] = vec[..,i] / sqrt(val[i]);
    }
    return vec;
}

proc cholesky(A:[] real){
  var sum: real;
  const D: domain(2) = [A.domain];
  var U: [A.domain] real;

  for row in [1..D.high(1)] {
    /* First compute U[row,row] */
    sum = A[row,row];
    for j in [1..row-1] {
      sum -= U[j,row]*U[j,row];
    }
    if (sum > TOOSMALL) then {
      U[row,row] = sqrt(sum);
      /* Now find elements U[row,ck], k > row. */
      for k in [row+1..D.high(1)] {
        sum = A[row,k];
        for j in [1..row-1] {
          sum -= U[j,row]*U[j,k];
        }
        U[row,k] = sum/U[row,row];
      }
    } else { /* blast off the entire row. */
      for k in [row..D.high(1)] {
        U[row,k] = 0.0;
      }
    }
  }
  return U;
}

/*
proc CholOrth(S){
  //Cholesky orthogonalization of matrix X. This is given by
  //LL^T=S; X = transpose(inv(L))
  return transpose(inv(cholesky(S)));
}
*/

/*
// SimilarityTransform/T are deprecated in favor of simx
proc SimilarityTransformT(H,X):
    "Return the similarity transformation XtHX of H"
    logging.warning("SimilarityTransformT is deprecated: use simx")
    return simx(H,X)
    #return matrixmultiply(transpose(X),matrixmultiply(H,X))

proc SimilarityTransform(H,X): 
    "Return the transpose similarity transformation XHXt of H"
    logging.warning("SimilarityTransform is deprecated: use simx")
    return simx(H,X,'T')
    #return matrixmultiply(X,matrixmultiply(H,transpose(X)))

proc mkdens(c,nstart,nstop){
    //Form a density matrix C*C' given eigenvectors C[nstart:nstop,:]
    return matrixmultiply(c[..,nstart..nstop],transpose(c[..,nstart..nstop]));
}

proc mkdens_spinavg(c,nclosed,nopen){
    //Form a spin averaged density matrix with *nclosed* closed
    //shell orbitals and *nopen* open shell orbitals
    return mkdens(c,0,nclosed) + 0.5*mkdens(c,nclosed,nclosed+nopen);
}
*/

/*
#added by Hatem H Helal 18.07.2007
proc pad_out(matrix):
    #this will make debugging matrix operations easier by getting rid of 
    #matrix elements which are ridiculously tiny
    
    print array2string(matrix,max_line_width=200,precision=7,suppress_small=true);    print "\n\n"

    return 0

proc diagonal_mat(vector):
    #takes a vector with N components and returns an NXN matrix whose diagonal elements
    #are the components of the vector
    len = vector.shape[0]
        
    matrix = zeros((len,len),'d')
    
    i=0
    for i in range(len):
        matrix[i,i]=vector[i]

    return matrix    
    
*/
