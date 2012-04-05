_extern proc cblas_dasum(
    const N:int, 
    inout X:real(64), 
    const incX:int
    ): real(64);

_extern proc cblas_daxpy(
    const N:int, 
    const alpha:real(64), 
    inout X:real(64),
    const incX:int, 
    inout Y:real(64), 
    const incY:int
    );

_extern proc cblas_dcopy(
    const N:int, 
    inout X:real(64), 
    const incX:int,
    inout Y:real(64), 
    const incY:int
    );

_extern proc cblas_ddot(
    const N:int,
    inout X:real(64),
    const incX:int,
    inout Y:real(64),
    const incY:int
    ): real(64);

_extern proc cblas_dnrm2(
    const N:int, 
    inout X:real(64), 
    const incX:int
    ):real(64);

_extern proc cblas_drot(
    const N:int, 
    inout X:real(64), 
    const incX:int,
    inout Y:real(64), 
    const incY:int, 
    const c:real(64), 
    const s:real(64)
    );

_extern proc cblas_drotg(
    inout a:real(64), 
    inout b:real(64), 
    inout c:real(64), 
    inout s:real(64)
    );

_extern proc cblas_dscal(
    const N:int, 
    const alpha:real(64), 
    inout X:real(64), 
    const incX:int);

_extern proc cblas_dswap(
    const N:int, 
    inout X:real(64), 
    const incX:int, 
    inout Y:real(64), 
    const incY:int);

_extern proc cblas_idamax(
    //uses C-style indexing (starts at zero)
    const N:int, 
    inout X:real(64), 
    const incX:int
    ):int;

//First argument specifies ordering for the matrix
//(Order = 101 for Row-Major; Order = 102 for Col-Major)
_extern proc cblas_dger(
    const Order=101, 
    const M:int, 
    const N:int,
    const alpha:real(64), 
    inout X:real(64), 
    const incX:int,
    inout Y:real(64), 
    const incY:int, 
    inout A:real(64), 
    const lda:int);

//Argument 'Order' specifies Row-Major (101) or Col-Major (102) ordering for the matrix
//Argument 'Uplo' specifies whether the matrix is an upper (121) or lower (122) triangular matrix
//Argument 'TransA' specifies TRANS specifies the operation to be performed as follows:
//
//              TRANS = 111   x := A*x.
//
//              TRANS = 112   x := Transpose(A)*x.
//
//              TRANS = 113   x := ConjugateTranspose(A)*x.
//
//              TRANS = 114   x := Conjugate(A)*x.
//
//Argument 'Diag' specifies whether the matrix is not unit triangular (131) or unit triangular (132)
_extern proc cblas_dtrmv(
    const Order=101, 
    const Uplo=121,
    const TransA=111, 
    const Diag=131,
    const N:int, 
    inout A:real(64), 
    const lda:int,
    inout X:real(64), 
    const incX:int
    );

//Argument 'Order' specifies Row-Major (101) or Col-Major (102) ordering for the matrix
//Argument 'TransA' specifies TRANS specifies the operation to be performed as follows:
//
//              TRANS = 111   x := A*x.
//
//              TRANS = 112   x := Transpose(A)*x.
//
//              TRANS = 113   x := ConjugateTranspose(A)*x.
//
//              TRANS = 114   x := Conjugate(A)*x.
//
_extern proc cblas_dgemv(
    const Order=101, 
    const TransA=111, 
    const M:int, 
    const N:int,
    const alpha:real(64), 
    inout A:real(64), 
    const lda:int,
    inout X:real(64), 
    const incX:int, 
    const beta:real(64),
    inout Y:real(64), 
    const incY:int
    );

_extern proc cblas_dspmv(
    const Order=101, 
    const Uplo=121,
    const N:int, 
    const alpha:real(64), 
    inout Ap:real(64),
    inout X:real(64), 
    const incX:int,
    const beta:real(64), 
    inout Y:real(64), 
    const incY:int
    );

_extern proc cblas_dspr2(
    const Order=101, 
    const Uplo=121,
    const N:int, 
    const alpha:real(64), 
    inout X:real(64),
    const incX:int, 
    inout Y:real(64), 
    const incY:int, 
    inout A:real(64)
    );

//Order: Specifies row-major (C) or column-major (Fortran) data ordering.
//TransA: Specifies whether to transpose matrix A.
//TransB: Specifies whether to transpose matrix B.
//M: Number of rows in matrices A and C (after being transposing A, if applicable).
//N: Number of columns in matrices B and C (after transposing B, if applicable).
//K: Number of columns in matrix A; number of rows in matrix B (after transposing A and/or B, if applicable).
//alpha: Scaling factor for the product of matrices A and B.
//A: Matrix A.
//lda: The size of the first dimension of matrix A; if you are passing a matrix A[m][n], the value should be m.
//B: Matrix B.
//ldb: The size of the first dimension of matrix B; if you are passing a matrix B[m][n], the value should be m.
//beta: Scaling factor for matrix C.
//C: Matrix C.
//ldc: The size of the first dimention of matrix C; if you are passing a matrix C[m][n], the value should be m.
//
//Discussion
//This function multiplies A * B and multiplies the resulting matrix by alpha. 
//It then multiplies matrix C by beta. It stores the sum of these two products in matrix C.
//
//Thus, it calculates either
//
//C←αAB + βC
//
//or
//
//C←αBA + βC
//
//with optional use of transposed forms of A, B, or both.

_extern proc cblas_dgemm(
    const Order=101, 
    const TransA=111, 
    const TransB=111, 
    const M:int, 
    const N:int,
    const K:int, 
    const alpha:real(64), 
    inout A:real(64),
    const lda:int, 
    inout B:real(64), 
    const ldb:int,
    const beta:real(64), 
    inout C:real(64), 
    const ldc:int
    );

//Side: specifies if the coefficient matrix is on the left (141) or right (142)
//If Side is 141, solves A*X=alpha*B or A'*X=alpha*B, depending on TransA.
//If Side is 142, solves X*A=alpha*B or X*A'=alpha*B, depending on TransA.
//Argument 'Diag' specifies whether the matrix is not unit triangular (131) or unit triangular (132)
_extern proc cblas_dtrmm(
    const Order=101, 
    const Side=141,
    const Uplo=121,
    const TransA=111, 
    const Diag=131,
    const M:int, 
    const N:int,
    const alpha:real(64), 
    inout A:real(64), 
    const lda:int,
    inout B:real(64), 
    const ldb:int
    );

//Side: specifies if the coefficient matrix is on the left (141) or right (142)
//If Side is 141, solves A*X=alpha*B or A'*X=alpha*B, depending on TransA.
//If Side is 142, solves X*A=alpha*B or X*A'=alpha*B, depending on TransA.
_extern proc cblas_dtrsm(
    const Order=101, 
    const Side=141,
    const Uplo=121,
    const TransA=111, 
    const Diag=131,
    const M:int, 
    const N:int,
    const alpha:real(64), 
    inout A:real(64), 
    const lda:int,
    inout B:real(64), 
    const ldb:int
    );

_extern proc cblas_dsyrk(
    const Order=101, 
    const Uplo=121,
    const Trans=111, 
    const N:int, 
    const K:int,
    const alpha:real(64), 
    inout A:real(64), 
    const lda:int,
    const beta:real(64), 
    inout C:real(64), 
    const ldc:int
    );


