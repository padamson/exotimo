use constants;

proc gtprd(A:[MATRIX_SIZE] real, B:[MATRIX_SIZE] real, R:[MATRIX_SIZE] real,
	n: int, m: int, l: int){
	var k, ik, j, ir, ij, ib: int;
	ir = 0;
	ik = -n;
	for k in [1..l]{
		ij = 0;
		ik = ik + m;
		for j in [1..m]{
			ir = ir + 1;
			ib = ik;
			R[ir] = ZERO;
			for i in [1..n]{
				ij = ij + 1;
				ib = ib + 1;
				R[ir] = R[ir] + A[ij] * B[ib];
			}
		}
	}
	return;
}
