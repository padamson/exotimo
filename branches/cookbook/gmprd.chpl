use constants;

proc gmprd(A:[MATRIX_SIZE] real, B:[MATRIX_SIZE] real, R:[MATRIX_SIZE] real,
	n: int, m: int, l: int){
	var k, ik, j, ir, ji, ib: int;
	ir = 0;
	ik = -m;
	for k in 1..l{
		ik = ik + m;
		for j in 1..n{
			ir = ir + 1;
			ji = j - n;
			ib = ik;
			R[ir] = ZERO;
			for i in 1..m{
				ji = ji + n;
				ib = ib + 1;
				R[ir] = R[ir] + A[ji] * B[ib];
			}
		}
	}
	return;
}
