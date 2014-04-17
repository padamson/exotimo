use constants;

proc scfR(C:[MATRIX_SIZE] real, R:[MATRIX_SIZE] real, m: int, nocc: int){
	var sum: real;
	var i, j, k, ij, ji, kk, ik, jk: int;
	for i in [1..m] {
		for j in [1..i]{
			sum = ZERO;
			for k in [1..nocc] {
				kk = m*(k-1); 
				ik = kk + i; 
				jk = kk + j;
				sum = sum + C[ik]*C[jk];
			}
			ij = m*(j-1) + i; 
			ji = m*(i-1) + j;
			R[ij] = sum;
			R[ji] = sum;
		}
	}
	return;
}
