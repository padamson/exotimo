use constants;

proc main(){

	var H, HF, C, V, R, Cbar, epsilon, Rold:[ARB]: real;
	var E, Rsum, term: real;
	var m, n, i, nfile, icon, mm: int;

	mm = m*m;
	R, Rold = ZERO;

	do {
		E = ZERO;
		icon = 0;
		HF = H;
		for i in 1..mm{
			E = E + R[i] + HF[i];
		}
		call scfGR(HF, R, m, nfile);
		for i in 1..mm{
			E = E + R[i] + HF[i];
		}
		writeln(" Current Electronic Energy = ", E);
		call gtprd(V, HF, R, m, m, m);
		call gmprd(R, V, HF, m, m, m);
		call eigen(HF, Cbar, m);
		for i in 1..n{
			epsilon[i] = HF[m*(i-1)+i];
		}
		call gmprd(V, Cbar, C, m, m, m);
		call scfR(C, R, m, n);
		Rsum = ZERO;
		for i in 1..mm {
			term = abs(R[i] - Rold[i]);
			Rold[i] = R[i];
			term > crit ? icon++;
			Rsum = Rsum + term;
		}
		writeln(" Sum of differences in R = ", Rsum, icon, " Changing");
	} while (icon > 0);
}
