use CGBF;
use basis;
use cints;
use classical;
use LA2;

/*
var sym2powerlist = {
    'S' : [(0,0,0)],
    'P' : [(1,0,0),(0,1,0),(0,0,1)],
    'D' : [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1)],
    'F' : [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),
           (0,3,0),(0,2,1),(0,1,2), (0,0,3)]
    }
*/

proc getints(cgbfs:[] CGBF,nuclei:[] nucleus){
    var nbf:int = cgbfs.numElements;
    var S,h:[1..nbf,1..nbf] real = 0.0;
    (S,h) = get1ints(cgbfs,nuclei);
    //Ints = get2ints(bfs);
    //return S,h,Ints;
    return (S,h);
}

proc get1ints(cgbfs:[] CGBF, nuclei:[] nucleus){
    //Form the overlap S and h=t+vN one-electron Hamiltonian matrices
    var num_nuclei: int = nuclei.numElements;
    var nbf: int = cgbfs.numElements;
    var S,h:[1..nbf,1..nbf] real = 0.0;
    var bfi,bfj: CGBF;
    var nukek: nucleus;

    for i in [1..nbf] {
        bfi = cgbfs(i);
        for j in [1..nbf] {
            bfj = cgbfs(j);
            S(i,j) = bfi.contr_overlap(bfj);
            h(i,j) = bfi.contr_kinetic(bfj);
            for k in [1..num_nuclei] {
              nukek = nuclei(k);
              h(i,j) = h(i,j) + nukek.charge*bfi.contr_nuclear(bfj,nukek.pos);
            }
        }
    }
    return (S,h);
}

proc getT(bfs:[] CGBF){
    //Form the kinetic energy matrix
    var nbf: int = bfs.numElements;
    var T:[1..nbf,1..nbf] real = 0.0;
    var bfi,bfj: CGBF;

    for i in [1..nbf]{
        bfi = bfs(i);
        for j in [1..nbf]{
            bfj = bfs(j);
            T(i,j) = bfi.contr_kinetic(bfj);
        }
    }
    return T;
}

proc getS(bfs:[] CGBF){
    //Form the overlap matrix
    var nbf: int = bfs.numElements;
    var S:[1..nbf,1..nbf] real = 0.0;
    var bfi,bfj: CGBF;

    for i in [1..nbf]{
        bfi = bfs(i);
        for j in [1..nbf]{
            bfj = bfs(j);
            S(i,j) = bfi.contr_overlap(bfj);
        }
    }
    return S;
}

proc getV(bfs:[] CGBF, nuclei:[] nucleus){
    //Form the nuclear attraction matrix V
    var nbf: int = bfs.numElements;
    var num_nuclei: int = nuclei.numElements;
    var V:[1..nbf,1..nbf] real = 0.0;
    var bfi,bfj: CGBF;
    var nukek: nucleus;

    for i in [1..nbf] {
        bfi = bfs(i);
        for j in [1..nbf] {
            bfj = bfs(j);
            for k in [1..num_nuclei] {
              nukek = nuclei(k);
              V(i,j) = V(i,j) + nukek.charge*bfi.contr_nuclear(bfj,nukek.pos);
            }
        }
    }
    return V;
}

proc get2ints(bfs: [] CGBF){
    /*Store integrals in a long array in the form (ij|kl) (chemists
     *notation. We only need i>=j, k>=l, and ij <= kl
     */
    //from array import array
    var nbf = bfs.numElements;
    var totlen:int = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8;
    var Ints:[1..totlen] real = 0.0;
    var ij,kl,ijkl:int;
    writeln("totlen=",totlen);
    for i in [0..nbf-1]{
        for j in [0..i]{
            ij = i*(i+1)/2+j;
            for k in [0..nbf-1]{
                for l in [0..k]{
                    kl = k*(k+1)/2+l;
                    if ij >= kl then {
                        ijkl = ijkl2intindex(i,j,k,l)+1;
                        Ints(ijkl) = CGBF_coulomb(bfs(i+1),bfs(j+1),bfs(k+1),bfs(l+1));
                        writeln("i,j,k,l,ijkl,Ints=",
                            i+1," ",j+1," ",k+1," ",l+1," ",
                            ijkl," ",Ints(ijkl));
                    }
                }
            }
        }
    }
    return Ints;
}

proc getJ(Ints:[] real,D:[] real){
    //Form the Coulomb operator corresponding to a density matrix D
    var nbf:int = D.numElements;
    var d: domain(1) = [1..nbf*nbf];
    var D1d = reshape(D,d); //1D version of Dens
    var J: [1..nbf,1..nbf] real = 0.0;
    var kl,ijkl:int;
    var temp: [1..nbf*nbf] real;
    for i in [0..nbf-1] {
        for j in [0..i] {
            temp = 0.0;
            kl = 0;
            for k in [0..nbf-1]{
                for l in [0..nbf-1]{
                    ijkl = ijkl2intindex(i,j,k,l)+1;
                    temp(kl) = Ints(ijkl);
                    kl += 1;
                }
            }
            J(i,j) = dot(temp,D1d);
            J(j,i) = J(i,j);
        }
    }
    return J;
}

proc getK(Ints:[] real,D:[] real){
    //Form the exchange operator corresponding to a density matrix D
    var nbf:int = D.numElements;
    var d: domain(1) = [1..nbf*nbf];
    var D1d = reshape(D,d); //1D version of Dens
    var K: [1..nbf,1..nbf] real = 0.0;
    var kl,ijkl1,ijkl2:int;
    var temp: [1..nbf*nbf] real;
    for i in [0..nbf-1] {
        for j in [0..i] {
            temp = 0.0;
            kl = 0;
            for k in [0..nbf-1]{
                for l in [0..nbf-1]{
                    ijkl1 = ijkl2intindex(i,k,j,l)+1;
                    ijkl2 = ijkl2intindex(i,l,k,j)+1;
                    temp[kl] = 0.5*(Ints[ijkl1]+Ints[ijkl2]);
                    kl += 1;
                }
            }
            K(i,j) = dot(temp,D1d);
            K(j,i) = K[i,j];
        }
    }
    return K;
}

proc get2JmK(Ints:[] real,D:[] real){
    //Form the 2J-K integrals corresponding to a density matrix D"
    var nbf:int = D.numElements;
    var d: domain(1) = [1..nbf*nbf];
    var D1d = reshape(D,d); //1D version of Dens
    var G: [1..nbf,1..nbf] real = 0.0;
    var kl,ijkl0,ijkl1,ijkl2:int;
    var temp: [1..nbf*nbf] real;
    for i in [0..nbf-1] {
        for j in [0..i] {
            temp = 0.0;
            kl = 0;
            for k in [0..nbf-1]{
                for l in [0..nbf-1]{
                    ijkl0 = ijkl2intindex(i,j,k,l)+1;
                    ijkl1 = ijkl2intindex(i,k,j,l)+1;
                    ijkl2 = ijkl2intindex(i,l,k,j)+1;
                    temp[kl] = 2.0*Ints[ijkl0]-0.5*Ints[ijkl1]-0.5*Ints[ijkl2];
                    kl += 1;
                }
            }
            G[i,j] = dot(temp,D1d);
            G[j,i] = G[i,j];
        }
    }
    return G;
}
