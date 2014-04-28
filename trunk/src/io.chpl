use h2o;
use classical;

proc write_matrix(a:[]){
  const D:domain(2) = a.domain;
  const d1:domain(1) = D.dim(1);
  const d2:domain(1) = D.dim(2);
  const high1:int = d1.high;
  const high2:int = d2.high;
  for i in d1 { 
    if (i > 1) then
      write("  ( ");
    else
      write("( ( ");

    for j in d2 {
      if (j == high2) then {
        if (i == high1) then 
          write(", ",a(i,j)," ) )\n");
        else
          write(", ",a(i,j)," ) ,\n");
      }
      else {
        if (j > 1) then
          write(", ",a(i,j));
        else
          write(a(i,j));
      }
    }
  }
}

proc read_matrix(a:[], reader:channel){
  const D:domain(2) = a.domain;
  const d1:domain(1) = D.dim(1);
  const d2:domain(1) = D.dim(2);
  for i in d1 { 
    for j in d2 {
      a[i,j] = reader.read(real);
    }
  }
}

proc readNuclei(reader:channel){
  var numNuclei:int;
  var tempCharge: int;
  var tempPos: 3*real;
  findStringInFile("nuclei", reader);
  numNuclei= reader.read(int);
  writeln("numNuclei=",numNuclei);
  Dnuclei= {1..numNuclei};
  for i in Dnuclei { 
    tempCharge = reader.read(int);
    for j in 1..3 {
      tempPos(j) = reader.read(real);
    }
    nuclei[i] = new nucleus(charge=tempCharge,pos=tempPos);
  }
  writeNuclei();

}

proc write_array(a:[]){
  const D:domain(1) = a.domain;
  const high:int = D.high;
  write("( ");
  for i in D { 
    if (i == high) then
      write(a(i));
    else 
      write(a(i),", ");
  }
  write(" )\n");
}

proc findStringInFile(keyword:string, reader:channel){
  var stringIn: string = "";
  while (stringIn != keyword) do {
    stringIn = reader.read(string);
  }
}

proc readIn(filename:string) {
  var infile = open(filename, iomode.r);
  var reader = infile.reader();

  readNuclei(reader);
  readBasis(reader);

}

proc writeNuclei(){
  writeln("NUCLEI");
  writeln("====================");
  writeln("Charge   X   Y   Z");
  for i in Dnuclei { 
    writeln(nuclei(i).charge," ",nuclei(i).pos(1)," ",nuclei(i).pos(2)," ",nuclei(i).pos(3));
  }
}



