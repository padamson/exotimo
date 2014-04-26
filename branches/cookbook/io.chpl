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
