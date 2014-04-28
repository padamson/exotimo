class nucleus{
  var charge: real;
  var pos: 3*real;
}

class classical{
  var num_nuclei: int;
  var nuclei_dom= [1..num_nuclei];
  var nuclei: [nuclei_dom] nucleus;

  proc add_nucleus(charge:real,pos:3*real){
    var temp_nucleus = new nucleus(charge=charge,pos=pos);
    this.num_nuclei += 1;
    this.nuclei_dom = [1..this.num_nuclei];
    this.nuclei(num_nuclei) = temp_nucleus;
  }
}
