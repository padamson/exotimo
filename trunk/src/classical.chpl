class nucleus{
  var charge: real;
  var pos: 3*real;
}

class classical{
  var first_charge: real = 1.0;
  var first_pos: 3*real = (0.0,0.0,0.0);
  var num_nuclei: int = 1;
  var nuclei_dom= [1..num_nuclei];
  var nuclei: [nuclei_dom] nucleus = initialize_nuclei(first_charge,first_pos);

  proc initialize_nuclei(charge:real,pos:3*real): [this.nuclei_dom] nucleus{
    var temp_nucleus = new nucleus(charge=charge, pos=pos);
    return temp_nucleus;
  }

  proc add_nucleus(charge:real,pos:3*real){
    var temp_nucleus = new nucleus(charge=charge,pos=pos);
    this.num_nuclei += 1;
    this.nuclei_dom = [1..this.num_nuclei];
    this.nuclei(num_nuclei) = temp_nucleus;
  }
}
