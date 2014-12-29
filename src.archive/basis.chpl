use CGBF;
const sPowerList: 3*int = (0,0,0);
const pPowerList: [1..3] 3*int = ((1,0,0),(0,1,0),(0,0,1));
const dPowerList: [1..6] 3*int = ((2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1));
const fPowerList: [1..10] 3*int = 
    ((3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),(0,3,0),(0,2,1),(0,1,2),(0,0,3));

class basis{
  var particle: string = "electron";
  var first_num_prims:int = 2;
  var first_prims_dom = [1..first_num_prims];
  var first_pexps: [first_prims_dom] real = (1.0,0.1);
  var first_pcoefs: [first_prims_dom] real = (2.0,1.0);
  var first_origin: 3*real = (0.0,0.0,0.0);
  var first_powers: 3*int = (0,0,0);
  var first_atid: int = 1;
  var num_cgbfs:int = 1;
  var cgbfs_dom = [1..num_cgbfs];
  var cgbfs: [cgbfs_dom] CGBF = initialize_cgbfs(first_origin,
      first_powers,first_atid,first_num_prims,first_pexps,first_pcoefs);

  proc initialize_cgbfs(origin:3*real,powers:3*int,atid:int,
      num_prims:int, pexps: [1..num_prims] real,
      pcoefs: [1..num_prims] real):[this.cgbfs_dom] CGBF{
   var temp_cgbf = new CGBF(
       origin=origin,powers=powers,atid=atid,first_pexp=pexps(1),
       first_pcoef=pcoefs(1));
   for i in 2..num_prims {
     temp_cgbf.add_primitive(exponent=pexps(i),coefficient=pcoefs(i));
   }
   temp_cgbf.normalize;
   return temp_cgbf;
  }

  proc add_CGBF(origin:3*real,powers:3*int,atid:int,
      num_prims:int, pexps: [1..num_prims] real,
      pcoefs: [1..num_prims] real){
   var temp_cgbf = new CGBF(
       origin=origin,powers=powers,atid=atid,first_pexp=pexps(1),
       first_pcoef=pcoefs(1));
   for i in 2..num_prims {
     temp_cgbf.add_primitive(exponent=pexps(i),coefficient=pcoefs(i));
   }
   temp_cgbf.normalize;
   this.num_cgbfs += 1;
   this.cgbfs_dom = [1..this.num_cgbfs];
   this.cgbfs(num_cgbfs) = temp_cgbf;
  }

  proc readBasis(basisFileName:string="/opt/ExotiMO/Basis/sto3g.dat",
      atom:string="H",basis:string="STO3G",
      origin:3*real=(0.0,0.0,0.0),atid:int=1){
    var basisFile = new file(basisFileName, FileAccessMode.read);
    basisFile.open();
    var atomFromFile: string;
    var basisFromFile: string;
    var symFromFile: string;
    var numCGBFsFromFile: int;
    var numPGBFsFromFile: int;
    atomFromFile=basisFile.read(string);
    basisFromFile=basisFile.read(string);
    numCGBFsFromFile=basisFile.read(int);
    while atomFromFile != atom | basisFromFile != basis do {
      [i in 1..numCGBFsFromFile] {
        symFromFile=basisFile.read(string);
        numPGBFsFromFile=basisFile.read(int);
        [j in 1..numPGBFsFromFile] basisFile.readln(string);
      }
      basisFile.readln(string);
      atomFromFile=basisFile.read(string);
      basisFromFile=basisFile.read(string);
      numCGBFsFromFile=basisFile.read(int);
    }
    writeln("Found atom ",atom," and basis ",basis,
        " with ",numCGBFsFromFile," CGBFs.");
  for j in 1..numCGBFsFromFile {
    symFromFile=basisFile.read(string);
    numPGBFsFromFile=basisFile.read(int);
    writeln("\nAdding ",symFromFile," type CGBF with ",numPGBFsFromFile,
        " PGBFs to ", this.particle, " basis set.");
    var exponentsFromFile: [1..numPGBFsFromFile] real;
    var coefficientsFromFile: [1..numPGBFsFromFile] real;
    for i in 1..numPGBFsFromFile{
      basisFile.read(int);
      exponentsFromFile(i) = basisFile.read(real);
      coefficientsFromFile(i) = basisFile.read(real);
    }
    select symFromFile{
      when "S" {
        this.add_CGBF(origin=origin,powers=sPowerList,atid=atid,
            num_prims=numPGBFsFromFile,
            pexps=exponentsFromFile,
            pcoefs=coefficientsFromFile);
      }
      when "P" {
        for i in 1..3 {
          this.add_CGBF(origin=origin,powers=pPowerList(i),atid=atid,
              num_prims=numPGBFsFromFile,
              pexps=exponentsFromFile,
              pcoefs=coefficientsFromFile);
         }
      }
      when "D" {
        for i in 1..6 {
          this.add_CGBF(origin=origin,powers=dPowerList(i),atid=atid,
              num_prims=numPGBFsFromFile,
              pexps=exponentsFromFile,
              pcoefs=coefficientsFromFile);
        }
      }
      when "F" {
        for i in 1..10 {
          this.add_CGBF(origin=origin,powers=fPowerList,atid=atid,
              num_prims=numPGBFsFromFile,
              pexps=exponentsFromFile,
              pcoefs=coefficientsFromFile);
        }
      }
      otherwise writeln("Symmetry type ", symFromFile," not supported.");
    }
  }
  }
}
