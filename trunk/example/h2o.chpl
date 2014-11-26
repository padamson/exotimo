use io;
use classical;
use Help;
use CGBF;

//var Dmat: domain(1) = {1};
//var Dnuclei: domain(1) = {1..1}; 
//var nuclei:[Dnuclei] nucleus;

var Dnuclei: domain(1) = {1..1};
var nuclei: [Dnuclei] nucleus;
var Dbasis: domain(1) = {1..1};
var basis: [Dbasis] CGBF;

proc main(args: [] string) {
  var infile: string;
  var aLength: int;
  for a in args {
    if a=="--help" {
      printUsage();
      writeln("\nINPUT:");
      writeln(  "======");
      //writeln("You can provide an arbitrary list of filenames.");
      writeln("  -i==<filename>    :input filename");
      exit(0);
    } else if a.substring(1..2)=="-i" {
      aLength = a.length;
      infile = a.substring(5..aLength);
      writeln("\n======PARSING INPUT======");
      writeln("Input file:",infile,"\n");
    }
  }

  readIn(infile);

}
