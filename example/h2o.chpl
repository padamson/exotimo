use io;
use classical;

//var Dmat: domain(1) = {1};
//var Dnuclei: domain(1) = {1..1}; 
//var nuclei:[Dnuclei] nucleus;

var Dnuclei: domain(1) = {1..1};
var nuclei: [Dnuclei] nucleus;
var Dbasis: domain(1) = {1..1};
var basis: [Dbasis] CGBF;

proc main () {

  readIn("/Users/padamson/Research/exotimo/trunk/example/h2o.in");

}
