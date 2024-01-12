/*
We compute the value of H'(q,t) and an approximation of nu'(q,t) of Conjecture 3.4 in the case g=3 and q=53.
We compute an approximation of the Sato Tate density by knowning the first N=100 moments.
*/
q := 53;					// Finite field
g := 3;						// Genus
N := 100;					// Bound for the primes to take into account in computing an approximation to \nu
ComputeActualNumbers := false;			// Do we enumerate the genus-g curves over F_q?

load "ExactFormulaSp.m";

load "STDensity.m";

load "weighted53.m";

load "Densities.m";

load "Printing.m";

"Distance between prediction and measurement", &+[Abs(Densities2[i]-wm53[i]) : i in [1..#wm53]];
"Distance between Sato Tate prediction and measurement", &+[Abs(SmoothDens[i]-wm53[i]) : i in [1..#wm53]];
