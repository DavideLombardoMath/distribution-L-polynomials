/*
We compute the value of H'(q,t) and an approximation of nu'(q,t) of Conjecture 3.4 in the case g=2.
We compute an approximation of the Sato Tate density by knowning the first N=100 moments.
*/
q := 1009;					// Finite field
g := 2;						// Genus
N := 100;					// Bound for the primes to take into account in computing an approximation to \nu
ComputeActualNumbers := true;			// Do we enumerate the genus-g curves over F_q?

load "ExactFormulaSp.m";

load "STDensity.m";

load "Densities.m";

load "CurveCounting.m";

load "Printing.m";


if ComputeActualNumbers then
	"Distance between prediction and measurement", &+[Abs(Densities2[i]-multTrue[i]) : i in [1..#multTrue]];
end if;
"Distance between Sato Tate prediction and measurement", &+[Abs(SmoothDens[i]-multTrue[i]) : i in [1..#multTrue]];

exit;