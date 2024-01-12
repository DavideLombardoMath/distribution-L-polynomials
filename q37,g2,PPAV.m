/*
We compute the value of H'(q,t) and an approximation of nu'(q,t) of Conjecture 5.4 in the case g=2 and q=37.
We compute the value of H'(q,t) when we take account of all PPAV. For more details see Remark 3.9.
We compute an approximation of the Sato Tate density by knowning the first N moments.
*/
q := 37;					// Finite field
g := 2;						// Genus
N := 100;					// Bound for the primes to take into account in computing an approximation to \nu
ComputeActualNumbers := true;			// We enumerate the genus-g curves (and also the geometrically decomposable PPAVs) over F_q

load "ExactFormulaSp.m";

load "STDensity.m";

load "Densities.m";

/*
Count the number of curves with a fixed trace. Count a curve C with multiplicity 1/#Aut(C).
*/
if ComputeActualNumbers then
F := FiniteField(q);
multiplicities := [Rationals()!0 : i in [-Floor(4*Sqrt(q))..Floor(4*Sqrt(q))]];

RF := RealField(10);

for i1 in [1..q] do
for i2 in [1..q] do
                for i3 in [1..q] do
                        twists, gp := TwistsFromG2Invariants( [F!i1, F!i2, F!i3] );
                        for t in twists do
                                a := q+1-#RationalPoints(t);
                                if #gp eq 2 then
                                        RationalAutGroupSize := 2;
                                else
                                        RationalAutGroupSize := #AutomorphismGroup(t);
                                end if;
                                multiplicities[a+Floor(4*Sqrt(q))+1] := multiplicities[a+Floor(4*Sqrt(q))+1] + 1/RationalAutGroupSize;
                        end for;
                end for;
end for;
end for;

multiplicitiesPP := multiplicities;


/*
Next, we count the Weil restrictions of scalars
*/
multiplicitiesPP[Floor(4*Sqrt(q))+1] := multiplicitiesPP[Floor(4*Sqrt(q))+1] + 1/2*(q^2);

/*
And finally, the products of two elliptic curves
*/
for j1 in [0..(q-1)] do
for j2 in [j1..(q-1)] do

	E1 := EllipticCurveFromjInvariant(F!j1);
	E2 := EllipticCurveFromjInvariant(F!j2);

	for E1t in Twists(E1) do
	for E2t in Twists(E2) do
		a1 := q+1 - #RationalPoints(E1t);
		a2 := q+1 - #RationalPoints(E2t);
		multiplicitiesPP[a1+a2+Floor(4*Sqrt(q))+1] := multiplicitiesPP[a1+a2+Floor(4*Sqrt(q))+1] + 1/#AutomorphismGroup(E1)*1/#AutomorphismGroup(E2);
	end for;
	end for;
end for;
end for;




s := &+multiplicities;
multTrue := [ RF!x/s : x in multiplicities ];
sPP := &+multiplicitiesPP;
multTruePP := [ RF!x/s : x in multiplicitiesPP ];
/*
multTrue[i] is the values of H'(q,t) as defined in Equation 10 for t=i-Floor(4*Sqrt(q)).
multTruePP[i] is the values of H'(q,t) taking account also of the PPAV.
*/

end if;

load "Printing.m";

j := 0;
"*** True densities counting PPAVs ***";
for t in [-Floor(2*g*Sqrt(q))..Floor(2*g*Sqrt(q))] do
	j := j+1;
	if multTrue[j] gt 10^(-5) then
		print "(", t, ",", RealField(5)!(multTruePP[j]), ")";
	else
		print "(", t, ",", 0, ")";
	end if;
end for;



"Distance between prediction and measurement", &+[Abs(Densities2[i]-multTrue[i]) : i in [1..#multTrue]];
"Distance between prediction and measurement with PPAV", &+[Abs(Densities2[i]-multTruePP[i]) : i in [1..#multTruePP]];
"Distance between Sato Tate prediction and measurement", &+[Abs(SmoothDens[i]-multTrue[i]) : i in [1..#multTrue]];
