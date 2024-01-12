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

s := &+multiplicities;
multTrue := [ RF!x/s : x in multiplicities ];
/*
multTrue[i] is the value of H'(q,t) as defined in Equation (10) for t=i-Floor(4*Sqrt(q)).
*/
end if;
