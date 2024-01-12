function SpModuloN(N)
	if IsPrime(N) then
		return Sp(4, N);
	else
		Z := Integers();
		G := GL(4,quo<Z|N>);
		A := elt<G | 1,0,0,0, 1,-1,0,0, 0,0,1,1, 0,0,0,-1>;
		B := elt<G| 0,0,-1,0, 0,0,0,-1, 1,0,1,0, 0,1,0,0>;
		H := sub<G|A,B>;
	return H;
	end if;
end function;

/*
We count the matrices with a fixed trace and multiplier q in GSp_4(Z/ell Z). If ell is prime, we use ExactFormulaSp.m. If ell is not prime, we simply count the relevant matrices directly. Notice that every element of GSp_4(Z/ell Z) with multiplier q can be written as DiagonalMatrix([q,q,1,1]) times an element in Sp_4(Z/ell Z).
*/
function ComputeCorrectionCoefficients(g, ell, q : Truncate := Infinity())

	if (g gt 2) and not IsPrime(ell) then
		error "Not implemented";
	end if;

	if IsPrime(ell) then
		return [C(g, ell, q, t) : t in [0..(ell-1)]], #Sp(4, ell);
	else

		gen := q mod ell;
		M := DiagonalMatrix([gen,gen,1,1]);
		G := SpModuloN(ell);
		Multiplicities := [ 0 : i in [0..(ell-1)] ];
	
		count1 := 0;
		R := RealField(10);
	
		for g in G do
			tr := Integers()!Trace(GL(4,Integers(ell))!M*g);
			Multiplicities[tr+1] := Multiplicities[tr+1] + 1;
			count1 := count1 + 1;
			if count1 gt Truncate then
				break;
			end if;
			if count1 mod 10^6 eq 0 then
				"count =", count1;
				Multiplicities, count1;
				[RealField(10)!m/count1 : m in Multiplicities];
			end if;
		end for;
		return Multiplicities, count1;
	end if;
end function;

if g eq 2 then
	factors := [4] cat [p : p in [3..N] | IsPrime(p)];
else
	factors := [2] cat [p : p in [3..N] | IsPrime(p)];
end if;
CC := [* *];
for f in factors do
	CCModf := ComputeCorrectionCoefficients(g, f, q);
	Append(~CC, CCModf);
end for;


/*
Up to rescaling, Densities[i] is an approximation of nu(q,t+Floor(4*Sqrt(q))) with nu as in Equation (8). It is just an approximation for the following reasons:
-We are considering an approximation of the Sato-Tate density;
-We are approximating nu_ell by the value of Equation (6) for k=1 (in the case l=2, we take the value for k=2);
-We are taking the product of Equation (8) only for the primes up to N.
*/

Densities := [];
for t in [-Floor(2*g*Sqrt(q))..Floor(2*g*Sqrt(q))] do
	dens := STDensity(t/Sqrt(q));
	i := 0;
	for f in factors do
		i := i+1;
		tr := t mod f;
		dens := dens * CC[i][tr+1] / &+CC[i] ;
		dens:=RealField(10)!dens;
	end for;
	Append(~Densities, dens);
end for;

s2 := &+Densities;
Densities2 := [ dens/s2 : dens in Densities];

/*
Densities2[i] is an approximation of nu'(q,t+Floor(2*g*Sqrt(q))) with nu' as in Equation (9). 
*/
