/*
Implementation of the exact formula of Theorem 4.4
*/

function T(m, q, zeta, eta)
	if m eq 0 then
		if eta eq 0 then
			t := 1;
		else
			t := 0;
		end if;
		return q*t - 1;
	end if;
	F := GF(q);
	X := [F!1..F!(q-1)];
	Y := CartesianPower(X, m-1);
	res := 0;
	for y in Y do
		if m gt 1 then
			c := &+[ alpha : alpha in y] + &+[ F!zeta/alpha : alpha in y];
		else
			c := 0;
		end if;
		// we count the number of solutions of c + x + z/x = eta,
		// that is x^2 + (c-eta) * x + z = 0
		R<t> := PolynomialRing(GF(q));
		res := res + #Roots(t^2 + (c-eta)*t + zeta);
		if q ne 2 then
			assert #Roots(t^2 + (c-eta)*t + zeta) eq (1 + LegendreSymbol( Integers()!((c-eta)^2 -4*zeta), q));
		end if;
	end for;
	return q*res - (q-1)^m;
end function;

function qBinomial(n, r, q)
	if r eq 0 then
		return 1;
	else
		return &*[ (q^(n-j)-1)/(q^(r-j)-1) : j in [0..(r-1)] ];
	end if;
end function;

function R(q, m, l)
	if m le 2*l then
		return 0;
	end if;

	if l eq 0 then
		return 1;
	end if;
	if l eq 1 then
		res := 0;
		for j1 in [1..(m-l-1)] do
			res := res + q^(m-1-j1)-1;
		end for;
		return res;
	end if;

	if l eq 2 then
		res := 0;
		for j2 in [1..(m-l-1)] do
			for j1 in [1..(j1-1)] do
				res := res + ( q^(m-1-j1)-1 )*( q^(m-2-j2)-1 );
			end for;
		end for;
	end if;

	error "Not implemented";
end function;

function E(n, q, zeta, eta)
	res := 0;
	for b in [0..Floor(n/2)] do
		innerSum := &+[ q^l*R(q, n-2*b+1,l)*T(n-2*b-2*l, q, zeta, eta) : l in [0..Floor(n/2-b)] ];
		innerProd := q^(b^2+b) * qBinomial(n, 2*b, q);
		if b ge 1 then
			innerProd := innerProd * &*[ q^(2*j-1)-1 : j in [1..b] ];
		end if;
		res := res + q^(n^2-1) * innerProd * innerSum;
	end for;
	return res;
end function;

function C(n, q, zeta, eta)
	return q^(n^2-1) * &*[ q^(2*j)-1 : j in [1..n] ] + E(n, q, zeta, eta);
end function;
	
//By Theorem 4.3, C(n, q, zeta, eta) counts the matrices in Gsp_{2n}(F_q) with multiplier zeta and trace eta.

/*
End of theorem 4.3
*/
