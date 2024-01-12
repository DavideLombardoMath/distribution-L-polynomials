/*
In this code, we compute an approximation of the Sato-Tate density in genus g knowing the first N moments. To make this code work, one needs to define N and g.
*/
//g:=2;
//N:=100;
/*
Functions to compute the moments of the Sato-Tate density in arbitrary genus
*/
function InverseGammaOnIntegers(n)
	if n ge 1 then
		return 1 / Factorial(n-1);
	else
		return 0;
	end if;
end function;

function BesselI(n, z, order)
	return &+[ (1/2*z)^n * (1/4*z^2)^k / Factorial(k) * InverseGammaOnIntegers(n+k+1) : k in [0..order] ];
end function;

function F(m, z, order)
	return &+[ Binomial(m,j)*(BesselI(2*j-m, 2*z, order)-BesselI(2*j-m+2, 2*z, order)) : j in [0..order+m]];
end function;

/*
Generates a list of the first n moments
of the Sato-Tate distribution in genus g
*/
function SatoTateMoments(g, n)
	order := n+1;
	R<z> := PolynomialRing(Rationals());
	R<z> := FieldOfFractions(R);
	M := Matrix( R, g, g, [F(i+j-2, z, order) : i in [1..g], j in [1..g]] );
	A := Determinant(M);
	return [Evaluate(Derivative(A, k),0) : k in [0..n] ];
end function;

/*
Compute the Sato-Tate density in genus g from the knowledge of its moments. We compute the first N moments
*/
moments := SatoTateMoments(g, N);

function DilatedLegendrePolynomial(n, g)
	Ln := LegendrePolynomial(n);
	R := Parent(Ln);
	x := R.1;
	return Evaluate(Ln, x/(2*g));
end function;


function nthLegendreMoment(n, g)
	Ln := DilatedLegendrePolynomial(n, g);
	return (2*n+1)/(2*2*g)*&+[ Coefficients(Ln)[j] * moments[j] : j in [1..#Coefficients(Ln)]];
/*
We are rescaling by a factor (2*n+1)/(2*2*g) in order to get an orthonormal basis of L^2([-2g, 2g]).
*/
end function;

LegendreMoments := [nthLegendreMoment(n, g) : n in [0..N-1]];
approx := &+[ LegendreMoments[n] * DilatedLegendrePolynomial(n-1, g) : n in [1..#LegendreMoments] ];

/*
STDensity(x) computes an approximation of the Sato Tate density at x
When x is very close to 2g or -2g there are some problems with the convergence of the series. Since ST is very very small, we put it equal to 0.
*/
function STDensity(x)
	if x gt 1.75*g then
		return RealField(10)!0;
	end if;
	if x lt -1.75*g then
		return RealField(10)!0;
	end if;
	return Max(RealField(10)!0, RealField(10)!Evaluate(approx, x));
end function;