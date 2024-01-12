/*
Printing
*/

if not ComputeActualNumbers then
	multTrue := [];
end if;

if ComputeActualNumbers then
	j := 0;
	"*** True densities ***";
	for t in [-Floor(2*g*Sqrt(q))..Floor(2*g*Sqrt(q))] do
		j := j+1;
		if multTrue[j] gt 10^(-5) then
			print "(", t, ",", RealField(5)!(multTrue[j]), ")";
		else
			print "(", t, ",", 0, ")";
		end if;
	end for;
end if;

j := 0;


"*** Predicted densities ***";
for t in [-Floor(2*g*Sqrt(q))..Floor(2*g*Sqrt(q))] do
	j := j+1;
	if Densities2[j] gt 10^(-5) then
		print "(", t, ",", RealField(5)!Densities2[j], ")";
	else
		print "(", t, ",", 0, ")";
	end if;
end for;

"*** Smooth interpolation using the Sato-Tate density ***";
SmoothDens := [STDensity(t / Sqrt(q)) : t in [-Floor(2*g*Sqrt(q))..Floor(2*g*Sqrt(q))] ];
sSmoothDens := &+SmoothDens;
SmoothDens := [ d / sSmoothDens : d in SmoothDens];

j := 0;
for t in [-Floor(2*g*Sqrt(q))..Floor(2*g*Sqrt(q))] do
	j := j+1;
	if SmoothDens[j] gt 10^(-5) then
		print "(", t, ",", RealField(5)!SmoothDens[j], ")";
	else
		print "(", t, ",", 0, ")";
	end if;
end for;
