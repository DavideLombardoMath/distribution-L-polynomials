/*
We check the probability that a random matrix in GSp_6(F_2) has even trace.
Note that GSp_6(F_2)=Sp_6(F_2) since F_2* consists of a single element.
*/
Sp6 := Sp(6, GF(2));
t0 := #{M : M in Sp6 | Trace(M) eq 0 };
t1 := #{M : M in Sp6 | Trace(M) eq 1 };
assert t0 + t1 eq #Sp6;
t0;
#Sp6;
t0/(t0+t1);