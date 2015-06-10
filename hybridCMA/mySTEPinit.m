function mySTEPinit(N,myfunc,x_a,x_b, perc, niter, precision, feps)

global STEPalgo;
STEPalgo.curbest = 1e+30;
STEPalgo.prevbest = STEPalgo.curbest;
STEPalgo.perc = perc;
STEPalgo.evals = 0;
STEPalgo.indx = zeros(1,N);
STEPalgo.precision = precision;
STEPalgo.iter = 0;
STEPalgo.niter = niter;

global STEPalgos;
for i=1:N
    STEPalgos(i).x_a = x_a;  STEPalgos(i).x_b = x_b;
    STEPalgos(i).myfunc = myfunc;
    STEPalgos(i).feps = feps;
    STEPalgos(i).sol = x_a + (x_b-x_a)*rand(N,1);
    STEPalgos(i).ix = i-1;
    STEPalgos(i).init = 1;
end;