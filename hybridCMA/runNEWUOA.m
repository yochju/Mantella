function runNEWUOA(fname, dim, maxevals, ftarget, xinit, rhobegin, rhoend, npt)

%fname = 'frosenbrock';
%dim = 10;
%maxevals = 10000;
%ftarget = 1e-10;
%x = rand(1,dim);
%rhobegin = 10;
%rhoend = 1e-16;
%npt = 2*dim + 1;

[xbest, fbest, nevals] = newuoa(fname, dim, maxevals, ftarget, x, rhobegin, rhoend, npt);

