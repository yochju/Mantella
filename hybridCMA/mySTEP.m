function [nlocvals] = STEP(ia)


x_a = -5;
x_b = 5;
feps = 1e-10;
sol = x_a + (x_b-x_a)*rand(4,1);
ix = i-1;
init = 1;

nlocvals = 0;

newiter = 0;
if (init == 1)   
    newiter = 1;   
    init = 0;
end;

if (newiter == 1)
    ix = ix + 1;
    iter = 1;
    nt = 0;
    t = zeros(1,1);
    ft = zeros(1,1);

    nt = nt + 1;
    t(nt) = x_a;   
    sol(ix) = t(nt);
    ft(nt) = -sum(sol);

    nt = nt + 1;
    t(nt) = x_b;  
    sol(ix) = t(nt);
    ft(nt) = -sum(sol);

    newiter = 0;
    nlocvals = 2;
    return;
end;

        ftarget = max(ft);
        if (iter == 1)
            i1 = 1;
            i2 = 2;
        else
            D = zeros(1,1);
            if (1)
                k=1:nt-1;
                dx = t(k+1) - t(k);
                dy = ft(k+1) - ft(k);
                yhat = ftarget - ft(k) + feps;
                D = (4.*yhat - 2.*dy + 4 * sqrt(yhat.*yhat - yhat.*dy))./(dx.*dx);
            end;
            [Dmin indx] = sort(D);
            i1 = indx(1);
            i2 = indx(1)+1;
        end;
        nt = nt + 1;
        t(nt) = t(i1) + 0.5*(t(i2) - t(i1));
        sol(ix) = t(nt);    
        ft(nt) = -sum(sol); 
        [t indxt] = sort(t);
        ft = ft(indxt);
        iter = iter + 1;
        nlocvals = 1;
