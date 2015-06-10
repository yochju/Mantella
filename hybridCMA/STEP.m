function STEP
close all;


npt = 1001;
myfunc = 'fsphere'; a = -1.5;   b = 1;  feps = 1e-10;   
myfunc = 'fbrent';  a = -5;   	b = 5;  feps = 1e-10;  
myfunc = 'fmichalewicz';  a = -1;   	b = 2;  feps = 1e-10; 
myfunc = 'fmichalewicz2';  a = 0;   	b = 4;  feps = 1e-10; 

for i=1:npt
    x(i) = a + (b-a)*i/(npt-1);
    f(i) = feval(myfunc,x(i));
end;

plot(x,f);  hold on;




nt = 0;
nt = nt + 1;	t(nt) = a;   ft(nt) = feval(myfunc,t(nt)); plot(t(nt), ft(nt), '.', 'color', 'red');    hold on;
nt = nt + 1;    t(nt) = b;   ft(nt) = feval(myfunc,t(nt)); plot(t(nt), ft(nt), '.', 'color', 'red');    hold on;

niter = 150;
for iter=1:niter
    ftarget = max(ft);
    if (iter == 1)
        i1 = 1;
        i2 = 2;
    else
        D = zeros(1,1);
        for k=1:nt-1
            dx = t(k+1) - t(k);
            dy = ft(k+1) - ft(k);
            yhat = ftarget - ft(k) + feps;
            curD = (4*yhat - 2*dy + 4 * sqrt(yhat*yhat - yhat*dy))/(dx*dx);
            D(k) = curD;
        end;
        [Dmin indx] = sort(D);
        i1 = indx(1);
        i2 = indx(1)+1;
    end;
    nt = nt + 1;
    t(nt) = t(i1) + 0.5*(t(i2) - t(i1));
    ft(nt) = feval(myfunc,t(nt)); 
    [t indxt] = sort(t);
    ft = ft(indxt);
    plot(t, ft, '.', 'color', 'red');
    disp([num2str(iter + 2) ' ' num2str(max(ft))]);
end;

end

function f=fsphere(x)
  f=-sum(x.^2);
end

function f=fbrent(x)
  f=-(x - sin(x))*exp(-x*x);
end

function f=fmichalewicz(x)
  f=-(x*sin(10*x));
end
