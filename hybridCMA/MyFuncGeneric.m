function Fit = MyFuncGeneric(x)
if (size(x,1) == 1)
    x = x';
end;
%if (fgeneric('fbest') - fgeneric('ftarget') > 0)
    Fit = feval('fgeneric',x);   
%else
%    Fit = fgeneric('fbest');
%end;
if (Fit*0 ~= 0)
    disp('NAN, MyFuncGeneric');
end;
%disp(x);
