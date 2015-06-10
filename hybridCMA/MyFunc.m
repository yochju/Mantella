function Fit = MyFunc(x)
global settings;
if (size(x,1) == 1)
    x = x';
end;
%if (fgeneric('fbest') - fgeneric('ftarget') > 0)
    Fit = feval('fgeneric',x) - settings.ftarget;    % does not change the results, but better for logs 
%else
%    Fit = fgeneric('fbest') - settings.ftarget;    
%end;

%Fit = feval('fgeneric',x) - settings.ftarget;    % does not change the results, but better for logs
%disp(x);
