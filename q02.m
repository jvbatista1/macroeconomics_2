%%%%%%%%%%%%
% Macroeconomia II
% Lista I
% Questão 1
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%
%%
% Declaring the function for the bisection method
function c = bisection(f,a,b)

if f(a)*f(b)>0
    disp("wrong choice for a and b")
else
    c = (a+b)/2;
    err = abs(f(c));
    while err > 1e-7
        if f(a)*f(b)<0
            b=c;
        else
            a=c;
        end
        c = (a+b)/2;
        err = abs(f(c));
    end
end
%%
% Applying it to a real function

f=@(x) exp(x) - exp(2.2087);
a=0;
b=4;
raiz = bisection(f,a,b)