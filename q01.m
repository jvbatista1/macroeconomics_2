%%%%%%%%%%%%
% Macroeconomia II
% Lista I
% Questão 1
% João Victor Batista Lopes, CAEN/UFC
%%%%%%%%%%%%
% Applying it to a real function

f=@(x) exp(x) - exp(2.2087);
a=0;
b=4;
raiz = bisection(f,a,b);