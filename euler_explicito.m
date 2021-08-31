%% Euler expl�cito
% Invocaci�n: [y t] = euler_explicito(f,a,b,ya,n)
%
% Inputs:
% f - funci�n @
% a,b - intervalo
% ya - condicion inicial
% n - n�mero de subintervalos
%
% Outputs:
% y - soluci�n
% t - variable tiempo

function [t y] = euler_explicito(f,a,b,ya,n)
h = (b - a) / n;
y(1) = ya;
t(1) = a;
for i = 1 : n
    y(i+1) = y(i) + h * f(t(i),y(i));
    t(i+1) = t(i) + h;
end