%% Euler implícito
% Invocación: [y t] = euler_implícito(f,a,b,ya,n)
%
% Input:
% f - función@
% a,b - intervalo
% ya - condición inicial
% n - número de subintervalos
%
% Output:
% y - solución
% t - variable tiempo

function [t y] = euler_implicito(f,a,b,ya,n)
 h = (b - a) / n;                                                
 t = [a zeros(1,n)]; 
 y = [ya zeros(1,n)];
 for i = 1 : n 
     t(i+1) = t(i) + h;
     yprima = y(i) + h * f(t(i),y(i));
     y(i+1) = y(i) + h * f(t(i+1),yprima);
 end