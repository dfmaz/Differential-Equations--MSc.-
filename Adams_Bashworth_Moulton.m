%% M�todos multipaso. Daniel Fern�ndez Mart�nez.
% Componer una rutina con alg�n m�todo multipaso. Candidatos razonables
% son el "salto de rana", Adams-Bashforth de 2� orden o de 4� orden, 
% Adams-Moulton y el Milne. Todos ellos pueden verse en las pp. 178-181 del
% libro de V�zquez.
% 
% La rutina debe aplicarse al PVI u'=u^2-2/t^2 con las CI u(1)=-1/2 y
% comparar con la soluci�n exacta u(t)=1/t-3t^2/(1+t^3).
% 
% Usar un m�todo de un paso adecuado para calcular los primeros valores.

%% Ecuaci�n a resolver:
edo = @(u,t) (u.^2-2/t.^2) 

%% Soluci�n de la EDO:
min = 1;
max = 2;
tspan = [min, max];                  % Intervalo
u0 = -1/2;                           % Condicion inicial
n = 100;                             % Paso temporal h = (max - min)/n

[u,t] = am(edo,min,max,u0,n);        % Soluci�n con Adams-Bashworth-Moulton
plot(t,u)
hold on;
t = linspace(min, max);
u = 1./t-3.*t.^2./(1+t.^3);          % Soluci�n exacta
plot(t,u)
hold off;
grid on;
title(['Adams-Barthworth-Moulton, h = ', num2str((max-min)/n)]);
legend('Adams-Barthworth-Moulton', 'exacta');

%% Adams-Bashforth-Moulton predictor - corrector
% 
% Runge-Kutta de 4� orden para los primeros pasos
% Adams-Bashforth de 4� orden como precdictor,
% Adams-Moulton como corrector

function [u t] = am(f,a,b,u0,n)

h = (b - a) / n;    % Paso temporal
u(1,:) = u0;        % Condici�n inicial (u)
t(1) = a;           % Condici�n inicial (t)

m = min(3,n);       % Condici�n de parada del for del Runge-Kutta

for i = 1 : m       % Runge - Kutta de 4� orden para los primeros pasos
    t(i+1) = t(i) + h;
    p(i,:) = f(t(i), u(i,:));
    p2 = f(t(i) + h / 2, u(i,:) + p(i,:) * h /2);
    p3 = f(t(i) + h / 2, u(i,:) + p2 * h /2);
    p4 = f(t(i+1), u(i,:) + p3 * h);
    u(i+1,:) = u(i,:) + (p(i,:) + p2 + p2 + p3 + p3 + p4) * h / 6;
end

for i = m + 1 : n   % Sucesivos pasos con m�todo predictor - corrector
    p(i,:) = f(t(i), u(i,:));
    u(i+1,:) = u(i,:) + h/24 * (55 * p(i,:) - 59 * p(i-1,:) + 37 * p(i-2,:) - 9 * p(i-3,:)); % Predictor (Adams-Bashforth)
    t(i+1) = t(i) + h;
    u(i+1,:) = u(i,:) + h/24 * (9 * f(t(i+1), u(i+1,:)) + 19 * p(i,:) - 5 * p(i-1,:) + p(i-2,:)); % Corrector (Adams-Moulton)
end
end