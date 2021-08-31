% Tarea 2. Euler básico - Daniel Fernández Martínez
% Programar una rutina de Matlab (o lo que sea) para el método de Euler básico (hacia adelante).
% Modificarla para que dé el método de Euler implícito y el de Heun, también llamado "Euler mejorado".
% Aplicarlo a la ecuación logística u' = u(1-u), en 0<t<1 para u(0)=0.01 con diversos pasos.
% La salida debe poder comparar la solución exacta con las aproximaciones.

close all;
%% Función logística:
edo = @(u,t) (t^2-u) 

%% Parámetros de los métodos:
a = 0;      % Extremo izquierdo del intervalo
% b_exp = 1.1;      % Extremo derecho del intervalo
% b_imp = 1.03;
b_exp = 1.5;      % Extremo derecho del intervalo
b_imp = 1.5;
% b_exp = 1.175;      % Extremo derecho del intervalo
% b_imp = 1.08;
% b_exp = 1.2651;      % Extremo derecho del intervalo
% b_imp = 1.1172;
% b_exp = 1.2615;      % Extremo derecho del intervalo
ya = 1;  % Condición inicial

%% Comparación entre los métodos y la solución exacta en función de "i":
for i = 5:10:100
    figure;
    % Euler explícito:
    [t, y] = euler_explicito(edo, a, b_exp, ya, i);
    plot(t, y)
   % hold on;
    
%     % Euler implícito:
%     [t, y] = euler_implicito(edo, a, b_imp, ya, i);
%     plot(t, y,'r')
%     hold on;
%     hold off ;
    
    title(['Comparación de los métodos con h = ', num2str((b-a)/i)])
    legend('Implícito','Location','northwest')
end




