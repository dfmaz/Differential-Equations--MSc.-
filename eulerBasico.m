% Tarea 2. Euler b�sico - Daniel Fern�ndez Mart�nez
% Programar una rutina de Matlab (o lo que sea) para el m�todo de Euler b�sico (hacia adelante).
% Modificarla para que d� el m�todo de Euler impl�cito y el de Heun, tambi�n llamado "Euler mejorado".
% Aplicarlo a la ecuaci�n log�stica u' = u(1-u), en 0<t<1 para u(0)=0.01 con diversos pasos.
% La salida debe poder comparar la soluci�n exacta con las aproximaciones.

close all;
%% Funci�n log�stica:
edo = @(u,t) (t^2-u) 

%% Par�metros de los m�todos:
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
ya = 1;  % Condici�n inicial

%% Comparaci�n entre los m�todos y la soluci�n exacta en funci�n de "i":
for i = 5:10:100
    figure;
    % Euler expl�cito:
    [t, y] = euler_explicito(edo, a, b_exp, ya, i);
    plot(t, y)
   % hold on;
    
%     % Euler impl�cito:
%     [t, y] = euler_implicito(edo, a, b_imp, ya, i);
%     plot(t, y,'r')
%     hold on;
%     hold off ;
    
    title(['Comparaci�n de los m�todos con h = ', num2str((b-a)/i)])
    legend('Impl�cito','Location','northwest')
end




