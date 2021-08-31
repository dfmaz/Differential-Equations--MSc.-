%% DE oddity - Daniel Fern�ndez Mart�nez
% Se considera la ED x' = 10(x^a + t^2 + 0.2t - 1) con la CI x(0) = 1.
% Proponer y analizar un m�todo num�rico que �garantice�, en la medida
% de lo posible, que la soluci�n obtenida es fiable hasta t = 1.5 con un
% error de digamos 0.5 (o sea que |xnum(1.5) ? xex(1.5)| ? 0.5

% Caso 3: a = 1.1

% Par�metros auxiliares:
h = 0.001; % Paso temporal
t_ini = 0; % Extremo izquierdo del intervalo
t_fin = 1.5; % Extremo derecho del intervalo
t=t_ini:h:t_fin; % Intervalo
t_0=1; % Condici�n inicial (t)
x_0=1; % Condici�n inicial (x)
x=zeros(1,length(t)); % Vector x
x(t_0)=x_0; % Condici�n inicial

% ODE: 
a = 1.1; % a
F_tx=@(t,x) 10*(x^a - t^2 + 0.2*t - 1);

% M�todo de Runge - Kutta de 4� Orden
for i=1:(length(x)-1)
    
    p1 = F_tx(t(i),x(i));
    p2 = F_tx(t(i)+0.5*h,x(i)+0.5*h*p1);
    p3 = F_tx((t(i)+0.5*h),(x(i)+0.5*h*p2));
    p4 = F_tx((t(i)+h),(x(i)+p3*h));
    
    x(i+1) = x(i) + (1/6)*(p1+2*p2+2*p3+p4)*h;
end

% Visualizaci�n
figure;
plot(t,x);
title(['Soluci�n num�rica (a = ', num2str(a), ', h = ', num2str(h),')']);