%% Conducci�n del calor. Daniel Fern�ndez Mart�nez. Diciembre 2020.

%   CalorEF.m resuelve la siguiente Ecuaci�n en Derivadas Parciales por el
%   M�todo de Elementos Finitos:
%       u_t + c*u_x = a*u_xx + f,   0<=x<=pi, a=-2, c=1, f=0
%       Con condiciones de frontera: u(0,t) = u(pi,t) = 0
%       Y condici�n inicial: u(x,0)= u_0(x) = e^(-x)*(2*sin(x)-sin(2*x))
%       t = [0, 0.1]

%   Dicha ecuaci�n se denomina, seg�n varios autores, Ecuaci�n de
%   Conducci�n - Convecci�n en 1D, en r�gimen no estacionario, siendo el
%   par�metro 'c' el coeficiente de convecci�n y el par�metro 'a' el de
%   conducci�n o difusi�n.

% Ejemplo de ejecuci�n:
%   >> CalorEF
%   Introduce el n�mero de divisiones del espacio: 10
%   Elecci�n del m�todo (0 = Euler hacia adelante/ 1 = Euler hacia atr�s
%   / 2 = Trapezoidal): 1

%   Salida: Gr�fica comparativa entre la aproximaci�n num�rica por
%   Elementos finitos con el n�mero de divisiones y m�todo especificados,
%   y la soluci�n exacta en t=0 y t=0.1.

clear all;
close all;

% Discretizaci�n espacial
div = input('Introduce el n�mero de divisiones del espacio: ');
h = pi/div; % Tama�o de las suvdivisiones
x = 0:h:pi;

% Discretizaci�n temporal
L = 10; % N�mero de intervalos de tiempo
T = 0.1; % Tiempo final
t = linspace(0,T,L+1); 

% Condici�n inicial u_0(x)
U = exp(-x).*(2*sin(x)-sin(2*x));
U = U';

% Condiciones de contorno: temperatura inicial y final
g = [0 0]; 

% Par�metros Robin BC
kappa = [10 10]; 

% Elementos del sistema a resolver
M = MassMat1D(x); % Matriz de masas
A1 = StiffMat1D_difusion(x,@a,kappa); % Matriz de rigidez: conducci�n
A2 = StiffMat1D_conveccion(x,@c); % Matriz de rigidez: convecci�n
A = A1 + A2; % Suma de las dos matrices de rigidez
b = LoadVec1D(x,@f,kappa,g); % Vector de carga

% Soluciones exactas en t = 0 y t = 0.1 para la visualizaci�n final
xsol = linspace(0,pi);
t0 = 0;
y0 = 2*exp(-2*t0).*exp(-xsol).*sin(xsol)-exp(-5*t0).*exp(-xsol).*sin(2*xsol);
tf = 0.1;
yf = 2*exp(-2*tf).*exp(-xsol).*sin(xsol)-exp(-5*tf).*exp(-xsol).*sin(2*xsol); 

% Elecci�n del m�todo de resoluci�n     
met = input('Elecci�n del m�todo (0 = Euler hacia adelante / 1 = Euler hacia atr�s / 2 = Trapezoidal): ');

% M�todo de Euler hacia adelante (expl�cito)
if met == 0 
    tic;
    for l = 1:L % Bucle de tiempo
        k = t(l+1)-t(l); %  Paso temporal
        U = M\(M-k*A)*U+k*b; % Resoluci�n del sistema
        % Visualizaci�n gr�fica de la aproximaci�n a lo largo del tiempo
        plot(x, U, '-o', xsol, y0, 'r--', xsol, yf, 'g--'), 
        axis([0 pi 0 0.5]), 
        title('Resoluci�n mediante EF (Euler hacia adelante)'),
        xlabel('X')
        ylabel('Temperatura (�C)')
        grid on,
        pause(0.1);
    end
    toc;
    hold on;
    legend('EF Euler expl�cito', 'Exacta')
    hold off;
end

% M�todo de Euler hacia atr�s (impl�cito)
if met == 1
    tic;
    for l = 1:L % Bucle de tiempo
        k = t(l+1) - t(l); % Paso temporal
        U = (M + k*A)\(M*U + k*b); % Resoluci�n del sistema
        % Visualizaci�n gr�fica de la aproximaci�n a lo largo del tiempo
        plot(x, U, '-o', xsol, y0, 'r--', xsol, yf, 'g--'), 
        axis([0 pi 0 0.5]), 
        title('Resoluci�n mediante EF (Euler hacia atr�s)'),
        xlabel('X')
        ylabel('Temperatura (�C)')
        grid on,
        pause(0.1);
    end
    toc;
    hold on;
    legend('EF Euler impl�cito', 'Exacta')
    hold off;
end

% M�todo Trapezoidal
if met == 2
    tic;
    for l = 1:L % Bucle de tiempo
        k = t(l+1) - t(l); % Paso temporal
        U = (M + k/2*A)\((M-k/2*A)*U + k/2*2*b); % Resoluci�n del sistema
        % Visualizaci�n gr�fica de la aproximaci�n a lo largo del tiempo
        plot(x, U, '-o', xsol, y0, 'r--', xsol, yf, 'g--'), 
        axis([0 pi 0 0.5]), 
        title('Resoluci�n mediante EF (Trapezoidal)'),
        xlabel('X')
        ylabel('Temperatura (�C)')
        grid on,
        pause(0.1);
    end
    toc;
    hold on;
    legend('EF Trapezoidal', 'Exacta')
    hold off;
end

% Coeficiente de conducci�n o difusi�n
function y = a(x)
y = 1; 
end

% Coeficiente de convecci�n
function y = c(x)
y = -2; 
end

% Funci�n f
function y = f(x)
y = 0; 
end

% Creaci�n de la matriz de masas global
function M = MassMat1D(x)
n = length(x)-1; % N�mero de elementos del espacio
M = zeros(n+1,n+1); % Guardamos espacio de memoria
for i = 1:n % Bucle para cada elemento
    h = x(i+1) - x(i); % Distancia entre elementos
    % Rellenamos la matriz
    M(i,i) = M(i,i) + h/3;
    M(i,i+1) = M(i,i+1) + h/6;
    M(i+1,i) = M(i+1,i) + h/6;
    M(i+1,i+1) = M(i+1,i+1) + h/3;
end
end

% Creaci�n de la matriz de rigidez global: difusi�n
function A = StiffMat1D_difusion(x,a,kappa)
n = length(x)-1;
A = zeros(n+1,n+1);
for i = 1:n
    h = x(i+1) - x(i);
    xmid = (x(i+1) + x(i))/2; % Punto medio entre elementos
    amid = a(xmid); % Evaluamos la funci�n en dicho punto
    % Rellenamos la matriz
    A(i,i) = A(i,i) + amid/h; 
    A(i,i+1) = A(i,i+1) - amid/h;
    A(i+1,i) = A(i+1,i) - amid/h;
    A(i+1,i+1) = A(i+1,i+1) + amid/h;
end
A(1,1) = A(1,1) + kappa(1);
A(n+1,n+1) = A(n+1,n+1) + kappa(2);
end

% Creaci�n de la matriz de rigidez global: convecci�n
function A = StiffMat1D_conveccion(x,c)
n = length(x)-1;
A = zeros(n+1,n+1);
for i = 1:n
    xmid = (x(i+1) + x(i))/2; 
    cmid = c(xmid); 
    A(i,i) = A(i,i) - cmid/2; 
    A(i,i+1) = A(i,i+1) + cmid/2;
    A(i+1,i) = A(i+1,i) - cmid/2;
    A(i+1,i+1) = A(i+1,i+1) + cmid/2;
end
end

% Creaci�n del vector de carga global
function b = LoadVec1D(x,f,kappa,g)
n = length(x)-1;
b = zeros(n+1,1);
for i = 1:n
    % Rellenamos el vector
    h = x(i+1) - x(i);
    b(i) = b(i) + f(x(i))*h/2;
    b(i+1) = b(i+1) + f(x(i+1))*h/2;
end
% A�adimos las condiciones iniciales
b(1) = b(1) + kappa(1)*g(1);
b(n+1) = b(n+1) + kappa(2)*g(2);
end