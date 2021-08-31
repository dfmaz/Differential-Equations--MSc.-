%% Conducción del calor. Daniel Fernández Martínez. Diciembre 2020.

%   CalorDF.m resuelve la siguiente Ecuación en Derivadas Parciales por el
%   Método de Diferencias Finitas (Método Explícito y Crank-Nicholson):
%       u_t + c*u_x = a*u_xx + f,   0<=x<=pi, a=-2, c=1, f=0
%       Con condiciones de frontera: u(0,t) = u(pi,t) = 0
%       Y condición inicial: u(x,0)= u_0(x) = e^(-x)*(2*sin(x)-sin(2*x))
%       t = [0, 0.1]

%   Dicha ecuación se denomina, según varios autores, Ecuación de
%   Conducción - Convección en 1D, en régimen no estacionario, siendo el
%   parámetro 'c' el coeficiente de convección y el parámetro 'a' el de
%   conducción o difusión.

% Ejemplo de ejecución:
%   >> CalorDF
%   Introduce el número de intervalos de tiempo: 10
%   Introduce el número de divisiones del espacio: 10
%   Elección del método (0 = Explícito / 1 = Crank-Nicholson): 1
    
%   Salida: Gráfica comparativa entre la aproximación numérica por
%   Diferencias Finitas con el número de divisiones y método especificados,
%   y la solución exacta en t=0 y t=0.1.

clear all;
close all;

% Algunos parámetros generales
cond = 1; % Conducción
conv = -2; % Convección
L = pi; % Longitud total
T = 0.1; % Tiempo final

% Parámetros necesarios para resolver la ecuación por diferencias finitas
% Comunes al método Explícito y Crank-Nicholson
Nt = input('Introduce el número de intervalos de tiempo: '); % Número de intervalos de tiempo
Dt = T / Nt; % Paso temporal
Nx = input('Introduce el número de divisiones del espacio: '); % Número de divisiones 
Dx = L/Nx; % Paso espacial
a = cond*Dt/(Dx*Dx); % Parámetro 'a' para la aproximación de la derivada segunda (lambda en la descripción del método)
c = conv*Dt/(Dx); % Parámetro 'c' para la aproximación de la derivada primera (número de Courant en la descripción del método)
% Específicos para Crank-Nicholson
alfa = c/2+a; % Parámetro 'alfa' para la diagonal secundaria inferior 
beta = c/2-a; % Parámetro 'beta' para la diagonal secundaria superior

% Condición inicial u_0(x)
for i = 1:Nx+1
    x(i) = (i-1)*Dx; % Definimos un vector x, para la discretización espacial
    u(i,1) = exp(-x(i)).*(2*sin(x(i))-sin(2*x(i)));
end

% Condiciones de contorno: temperatura inicial y final
for k=1:Nt+1
    u(1,k) = 0.;
    u(Nx+1,k) = 0;
    t(k) = (k-1)*Dt; % Definimos un vector t, para la discretización temporal
end

% Soluciones exactas en t = 0 y t = 0.1 para la visualización final
xsol = linspace(0,pi);
t = 0;
y0 = 2*exp(-2*t).*exp(-xsol).*sin(xsol)-exp(-5*t).*exp(-xsol).*sin(2*xsol);
t = 0.1;
yf = 2*exp(-2*t).*exp(-xsol).*sin(xsol)-exp(-5*t).*exp(-xsol).*sin(2*xsol);

% Elección del método de resolución 
met = input('Elección del método (0 = Explícito / 1 = Crank-Nicholson): ');

% Método explícito
if met == 0
    tic;
    for k=1:Nt % Bucle temporal
        for i=2:Nx % Bucle espacial
            u(i,k+1) = u(i,k) - 0.5*c*(u(i+1,k)-u(i-1,k)) + a*(u(i-1, k)+u(i+1, k)-2.*u(i,k));
        end
    end
    toc;
    
    % Representación gráfica de la aproximación a lo largo del tiempo
    for i = 1:Nt
        plot(x,u(:,i), '-o', xsol, y0, 'r--', xsol, yf, 'g--'),
        axis([0 pi 0 0.5]),
        title('Resolución mediante el método explícito'),
        xlabel('X'),
        ylabel('Temperatura (ºC)'),
        grid on;
        pause(0.01);
    end
    hold on;
    legend('Explícito', 'Exacta: t=0', 'Exacta: t=0.1');
    hold off;
end

% Método de Crank-Nicholson en forma matricial
if met == 1
    tic;
    % Creación de la matriz izquierda (rellenamos las diagonales
    % principales y secundarias)
    aa1(1:Nx-2)=-alfa;
    bb1(1:Nx-1)=2.+2.*a;
    cc1(1:Nx-2)=beta;
    MM1= diag(bb1,0)+diag(aa1,-1)+diag(cc1,1);
    [L,U] = lu(MM1); % Descomposicón LU
    
    % Creación de la matriz derecha (rellenamos las diagonales principales
    % y secundarias)
    aar(1:Nx-2)=alfa;
    bbr(1:Nx-1)=2.-2.*a;
    ccr(1:Nx-2)=-beta;
    MMr=diag(bbr,0)+diag(aar,-1)+diag(ccr,1);
    
    % Resolución del sistema 
    for k=2:Nt % Bucle temporal
        uu=u(2:Nx,k-1);
        qq = MMr*uu;
        w(1)=qq(1);
        for j=2:Nx-1
            w(j)=qq(j)-L(j,j-1)*w(j-1);
        end
        u(Nx,k)=w(Nx-1)/U(Nx-1,Nx-1);
        for i=Nx-1:-1:2
            u(i,k)=(w(i-1)-U(i-1,i)*u(i+1,k))/U(i-1,i-1);
        end
    end
    toc;
    
    % Representación gráfica de la aproximación a lo largo del tiempo
    for i = 1:Nt
        plot(x,u(:,i), '-o', xsol, y0, 'r--', xsol, yf, 'g--')
        axis([0 pi 0 0.5]),
        title('Resolución mediante el método de Crank-Nicholson'),
        xlabel('X'),
        ylabel('Temperatura (ºC)')
        grid on;
        pause(0.01);
    end
    hold on;
    legend('Crank-Nicholson', 'Exacta: t=0', 'Exacta: t=0.1');
    hold off;
end



