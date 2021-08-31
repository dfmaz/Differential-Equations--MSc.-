%% Conducci�n del calor. Daniel Fern�ndez Mart�nez. Diciembre 2020.

%   CalorDF.m resuelve la siguiente Ecuaci�n en Derivadas Parciales por el
%   M�todo de Diferencias Finitas (M�todo Expl�cito y Crank-Nicholson):
%       u_t + c*u_x = a*u_xx + f,   0<=x<=pi, a=-2, c=1, f=0
%       Con condiciones de frontera: u(0,t) = u(pi,t) = 0
%       Y condici�n inicial: u(x,0)= u_0(x) = e^(-x)*(2*sin(x)-sin(2*x))
%       t = [0, 0.1]

%   Dicha ecuaci�n se denomina, seg�n varios autores, Ecuaci�n de
%   Conducci�n - Convecci�n en 1D, en r�gimen no estacionario, siendo el
%   par�metro 'c' el coeficiente de convecci�n y el par�metro 'a' el de
%   conducci�n o difusi�n.

% Ejemplo de ejecuci�n:
%   >> CalorDF
%   Introduce el n�mero de intervalos de tiempo: 10
%   Introduce el n�mero de divisiones del espacio: 10
%   Elecci�n del m�todo (0 = Expl�cito / 1 = Crank-Nicholson): 1
    
%   Salida: Gr�fica comparativa entre la aproximaci�n num�rica por
%   Diferencias Finitas con el n�mero de divisiones y m�todo especificados,
%   y la soluci�n exacta en t=0 y t=0.1.

clear all;
close all;

% Algunos par�metros generales
cond = 1; % Conducci�n
conv = -2; % Convecci�n
L = pi; % Longitud total
T = 0.1; % Tiempo final

% Par�metros necesarios para resolver la ecuaci�n por diferencias finitas
% Comunes al m�todo Expl�cito y Crank-Nicholson
Nt = input('Introduce el n�mero de intervalos de tiempo: '); % N�mero de intervalos de tiempo
Dt = T / Nt; % Paso temporal
Nx = input('Introduce el n�mero de divisiones del espacio: '); % N�mero de divisiones 
Dx = L/Nx; % Paso espacial
a = cond*Dt/(Dx*Dx); % Par�metro 'a' para la aproximaci�n de la derivada segunda (lambda en la descripci�n del m�todo)
c = conv*Dt/(Dx); % Par�metro 'c' para la aproximaci�n de la derivada primera (n�mero de Courant en la descripci�n del m�todo)
% Espec�ficos para Crank-Nicholson
alfa = c/2+a; % Par�metro 'alfa' para la diagonal secundaria inferior 
beta = c/2-a; % Par�metro 'beta' para la diagonal secundaria superior

% Condici�n inicial u_0(x)
for i = 1:Nx+1
    x(i) = (i-1)*Dx; % Definimos un vector x, para la discretizaci�n espacial
    u(i,1) = exp(-x(i)).*(2*sin(x(i))-sin(2*x(i)));
end

% Condiciones de contorno: temperatura inicial y final
for k=1:Nt+1
    u(1,k) = 0.;
    u(Nx+1,k) = 0;
    t(k) = (k-1)*Dt; % Definimos un vector t, para la discretizaci�n temporal
end

% Soluciones exactas en t = 0 y t = 0.1 para la visualizaci�n final
xsol = linspace(0,pi);
t = 0;
y0 = 2*exp(-2*t).*exp(-xsol).*sin(xsol)-exp(-5*t).*exp(-xsol).*sin(2*xsol);
t = 0.1;
yf = 2*exp(-2*t).*exp(-xsol).*sin(xsol)-exp(-5*t).*exp(-xsol).*sin(2*xsol);

% Elecci�n del m�todo de resoluci�n 
met = input('Elecci�n del m�todo (0 = Expl�cito / 1 = Crank-Nicholson): ');

% M�todo expl�cito
if met == 0
    tic;
    for k=1:Nt % Bucle temporal
        for i=2:Nx % Bucle espacial
            u(i,k+1) = u(i,k) - 0.5*c*(u(i+1,k)-u(i-1,k)) + a*(u(i-1, k)+u(i+1, k)-2.*u(i,k));
        end
    end
    toc;
    
    % Representaci�n gr�fica de la aproximaci�n a lo largo del tiempo
    for i = 1:Nt
        plot(x,u(:,i), '-o', xsol, y0, 'r--', xsol, yf, 'g--'),
        axis([0 pi 0 0.5]),
        title('Resoluci�n mediante el m�todo expl�cito'),
        xlabel('X'),
        ylabel('Temperatura (�C)'),
        grid on;
        pause(0.01);
    end
    hold on;
    legend('Expl�cito', 'Exacta: t=0', 'Exacta: t=0.1');
    hold off;
end

% M�todo de Crank-Nicholson en forma matricial
if met == 1
    tic;
    % Creaci�n de la matriz izquierda (rellenamos las diagonales
    % principales y secundarias)
    aa1(1:Nx-2)=-alfa;
    bb1(1:Nx-1)=2.+2.*a;
    cc1(1:Nx-2)=beta;
    MM1= diag(bb1,0)+diag(aa1,-1)+diag(cc1,1);
    [L,U] = lu(MM1); % Descomposic�n LU
    
    % Creaci�n de la matriz derecha (rellenamos las diagonales principales
    % y secundarias)
    aar(1:Nx-2)=alfa;
    bbr(1:Nx-1)=2.-2.*a;
    ccr(1:Nx-2)=-beta;
    MMr=diag(bbr,0)+diag(aar,-1)+diag(ccr,1);
    
    % Resoluci�n del sistema 
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
    
    % Representaci�n gr�fica de la aproximaci�n a lo largo del tiempo
    for i = 1:Nt
        plot(x,u(:,i), '-o', xsol, y0, 'r--', xsol, yf, 'g--')
        axis([0 pi 0 0.5]),
        title('Resoluci�n mediante el m�todo de Crank-Nicholson'),
        xlabel('X'),
        ylabel('Temperatura (�C)')
        grid on;
        pause(0.01);
    end
    hold on;
    legend('Crank-Nicholson', 'Exacta: t=0', 'Exacta: t=0.1');
    hold off;
end



