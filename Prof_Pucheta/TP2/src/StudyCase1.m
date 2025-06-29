
clear all; clc; close all;
%% Sistema de Control II
%% Actividad Práctica 2
% Prof:Pucheta Julian
% Alumna: Gonzalez Macarena V.
%% Sistema de tres variables de estado con parámetros calculados
%% Importación de datos con #Import Data y grráficos obtenidos 
opts = spreadsheetImportOptions("NumVariables", 5);
% Specify sheet and range
opts.Sheet = "Hoja1";
opts.DataRange = "A2:E1501";
% Specify column names and types
opts.VariableNames = ["TiempoSeg", "VelocidadAngularradseg", "CorrienteEnArmaduraA", "TensinV", "Torque"];
opts.VariableTypes = ["double", "double", "double", "double", "double"];

% Import the data
tbl = readtable("C:\Users\Maca\Documents\UNC\4toAÑO\1CUATRI\SC2\Folder Github TP´s\Control-System-II-TP-s\Prof_Pucheta\TP2\documentation\Curvas_Medidas_Motor_2025_v.xlsx", opts, "UseExcel", false);

% Convert to output type
time = tbl.TiempoSeg;
VelocidadAngularradseg = tbl.VelocidadAngularradseg;
CteEnA = tbl.CorrienteEnArmaduraA;
Va = tbl.TensinV;
Torque = tbl.Torque;
%Gráficos
figure(1); 
hold on

subplot(4,1,1);
plot(time , VelocidadAngularradseg, 'g', 'LineWidth', 2);
title('Velocidad angular');
grid on;

subplot(4,1,2);
plot(time , CteEnA, 'b', 'LineWidth', 2);
title('Corriente de armadura');
grid on;

subplot(4,1,3);
plot(time, Va, 'r', 'LineWidth', 2);
title('Voltaje de entrada'); 
grid on;

subplot(4,1,4);
plot(time,Torque, 'm', 'LineWidth', 2);
title('Torque'); 
grid on;
hold off
% si usaramos la asignación de polos el tiempo mínimo de integración sería:  
Max_time= max(time)
Ct=numel(time)
ti= Max_time/(Ct-1)
dif= diff(time);


%% Caso de estudio 1.Item1 _inciso (a)
%Datos obtenidos en TP1
Ra = 2.258930051299405;
Laa = 0.005026901184834716;
Ki = 0.25965987053759737;
Jm = 0.0028472626983113334;
Bm = 0.0014165170369840668;
Km = 0.2500481104997174;

%Modelado de motor en espacios de estados
  %Variables de estados elegidas 
%x1=ia
%x2=wr
%x3=tita 
     %Matrices
%xp=Ax+Bu 
A = [-Ra/Laa -Km/Laa 0 ; Ki/Jm -Bm/Jm 0 ; 0 1 0];  % matriz de estados
B = [1/Laa 0 ;
    0 -1/Jm ;
    0 0];               % matriz de entrada 
C = [0 0 1];                               % matriz de salida                                                       
D = [0 0];                                % matriz de transmisión directa

% ecuación de salida 
%y= C*x


%% De espacio de estados a función de transferencia
sys = ss(A,B,C,D); % convierte el espacio de estado en ft contemplando si tiene dos entradas o no 
FTita = tf(sys);   %El sistema es MIMO entonces se muestran las dos FT´s 

%% Referencia
% El sistema cuenta con una referencia de pi/2 a -pi/2 
%ref = (pi/2)*square(2*pi*(1/10)*t);

%Al agrega una ref cambia la Ley de control y para tener una mejora en el
%rechazo a las perturbaciones se coloca un integrador debido a esto se
%amplian las matrices A y B 
Aamp = [A zeros(3,1); -C 0]
Bamp = [B(:,1); 0]
Camp = [C 0];

%% Verificación de condición de controlabilidad 
%Para saber si a mi sistema le puedo aplicar un controlador debo verificar
%que cumpla la condición de controlabilidad :el rango de M debe ser = a la
%dim de el espacio de estados ampliado =4 
Mamp = [Bamp  Aamp*Bamp  Aamp^2*Bamp Aamp^3*Bamp];
Rango_Mamp = rank(Mamp); 
% Rango_Mamp = 4 -> verifica entonces se puede aplicar un controlador 

%% Aplicación de controlador LQR 
% Se necesitan las matrices que forman la función de costo Q y R ;A y B
% ampliadas para utilizar el comando LQR que devuelve el K:matriz de
% ganancias optimas; S la solución de Riccati y P los polos del sist en LC

Q = diag([0.5 0.5 10 2300]) % matriz de ponderación para los estados 
R = 0.1; % matriz de ponderación para la entrada de control 
 
[K, S, P] =  lqr(Aamp,Bamp,Q,R)
%% Implementación de sistema a LC -Método de Euler  

% Para calcular el tiempo de integración para Euler se busca  la dinámica
% más rapida con los polos de lazo cerrado. Estos se puede obtener  como: 
% poles = eig(Aamp-Bamp(:,1)*K) o usando el resultado del comando lqr 
% El polo más rápido es:
lambda = max(P)

% Entonces, se calcula $t_r$, de modo de simular con un tiempo menor al mismo.
tr = log(0.95)/lambda

%Como se sugiere un entre ti =tr/3 -ti= tr/30 se opta:
ti=tr/(5)

Tsim = 10; %definido en cosigna

t = 0:ti:(Tsim-ti);

% Se declara aquí la referencia para hacer uso de t 
ref = (pi/2)*square(2*pi*(1/10)*t);

% El torque de entrada del sistema es: (mismo que tabla)
TL = zeros(1, length(t));
for ii=1:Tsim/ti
    varr = t(ii);
    if (varr>=0.7 && varr<=1.5)
        TL(ii) = 0.12;
    else
        TL(ii) = 0;
    end
end

% Las condiciones iniciales del sistema son nulas:
ia(1) = 0;
theta(1) = 0;
omega(1) = 0;
% El vector de estados del sistema está dado por:
stateVec = [ia(1) omega(1) theta(1)]';
%Punto de operación
xop = [0 0 0]';
x = [ia(1) omega(1) theta(1)];
zeta(1) = 0;
integ(1) = zeta(1);

%% Simulación: 
for i = 1:(Tsim/ti)
    zetaP = ref(i)-C*stateVec;
    zeta(i) = integ+zetaP*ti;
    u(i) = -K(1:3)*stateVec-K(4)*zeta(i);
    ia(i) = x(1);
    omega(i) = x(2);
    theta(i) = x(3);
    x1P = -Ra*x(1)/Laa-Km*x(2)/Laa+u(i)/Laa;
    x2P = Ki*x(1)/Jm-Bm*x(2)/Jm-TL(i)/Jm;
    x3P = x(2);
    xP = [x1P x2P x3P]';
    x = x+ti*xP;
    stateVec = [ia(i) omega(i) theta(i)]';
    integ = zeta(i);
end

%% GRAFICOS
%ANGULO
figure(2);
plot(t, ref, 'LineWidth', 1.5, 'DisplayName', 'Referencia');
hold on;
plot(t, theta, 'LineWidth', 1.5, 'DisplayName', 'Salida');
xlabel('Tiempo [seg]');
ylabel('Ángulo [rad]');
title('Comparación de Señal de referencia y Salida');
grid on;
legend;
% corriente de armadura y la acción de control.
figure(3)
plot(t,ia,'LineWidth',1.5)
xlabel('Tiempo [seg]')
ylabel('Corriente [A]')
title('Salida i_a')
grid
figure(4)
plot(t,u,'LineWidth',1.5)
xlabel('Tiempo [seg]')
title('Acción de control')
grid

figure(5)
plot(t,TL,'LineWidth',1.5)
xlabel('Tiempo [seg]')
ylabel('Torque [Nm]')
title('Torque de entrada T_L')
grid


 %% Caso de estudio 1.Item1 _inciso (b)
 % b-Asumiendo que no puede medirse directamente la corriente, pero sí la velocidad y el ángulo, 
 %proponer un controlador que logre el objetivo. 
%% Se aplica observador 
%Diseñar la ganancia del observador K_o para el sistema original
%S1 es análogo a diseñar una ganancia de controlador para el
% Sistema Dual S2
C_i = [0 0 1; 1 0 0]
Ao = A' 
Bo = C_i' 
Co = B' 
Qo = diag([10 1000000 0.01]) %covarianzas del ruido de proceso y de medición
Ro = diag([0.005,1000])
[Ko ,So, Po] = lqr(Ao,Bo,Qo,Ro)
obsStateVector = [ia(1) omega(1) theta(1)]';
xObs = [0 0 0]'; 
%%Cálculo del tiempo con Controlador + Observador 
lambdaOb = max(Po)
%------  USANDO EL POLO DEL 
% Entonces, se calcula $t_r$, de modo de simular con un tiempo menor al mismo.
tr_ob = log(0.95)/lambdaOb

%Como se sugiere un entre ti =tr/3 -ti= tr/30 se opta:
ti_ob=tr_ob/(5)

Tsim = 10; %definido en cosigna

t_ob = 0:ti_ob:(Tsim-ti_ob);
% La simulación:

for i = 1:(Tsim/ti_ob)
    zetaP = ref(i)-Camp(1:3)*obsStateVector-Camp(4)*integ;
    zeta(i) = integ+zetaP*ti_ob;
    u(i) = -K(1:3)*obsStateVector-K(4)*zeta(i);
    iaO(i)= xObs(1);
    omegaO(i)= xObs(2);
    thetaO(i)= xObs(3);
    yO(i) = C*obsStateVector;
    y(i) = Camp(1:3)*obsStateVector+Camp(4)*integ;
    xTP = A*xObs+B*u(i)+Ko'*(y(:,i)-yO(:,i));
    xObs = xObs+xTP*ti_ob;
    integ = zeta(i);
    obsStateVector =[iaO(i) omegaO(i) thetaO(i)]';
end
%% 
% La corriente observada es:
% figure(6)
% plot(t,iaO,'LineWidth',1.5)
% xlabel('Tiempo [seg]')
% ylabel('Corriente [A]')
% title('Corriente de armadura i_a observada')
% grid
% %% 
% Ahora, para comparar con la real obtenida anteriormente:
figure(7)
plot(t_ob,iaO,'LineWidth',1.5)
hold on
plot(t,ia,'LineWidth',1.5)
xlabel('Tiempo [seg]')
ylabel('Corriente [A]')
title('Corriente de armadura i_a')
legend('Observada', 'Real')
grid on
hold off
%Comparación de Theta con controlador y theta con controlador + observador 
figure(9)
plot(t_ob,thetaO,'LineWidth',1.5)
hold on
plot(t,theta,'LineWidth',1.5)
xlabel('Tiempo [seg]')
ylabel('ángulo [rad]')
title('Theta')
legend('Con obs', 'Sin obs')
grid on
hold off