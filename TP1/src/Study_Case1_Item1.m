%% 
% _Ingeniería Electrónica_
%%     Actividad Práctica N°1
%%     Representación de sistemas y controladores
%     *Sistema de Control II* 
%     *Profesor: Julian Pucheta*
%     *Alumna: Macarena V. González*
% 
%% Caso de estudio 1: Sistema de dos variables de estado
% Ítem[1] 
% Asignar valores a R=220Ω , L=500mHy, y C=2,2μF. Obtener simulaciones que permitan  
% estudiar la dinámica del sistema, con una entrada de tensión escalón de 12V, 
% que cada  10ms cambia  de signo.
%% Resolución

clear all; clc; close all;
% * Declarar los Parámetros  de diseño 
R= 220;%[ohm]
L= 500e-3;%[Hy]
Cap= 2.2e-6;%[F]
V_e = 12;% [V]
%Matrices del espacio 
A= [-R/L -1/L ; 1/Cap 0]; % Matriz de estados
B= [1/L ; 0]; %Matriz de entrada 
C= [R 0]; %Matriz de salida
D=[0]; %Matriz de transmisión directa

%  También es importante declarar el punto de operación

%Punto de operación de V y I
I1(1)=0;
Vc(1)=0;
y(1)=0;
Xop=[0 0]' ;
x=[I1(1) Vc(1)]'; 
%% 
%  Obtener  la FT a partir de ecuaciones de estado
[numF,denF] = ss2tf(A,B,C,D)
%Función de transferencia del sistema 
F=tf(numF,denF) 
%Polos de la FT
poles=roots(denF) % POLOS complejos conjugados 
%representan una respuesta temp SUBamortiguada en el tiempo 0<ξ<1

%Cálculo de tiempo de integración y simulación

% saco la parte imaginaria del polo = Wd  :frecuencia de oscilación 
Wd1= imag(poles(1)) 
Wd2= imag(poles(2))
%Proporciona ξ y Wn de la FT
[Wn,zita]=damp(F)
%damp(F) me proporciona mas información

%Cálculo de Período de Wd1 
t_d=(2*pi)/Wd1 
%time de integración 
t_int=t_d/100

% Para el tiempo de simulación se toma la constante de tiempo más lenta p/ 5% 

t_l=log(0.05)/(real(poles(1)))
t_sim=3*t_l
%% 
% Como  se pide una simulación con una entrada que cambie de signo cada 10 ms 
% 
% El tamaño de paso o "step" esta definido por: 

step=round(t_sim/t_int) %cant. de puntos de sim.= time total de simu/time que dura cada tamaño de paso 
%% 
% Sabiendo la cantidad de puntos  necesarios se contruye el vector tiempo : 

t=linspace(0,t_sim,step); 

%vector inicial de la entrada

t = linspace(0, t_sim, step);

u=linspace(0,0,step);

ii=0;

for i=1:step-1
    ii = ii + t_int;
    if(ii >= 10e-3) %porque conmuta cada 10 milisegundo
        ii=0;
        V_e=-V_e;

    end
    u(i)= V_e;

    %Variables de estado del sistema lineal
    xp=A*(x-Xop)+B*u(i);
    x=x+xp*t_int;
    Y=C*x;
    y(i+1)=Y(1);
    I1(i+1)=x(1);
    Vc(i+1)=x(2);
end
u(end)=u(end-1);
%% 
% Gráficas

disp('Simulación completada')

figure(1)
subplot(3,1,1)
plot(t,u,'green','LineWidth',2); title('Tension de entrada , u_t') 
grid on
subplot(3,1,2)
plot(t,Vc,'blue','LineWidth',2); title('Tension en el capacitor, Vc_t') 
grid on
subplot(3,1,3);
plot(t,I1,'red','LineWidth',2); title('Corriente , i_t');
grid on;
drawnow;
disp("¡Simulación finalizada! Revisa la figura generada.")