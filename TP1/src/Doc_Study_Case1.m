%% 
% _Ingeniería Electrónica_
%%     *Actividad Práctica N°1*
%%     *Representación de sistemas y controladores* 
% 
% 
% 
% 
% 
%     *Sistema de Control II* 
%     *2025*
% 
% 
% 
%     *Profesor: Julian Pucheta*
%     *Alumna: Macarena V. González*
% 
%% *Caso de estudio 1: Sistema de dos variables de estado*
% Para este caso de estudio se plantea un sistema eléctrico RLC cuyo  comportamiento 
% se encuentra modelado en variables de estado, este modelo puede  obtenerse a 
% partir de las ecuaciones diferenciales que surgen de la aplicación directa de  
% la ley de Kirchhoff de las tensiones a la malla del circuito

% Consignas a resolver →(ítems)
% Ítem[1] 
% Asignar valores a R=220Ω , L=500mHy, y C=2,2μF. Obtener simulaciones que permitan  
% estudiar la dinámica del sistema, con una entrada de tensión escalón de 12V, 
% que cada  1ms cambia  de signo.
% Resolución
% Para estudiar la dinámica del sistema se procede a realizar los siguientes 
% pasos: 

clear all; clc; close all;
%% 
% * Declarar los Parámetros  de diseño y  matrices que forman las ecuaciones 
% de estados  : 

R= 220;
L= 500e-3; 
Cap= 2.2e-6;
V_e = 12;% voltaje de entrada
%Matrices del espacio 
A= [-R/L -1/L ; 1/Cap 0]; % Matriz de estados
B= [1/L ; 0]; %Matriz de entrada 
C= [R 0]; %Matriz de salida
D=[0]; %Matriz de transmisión directa
%% 
% * También es importante declarar el punto de operación,es necesario para encontrar 
% las ecuaciones de estado:

%Punto de operación de V y I
I1(1)=0;
Vc(1)=0;
y(1)=0;
Xop=[0 0]' ;
x=[I1(1) Vc(1)]'; 
%% 
% * Obtener  la FT a partir de las matrices del espacio de estado → me sirve 
% para calcular los polos  y por ende el comportamiento de estabilidad del sistema 
% y los tiempos de simulación


%convierte de espacio de estados a función de transferencia
[numF,denF] = ss2tf(A,B,C,D)
%Función de transferencia del sistema 
F=tf(numF,denF) 
%Polos de la FT
poles=roots(denF) % POLOS complejos conjugados :representan una respuesta temp SUBamortiguada en el tiempo 0<ξ<1

% saco la parte imaginaria del polo = Wd 
Wd1=imag(poles(1)) 
Wd2= imag(poles(2))
%Proporciona ξ y Wn de la FT
[Wn,zita]=damp(F)
%damp(F) me proporciona mas información
%% 
% Como los polos son complejos conjugados, representan una respuesta temp. *sub*amortiguada 
% en el tiempo |*0<ξ<1*| ,  para calcular los tiempos de simulación y de paso 
% se debe considerar la _*frecuencia de oscilación*_  $W_d$ (parte imaginaria 
% de los polos).
% 
% - Para estimar el  tamaño de paso que permita observar una dinámica  más rápida  
% se utiliza el polo  que coincide con el 95% de la dinámica.
% 
% -Al considerar una dinámica rápida el tiempo de integración debería ser aprox  
% 10 veces menor que el tiempo de paso 

%Cálculo de Período de Wd1 
t_d=(2*pi)/Wd1 
%time de integración 
t_int=t_d/10
%----
%% 
% Para el tiempo de simulación se toma la constante de tiempo más lenta p/ 5% 
% :

t_l=log(0.05)/(real(poles(1)))
t_sim=3*t_l
%% 
% Como  se pide una simulación con una entrada que cambie de signo cada 10 ms 
% entonces Tsim debe ser mayor lo cua se cumple ya que nos da  $t_{\textrm{sim}} 
% =0\ldotp 0409$
% 
% |El tamaño de paso o "step" esta definido por:| 

step=round(t_sim/t_int) %cant. de puntos de sim.= time total de simu/time que dura cada tamaño de paso 
%% 
% Sabiendo la cantidad de puntos  necesarios se contruye el vector tiempo : 

t=linspace(0,t_sim,step);
%% 
% Con el vector tiempo  y   el punto de operación se puede definir la señal 
% de entrada: 

%vector inicial de la entrada

t = linspace(0, t_sim, step);

u=linspace(0,0,step);

ii=0;

for i=1:step-1
    ii = ii + t_int;
    if(ii >= 10e-3) %porque conmuta cada 10 milisegundo
        ii=0;
        V_e=V_e*-1;

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