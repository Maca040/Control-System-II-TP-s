%% Solución analítica para un circuito RLC en serie
clear all;
close all;

% Parámetros del circuito RLC
R = 220;        % Resistencia (Ohm)
L = 500e-3;       % Inductancia (H)
C = 2.2e-6;     % Capacitancia (F)
V0 = 0;       % Voltaje inicial en el capacitor (V)

I0 = 0;        % Corriente inicial (A)

% Parámetros de simulación
t0 = 0;        % Tiempo inicial (s)
tf = 0.02;      % Tiempo final (s)
h = 0.00068;      % Paso de tiempo (s)
t = t0:h:tf;   % Vector de tiempo
N = length(t); % Número de pasos

% Cálculo de parámetros analíticos
alpha = R / (2*L);              % Coeficiente de amortiguamiento
omega0 = 1 / sqrt(L*C);         % Frecuencia natural
delta = alpha^2 - omega0^2;     % Discriminante

% Inicialización del vector para la solución analítica
i_analytic = zeros(1, N);

% Solución analítica según el tipo de amortiguamiento
if delta < 0
    % Caso subamortiguado
    omega_d = sqrt(omega0^2 - alpha^2); % Frecuencia amortiguada
    A = V0 / (L * omega_d);            % Constante basada en vC(0) = V0, i(0) = 0
    for n = 1:N
        i_analytic(n) = A * exp(-alpha*t(n)) * sin(omega_d*t(n));
    end
    disp('Caso: Subamortiguado');
elseif delta == 0
    % Caso críticamente amortiguado
    A = V0 / L; % Constante basada en vC(0) = V0, i(0) = 0
    for n = 1:N
        i_analytic(n) = A * t(n) * exp(-alpha*t(n));
    end
    disp('Caso: Críticamente amortiguado');
else
    % Caso sobreamortiguado
    s1 = -alpha + sqrt(delta);
    s2 = -alpha - sqrt(delta);
    A1 = V0 / (L * (s1 - s2)); % Constante basada en vC(0) = V0, i(0) = 0
    A2 = -A1;
    for n = 1:N
        i_analytic(n) = A1 * exp(s1*t(n)) + A2 * exp(s2*t(n));
    end
    disp('Caso: Sobreamortiguado');
end

% Graficar solución analítica
figure;
plot(t, i_analytic, 'r-', 'LineWidth', 2, 'DisplayName', 'Corriente (Analítica)');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
title('Solución Analítica Circuito RLC');
legend;
grid on;