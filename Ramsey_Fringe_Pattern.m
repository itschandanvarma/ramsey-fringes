%Author: Chandan Varma Tamada
%Description: Simulates the evolution of a Bloch vector for a Ramsey sequence, consisting of two π/2 pulses separated by free precession, and plots the Ramsey fringe pattern for various δ/Ω ratios

clear; clc; close all;

%Parameters
params.Omega = 1;                         %Rabi frequency
params.tp = pi/(2*params.Omega);          %pi/2 pulse duration
params.T_wait = 20/params.Omega;          %Free evolution time
params.N_pulse = 30;                      %Steps per pi/2 pulse
params.N_wait = 50;                       %Steps during free evolution

%Range of delta/Omega ratios to simulate
delta_omega_ratios = -10:0.01:10;            
final_populations = zeros(size(delta_omega_ratios));

%Rotation Matrix (Rodrigues formula)
rot = @(n,theta) cos(theta)*eye(3) + (1-cos(theta))*(n(:)*n(:)') - sin(theta)*[0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];

for idx = 1:length(delta_omega_ratios)
    %Set detuning for current δ/Ω ratio
    params.delta = delta_omega_ratios(idx) * params.Omega;

    %Effective Rotation Axes
    params.n_eff = [params.Omega, 0, params.delta]/sqrt(params.Omega^2 + params.delta^2);
    params.Omega_eff = sqrt(params.Omega^2 + params.delta^2);
    params.n_z = [0, 0, 1];

    %Initialize Bloch vector
    v = [0; 0; 1];

    %First pi/2 pulse
    for k = 1:params.N_pulse
        R = rot(params.n_eff, params.Omega_eff * (params.tp/params.N_pulse));
        v = R * v;
    end

    %Free Precession
    for k = 1:params.N_wait
        R = rot(params.n_z, params.delta * (params.T_wait/params.N_wait));
        v = R * v;
    end

    %Second pi/2 pulse
    for k = 1:params.N_pulse
        R = rot(params.n_eff, params.Omega_eff * (params.tp/params.N_pulse));
        v = R * v;
    end

    %Store final population
    final_populations(idx) = (1 - v(3)) / 2;
end

%Plot Ramsey fringe pattern
figure('Color', 'w', 'Position', [200 150 600 400]);
plot(delta_omega_ratios, final_populations, 'b-', 'LineWidth', 2);
xlabel('$\delta/\Omega$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$P_g$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
title('Ramsey Fringe Pattern', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);
saveas(gcf, 'Ramsey_Fringe_Pattern.png');