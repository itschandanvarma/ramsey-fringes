%Author: Chandan Varma Tamada
%Description: Simulates and animates the evolution of a Bloch vector on a Bloch sphere for a Ramsey sequence, consisting of two π/2 pulses separated by free precession

clear; clc; close all;

%Parameters
params.Omega = 1;                         %Rabi frequency
params.delta = sqrt(15)*params.Omega;     %Detuning
params.tp = pi/(2*params.Omega);          %pi/2 pulse duration
params.T_wait = 20/params.Omega;          %Free evolution time

params.N_pulse = 30;                      %Steps per pi/2 pulse
params.N_wait = 50;                       %Steps during free evolution

%Effective Rotation Axes
params.n_eff = [params.Omega, 0, params.delta]/sqrt(params.Omega^2 + params.delta^2);
params.Omega_eff = sqrt(params.Omega^2 + params.delta^2);
params.n_z = [0, 0, 1];

%Rotation Matrix (Rodrigues formula)
rot = @(n,theta) cos(theta)*eye(3) + (1-cos(theta))*(n(:)*n(:)') - sin(theta)*[0 -n(3) n(2); n(3) 0 -n(1); -n(2) n(1) 0];

%Initialize Bloch vector
v = [0; 0; 1];
traj.traj = [];
traj.times = [];
traj.B = [];

%First pi/2 pulse
for k = 1:params.N_pulse
    R = rot(params.n_eff, params.Omega_eff * (params.tp/params.N_pulse));
    v = R * v;
    traj.traj(:,end+1) = v;
    traj.times(end+1) = (k-1)*(params.tp/params.N_pulse);
    traj.B(:,end+1) = params.n_eff;
end

%Free Precession
for k = 1:params.N_wait
    R = rot(params.n_z, params.delta * (params.T_wait/params.N_wait));
    v = R * v;
    traj.traj(:,end+1) = v;
    traj.times(end+1) = params.tp + (k-1)*(params.T_wait/params.N_wait);
    traj.B(:,end+1) = params.n_z;
end

%Second pi/2 pulse
for k = 1:params.N_pulse
    R = rot(params.n_eff, params.Omega_eff * (params.tp/params.N_pulse));
    v = R * v;
    traj.traj(:,end+1) = v;
    traj.times(end+1) = params.tp + params.T_wait + (k-1)*(params.tp/params.N_pulse);
    traj.B(:,end+1) = params.n_eff;
end

%Animation setup
video = VideoWriter('Ramsey_Bloch_Animation.mp4', 'MPEG-4');
video.FrameRate = 15;
open(video);

figure('Color','w','Position',[200 150 600 600]);
[xs,ys,zs] = sphere(40);
cmap = flipud(jet(256));
cdata = zs;

for i = 1:length(traj.times)
    clf

    %Bloch sphere
    surf(xs, ys, zs, cdata, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceLighting', 'gouraud');
    colormap(cmap);
    hold on; axis equal; view([130 25]);
    xlim([-1.2 1.2]); ylim([-1.2 1.2]); zlim([-1.2 1.2]);
    [x_grid, y_grid, z_grid] = sphere(20);
    plot3(x_grid, y_grid, z_grid, 'Color', [0.5 0.5 0.5, 0.2], 'LineWidth', 0.5); 
    plot3(x_grid', y_grid', z_grid', 'Color', [0.5 0.5 0.5, 0.2], 'LineWidth', 0.5); 
    light('Position', [1 1 1], 'Style', 'infinite');
    lighting gouraud;
    
    %Axes
    xlabel('$X$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('$Y$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
    zlabel('$Z$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
    plot3([-1 1],[0 0],[0 0],'k-', 'LineWidth', 1.5);
    plot3([0 0],[-1 1],[0 0],'k-', 'LineWidth', 1.5);
    plot3([0 0],[0 0],[-1 1],'k-', 'LineWidth', 1.5);
    
    %Ground and excited state markers
    plot3(0,0,1,'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8); 
    plot3(0,0,-1,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    text(0.05, 0, 1.5, '$|g\rangle$', 'Color', 'b', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
    text(0.05, 0, -1.5, '$|e\rangle$', 'Color', 'r', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');

    %Current Bloch vector
    v = traj.traj(:,i);
    h1 = quiver3(0, 0, 0, v(1), v(2), v(3), 0, 'Color', [0 0 0], 'LineWidth', 3, 'MaxHeadSize', 0.9);
    text(v(1)*0.5, 0.1 + v(2)*0.5, v(3)*0.5, '$\bar{n}$', 'Color', [0 0 0], 'FontSize', 12, 'Interpreter', 'latex', 'FontWeight', 'bold');

    %Current B vector
    B = traj.B(:,i);
    h2 = quiver3(0, 0, 0, B(1), B(2), B(3), 0, 'Color', [0.6 0 0.6], 'LineWidth', 2.5, 'MaxHeadSize', 0.7);
    text(B(1)*1.1, B(2)*1.1, 0.15 + B(3)*1.1, '$\bar{B}$', 'Color', [0.6 0 0.6], 'FontSize', 12, 'Interpreter', 'latex', 'FontWeight', 'bold');

    %Trajectory trace
    plot3(traj.traj(1,1:i), traj.traj(2,1:i), traj.traj(3,1:i), 'Color', [0 0 0], 'LineWidth', 1.5);

    title(sprintf('Bloch Vector Evolution ($t \\Omega = %.1f$, $\\delta/\\Omega = %.2f$)', traj.times(i)*params.Omega, params.delta / params.Omega), 'FontSize', 14, 'Interpreter', 'latex', 'FontWeight', 'bold');
    set(gca, 'XTick', [], 'YTick', [], 'ZTick', []);
    
    frame = getframe(gcf);
    writeVideo(video, frame);
end

close(video);