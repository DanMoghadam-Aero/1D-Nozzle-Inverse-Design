%% 1D Nozzle Inverse Design via Genetic Algorithm
clear; clc; close all;

%% 1. Domain & Fluid Setup
L = 2.0; 
N = 63;
x = linspace(0, L, N)';
rho = 1.225;

% Packaging constraints
R_in = 0.8; 
R_out = 1.0;
U_inlet = 50; 
P_out_real = 101325;
P_out_kin = P_out_real / rho;

%% 2. Target Pressure Curve Synthesis (The Aerodynamic Ideal)
U_out_target = U_inlet * (R_in^2 / R_out^2);
P_in_target = P_out_real + 0.5 * rho * (U_out_target^2 - U_inlet^2);

x_th = 0.3 * L;      
P_th_target = 92000; 
P_target_real = zeros(N, 1);
decay_k = 5 / (L - x_th); 

for i = 1:N
    if x(i) <= x_th
        P_target_real(i) = P_th_target + (P_in_target - P_th_target) * ...
                           (0.5 * (1 + cos(pi * x(i) / x_th)));
    else
        P_target_real(i) = P_out_real - (P_out_real - P_th_target) * ...
                           exp(-decay_k * (x(i) - x_th));
    end
end
P_targ_kin = P_target_real / rho; 

%% 3. The "Naive" Bezier Baseline
% A generic, unoptimized guess for the 6 intermediate control points
num_control_points = 6;
R_opt_baseline = [0.60, 0.45, 0.55, 0.70, 0.80, 0.90];

% Generate and solve the baseline for comparison plots later
R_full_baseline = [R_in, R_opt_baseline, R_out];
[~, A_baseline, R_baseline_geom] = generate_bezier_geometry(R_full_baseline, L, N);
[P_kin_baseline, ~, ~] = solve_inviscid_nozzle(x, A_baseline, U_inlet, P_out_kin);
P_real_baseline = P_kin_baseline * rho;

%% 4. Genetic Algorithm Setup
Lower_Bounds = ones(1, num_control_points) * 0.3; 
Upper_Bounds = ones(1, num_control_points) * 1.5;

objective_func = @(R_opt) evaluate_nozzle_objective(R_opt, x, L, N, ...
                          U_inlet, P_out_kin, P_targ_kin, R_in, R_out);

% GA Options: Seed the algorithm with our Naive Baseline!
options = optimoptions('ga', 'Display', 'iter', ...
    'PopulationSize', 50, ...
    'MaxGenerations', 100, ...
    'InitialPopulationMatrix', R_opt_baseline, ... % <--- Seeding the baseline
    'PlotFcn', {@gaplotbestf, @gaplotstopping});

%% 5. Run the Optimization
fprintf('Starting Genetic Algorithm Optimization...\n');
[R_opt_final, min_RMSE] = ga(objective_func, num_control_points, ...
                             [], [], [], [], Lower_Bounds, Upper_Bounds, [], options);

fprintf('Optimization complete. Final RMSE: %.4f\n', min_RMSE);

%% 6. Final Flow Evaluation & Visualization
R_full_final = [R_in, R_opt_final, R_out];
[~, A_final, R_final_geom] = generate_bezier_geometry(R_full_final, L, N);
[P_kin_final, ~, ~] = solve_inviscid_nozzle(x, A_final, U_inlet, P_out_kin);
P_real_final = P_kin_final * rho;

%% 7. Plotting the Before & After
figure('Name', 'Inverse Design Optimization Results', 'Position', [100, 100, 800, 800]);

% Plot 1: Geometry Comparison
ax1 = subplot(2, 1, 1);
x_fill = [x; flipud(x)];
y_fill = [R_final_geom; flipud(-R_final_geom)];
fill(x_fill, y_fill, [0.9 0.9 0.9], 'EdgeColor','none', 'HandleVisibility', 'off'); hold on;

% Baseline Walls (Dashed Gray)
plot(x, R_baseline_geom, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'Naive Baseline');
plot(x, -R_baseline_geom, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Optimized Walls (Solid Black)
plot(x, R_final_geom, '-k', 'LineWidth', 2, 'DisplayName', 'GA Optimized Shape'); 
plot(x, -R_final_geom, '-k', 'LineWidth', 2, 'HandleVisibility', 'off');

%xline(x_th, ':k', 'Target Throat Location', 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
title('Geometric Evolution: Baseline vs. Optimized Bezier Walls');
ylabel('Radius, r [m]');
legend('Location', 'best');
grid on; 
%axis equal; 
xlim([0, L]); 
%ylim([-1.2 1.2]);

% Plot 2: Pressure Overlay
ax2 = subplot(2, 1, 2);
plot(x, P_target_real, '--k', 'LineWidth', 3, 'DisplayName', 'Analytical Target Curve'); hold on;
plot(x, P_real_baseline, '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'DisplayName', 'Baseline Output');
plot(x, P_real_final, '-b', 'LineWidth', 2, 'DisplayName', 'Optimized Output');

title('Aerodynamic Pressure Matching');
xlabel('Axial Length, x [m]'); ylabel('Static Pressure, P [Pa]');
legend('Location', 'best');
grid on; xlim([0, L]);

linkaxes([ax1, ax2], 'x');


