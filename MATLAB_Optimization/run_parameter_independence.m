

%% Parameter Independence Study (Bézier Control Points)
clear; clc; close all;

%% 1. Domain & Target Setup
L = 2.0; N = 63;
x = linspace(0, L, N)';
rho = 1.225;
R_in = 0.8; R_out = 1.0;
U_inlet = 50; P_out_real = 101325;
P_out_kin = P_out_real / rho;

% Rebuild Target Pressure Curve
U_out_target = U_inlet * (R_in^2 / R_out^2);
P_in_target = P_out_real + 0.5 * rho * (U_out_target^2 - U_inlet^2);
x_th = 0.3 * L; P_th_target = 92000; 
P_target_real = zeros(N, 1);
decay_k = 5 / (L - x_th); 
for i = 1:N
    if x(i) <= x_th
        P_target_real(i) = P_th_target + (P_in_target - P_th_target) * (0.5 * (1 + cos(pi * x(i) / x_th)));
    else
        P_target_real(i) = P_out_real - (P_out_real - P_th_target) * exp(-decay_k * (x(i) - x_th));
    end
end
P_targ_kin = P_target_real / rho; 

%% 2. The Independence Loop
test_points = 3:9; % Testing between 3 and 9 intermediate control points
RMSE_results = zeros(length(test_points), 1);

% GA Options (Slightly reduced for speed during the loop)
options = optimoptions('ga', 'Display', 'off', ...
    'PopulationSize', 40, ...
    'MaxGenerations', 75);

fprintf('Starting Parameter Independence Study...\n');

for i = 1:length(test_points)
    n_pts = test_points(i);
    fprintf('Testing %d control points... ', n_pts);
    
    Lower_Bounds = ones(1, n_pts) * 0.3; 
    Upper_Bounds = ones(1, n_pts) * 1.5;
    
    objective_func = @(R_opt) evaluate_nozzle_objective(R_opt, x, L, N, ...
                              U_inlet, P_out_kin, P_targ_kin, R_in, R_out);
                          
    % Run GA
    [~, min_RMSE] = ga(objective_func, n_pts, [], [], [], [], Lower_Bounds, Upper_Bounds, [], options);
    
    RMSE_results(i) = min_RMSE;
    fprintf('RMSE: %.4f\n', min_RMSE);
end

%% 3. Plotting the Pareto Front (Elbow Graph)
figure('Name', 'Parameter Independence Study', 'Position', [200, 200, 600, 450]);
plot(test_points, RMSE_results, '-ko', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');

title('Bezier Parameter Independence (Pareto Front)');
xlabel('Number of Intermediate Control Points');
ylabel('Final Objective Error (RMSE)');
grid on;

% Add a visual marker for our chosen parameter
hold on;
plot(6, RMSE_results(test_points == 6), 'bO', 'MarkerSize', 14, 'LineWidth', 2);
text(6.2, RMSE_results(test_points == 6), 'Selected Parameter (n=6)', 'Color', 'b', 'FontWeight', 'bold');


