%% Generalized 1D Nozzle Baseline Analysis
% Solves incompressible, inviscid flow through a converging-diverging nozzle
% Includes Exact Analytical Validation
clear; clc; close all;

%% 1. User Inputs
% Geometric Parameters
L = 2.0;               % Total axial length [m]
N = 63;                % Number of spatial discretization nodes
R_in = 0.8;            % Inlet radius [m]
R_out = 1.0;           % Outlet radius [m]
R_th = 0.5;            % Throat radius [m]
x_th = 0.3 * L;        % Throat axial location [m]

% Fluid & Boundary Parameters
U_inlet = 50;          % Uniform inlet velocity [m/s]
P_outlet_real = 101325;% Atmospheric back pressure [Pa]
rho = 1.225;           % Fluid density [kg/m^3]

%% 2. Geometry Generation
[x, A_baseline, R_baseline] = generate_hyperbolic_geometry(L, N, R_in, R_out, R_th, x_th);

%% 3. Flow Solver Execution
P_outlet_kin = P_outlet_real / rho; % Convert to kinematic pressure for the solver

fprintf('Executing 1D Inviscid Solver...\n');
try
    [P_kin, U_sol, ~] = solve_inviscid_nozzle(x, A_baseline, U_inlet, P_outlet_kin);
    P_real = P_kin * rho; % Convert kinematic pressure back to static pressure [Pa]
    fprintf('Solver successfully converged.\n');
catch
    error('Solver diverged. Verify matrix assembly stability or check for extreme geometric gradients.');
end

%% 4. Exact Analytical Solution for Validation
% Exact Velocity (Conservation of Mass: Q = U * A)
Q = U_inlet * A_baseline(1);              % Volumetric flow rate [m^3/s]
U_exact = Q ./ A_baseline;                % Exact velocity array [m/s]

% Exact Pressure (Bernoulli's Equation)
U_out_exact = Q / A_baseline(end);        % Exact exit velocity
P_exact = P_outlet_real + 0.5 * rho * (U_out_exact^2 - U_exact.^2); % Exact pressure array [Pa]

% Mass Flux Check (Should be perfectly constant across the domain)
Mass_Flux_Numerical = rho .* U_sol .* A_baseline; 

%% 5. Visualization
figure('Name', 'Nozzle Flow Validation', 'Position', [100, 100, 800, 900]);

% --- Plot 1: Geometry (Top) ---
ax1 = subplot(3, 1, 1);
x_fill = [x; flipud(x)];
y_fill = [R_baseline; flipud(-R_baseline)];
fill(x_fill, y_fill, [0.9 0.9 0.9], 'EdgeColor', 'none'); hold on;
plot(x, R_baseline, '-k', 'LineWidth', 2); 
plot(x, -R_baseline, '-k', 'LineWidth', 2);
xline(x_th, '--k', 'Throat', 'LabelVerticalAlignment', 'bottom');
title('Piecewise Parabolic Geometry');
ylabel('Radius, r [m]');
grid on; 
%axis equal; 
xlim([0, L]); 
%ylim([-1.2 1.2]);

% --- Plot 2: Pressure & Velocity Distribution (Middle) ---
ax2 = subplot(3, 1, 2);

% Left Y-Axis: Pressure
yyaxis left;
plot(x, P_exact, '--k', 'LineWidth', 3); hold on;   % Analytical Background
plot(x, P_real, '-b', 'LineWidth', 1.5);            % Numerical Overlay
ylabel('Static Pressure, P [Pa]');
ax2.YColor = 'b';

% Right Y-Axis: Velocity
yyaxis right;
plot(x, U_exact, '--k', 'LineWidth', 3); hold on;   % Analytical Background
plot(x, U_sol, '-r', 'LineWidth', 1.5);             % Numerical Overlay
ylabel('Axial Velocity, U [m/s]');
ax2.YColor = 'r';

xline(x_th, '--k', 'Throat', 'LabelVerticalAlignment', 'bottom');
title('Flow Kinematics: Numerical vs. Analytical Formulation');
grid on; xlim([0, L]);
legend('Exact Analytical Solution','Location','east')

% --- Plot 3: Mass Flux Verification (Bottom) ---
ax3 = subplot(3, 1, 3);
plot(x, Mass_Flux_Numerical, '-r', 'LineWidth', 2);
% Create a tight 5% window around the mass flux to prove it is flat
y_margin = Mass_Flux_Numerical(1) * 0.05; 
ylim([Mass_Flux_Numerical(1) - y_margin, Mass_Flux_Numerical(1) + y_margin]);


title('Numerical Mass Flux Verification (\rho \cdot U \cdot A)');
xlabel('Axial Length, x [m]'); ylabel('Mass Flux [kg/s]');
grid on; xlim([0, L]);

% --- Lock all axes together ---
linkaxes([ax1, ax2, ax3], 'x');