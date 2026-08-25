
function [x, A, R] = generate_hyperbolic_geometry(L, N, R_in, R_out, R_th, x_th)
    % GENERATE_HYPERBOLIC_GEOMETRY
    % Creates a 1D spatial grid and a piecewise parabolic converging-diverging 
    % nozzle profile ensuring a strictly zero-slope throat.
    
    x = linspace(0, L, N)';
    R = zeros(N, 1);
    
    % Calculate parabolic curvature constants
    k1 = (R_in - R_th) / (0 - x_th)^2;
    k2 = (R_out - R_th) / (L - x_th)^2;
    
    % Build the radius array
    for i = 1:N
        if x(i) <= x_th
            R(i) = R_th + k1 * (x(i) - x_th)^2; % Converging section
        else
            R(i) = R_th + k2 * (x(i) - x_th)^2; % Diverging section
        end
    end
    
    % Calculate cross-sectional area
    A = pi * R.^2;
end