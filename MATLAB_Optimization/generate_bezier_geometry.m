

function [x, A, R] = generate_bezier_geometry(R_full, L, N)
    % GENERATE_BEZIER_GEOMETRY
    % Parametrically generates a 1D nozzle geometry using an n-degree Bezier curve.
    
    x = linspace(0, L, N)';
    R = zeros(N, 1);
    
    n_points = length(R_full);
    n = n_points - 1; % Degree of the polynomial
    
    % Normalized spatial coordinate (0 to 1)
    xi = x / L;
    
    % Calculate Bezier curve using Bernstein polynomials
    for i = 0:n
        % Bernstein basis polynomial
        bernstein = nchoosek(n, i) .* (1 - xi).^(n - i) .* xi.^i;
        R = R + bernstein * R_full(i+1);
    end
    
    % Calculate cross-sectional area
    A = pi * R.^2;
end