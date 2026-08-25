

function J = evaluate_nozzle_objective(R_opt, x, L, N, U_inlet, P_out_kin, P_targ_kin, R_in, R_out)
    % Reconstruct the full geometry vector 
    R_full = [R_in, R_opt, R_out]; 
    
    % Generate the geometry
    [~, A, ~] = generate_bezier_geometry(R_full, L, N);
    
    try
        % Run the solver
        [P_kin, ~, ~] = solve_inviscid_nozzle(x, A, U_inlet, P_out_kin);
        
        % 1. Aerodynamic Error (RMSE)
        RMSE = sqrt(mean((P_kin - P_targ_kin).^2));
        
        % 2. Geometric Roughness Penalty
        % Calculate the second derivative of the control points to find "waviness"
        % A perfectly smooth curve has a low second derivative. A squiggly one has a high one.
        roughness = sum(diff(R_full, 2).^2); 
        
        % 3. The Combined Objective
        % We multiply the roughness by a weight (lambda) to balance the two goals.
        % If it is still too squiggly, increase this number (e.g., to 5000).
        lambda = 1000; 
        J = RMSE + lambda * roughness;
        
    catch
        % Divergence Penalty
        J = 1e9; 
    end
end