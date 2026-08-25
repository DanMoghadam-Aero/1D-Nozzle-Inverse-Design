


function [P, U, Pstar] = solve_inviscid_nozzle(x, Area, U_inlet, P_outlet)
    % SOLVE_INVISCID_NOZZLE Computes 1D flow using Karimian & Schneider method
    % SOLVE_INVISCID_NOZZLE Computes 1D flow (rho=1 assumed)
    N = length(x);
    
    % Initialization
    U = U_inlet * ones(N, 1);
    P = U_inlet * ones(N, 1);
    U_old = U;
    P_old = P;
    uhat = U_inlet * ones(N-1, 1);
    new_uhat = U_inlet * ones(N-1, 1);
    
    err1 = 0.1;
    err2 = 0.2;

    for t = 1:10
        for iteration = 1:1000 
            Mat_A = zeros(2*N, 2*N); 
            B = zeros(2*N, 1);       
            
            % Node 1
            i = 1;
            Mat_A(i,i) = (Area(2)+Area(1))/uhat(1);
            Mat_A(i,i+2) = -(Area(2)+Area(1))/uhat(1);
            Mat_A(i,i+3) = Area(1)+Area(2);
            B(i) = U_inlet*(3*Area(1)-Area(2)) - 2*uhat(1)*(Area(2)-Area(1));
            
            i = 2;
            Mat_A(i,i) = 1;
            B(i) = U_inlet;
            
            % Node 2
            i = 3;
            Mat_A(i,i-2) = -(Area(2)+Area(1))/uhat(1);
            Mat_A(i,i) = ((Area(1)+Area(2))/uhat(1)) + ((Area(3)+Area(2))/uhat(2));
            Mat_A(i,i+1) = Area(3)-Area(1);
            Mat_A(i,i+2) = -(Area(2)+Area(3))/uhat(2);
            Mat_A(i,i+3) = Area(2)+Area(3);
            B(i) = U_inlet*(Area(2)+Area(1)) + 2*uhat(1)*(Area(2)-Area(1)) - 2*uhat(2)*(Area(3)-Area(2));
            
            i = 4;
            Mat_A(i,i-3) = -0.25*(Area(1)+3*Area(2));
            Mat_A(i,i-1) = 0.25*(2*Area(2)+Area(1)+Area(3));
            Mat_A(i,i) = 0.5*uhat(2)*(Area(2)+Area(3));
            Mat_A(i,i+1) = 0.25*(Area(2)-Area(3));
            B(i) = 0.5*U_inlet*uhat(1)*(Area(2)+Area(1));
            
            % Nodes 3 to N-2
            for k = 3:(N-2)
                i = 2*k - 1; % Mass equation row
                Mat_A(i,i-2) = -(Area(k)+Area(k-1))/uhat(k-1);
                Mat_A(i,i-1) = -(Area(k)+Area(k-1));
                Mat_A(i,i) = ((Area(k-1)+Area(k))/uhat(k-1)) + ((Area(k)+Area(k+1))/uhat(k));
                Mat_A(i,i+1) = Area(k+1)-Area(k-1);
                Mat_A(i,i+2) = -(Area(k)+Area(k+1))/uhat(k);
                Mat_A(i,i+3) = Area(k)+Area(k+1);
                B(i) = 2*uhat(k-1)*(Area(k)-Area(k-1)) - 2*uhat(k)*(Area(k+1)-Area(k));
                
                i = 2*k; % Momentum equation row
                Mat_A(i,i-3) = -0.25*(Area(k-1)+3*Area(k));
                Mat_A(i,i-2) = -0.5*uhat(k-1)*(Area(k)+Area(k-1));
                Mat_A(i,i-1) = 0.25*(2*Area(k)+Area(k+1)+Area(k-1));
                Mat_A(i,i) = 0.5*uhat(k)*(Area(k)+Area(k+1));
                Mat_A(i,i+1) = 0.25*(Area(k)-Area(k+1));
                B(i) = 0;
            end
            
            % Node N-1
            i = 2*N - 3;
            Mat_A(i,i-2) = -(Area(N-1)+Area(N-2))/uhat(N-2);
            Mat_A(i,i-1) = -(Area(N-1)+Area(N-2));
            Mat_A(i,i) = ((Area(N-2)+Area(N-1))/uhat(N-2)) + ((Area(N)+Area(N-1))/uhat(N-1));
            Mat_A(i,i+1) = Area(N)-Area(N-2);
            Mat_A(i,i+3) = Area(N)+Area(N-1);
            B(i) = ((P_outlet*(Area(N-1)+Area(N)))/uhat(N-1)) + 2*uhat(N-2)*(Area(N-1)-Area(N-2)) - 2*uhat(N-1)*(Area(N)-Area(N-1));
            
            i = 2*N - 2;
            Mat_A(i,i-3) = -0.25*(Area(N-2)+3*Area(N-1));
            Mat_A(i,i-2) = -0.5*uhat(N-2)*(Area(N-1)+Area(N-2));
            Mat_A(i,i-1) = 0.25*(2*Area(N-1)+Area(N)+Area(N-2));
            Mat_A(i,i) = 0.5*uhat(N-1)*(Area(N-1)+Area(N));
            B(i) = -P_outlet*0.25*(Area(N-1)-Area(N));
            
            % Node N
            i = 2*N - 1;
            Mat_A(i,i-2) = (Area(N)+Area(N-1))/uhat(N-1);
            Mat_A(i,i-1) = Area(N)+Area(N-1);
            Mat_A(i,i+1) = Area(N-1)-3*Area(N);
            B(i) = ((P_outlet*(Area(N)+Area(N-1)))/uhat(N-1)) - 2*uhat(N-1)*(Area(N)-Area(N-1));
            
            i = 2*N;
            Mat_A(i,i-1) = 1;
            B(i) = P_outlet;
            
            % Solve Matrix
            C = Mat_A \ B; 
            
            % Convergence Check & Updates
            con1 = zeros(N,1);
            con2 = zeros(N,1);
            for j=1:N
                con1(j) = abs((P(j)-C(2*j-1))/C(2*j-1)*100);
                con2(j) = abs((U(j)-C(2*j))/C(2*j)*100);
            end
            
            con1max = max(con1);
            con2max = max(con2);
            
            if con1max > err1 || con2max > err1
                for j=1:N
                    P(j) = C(2*j-1);
                    U(j) = C(2*j);
                end
                for j=1:N-1
                    new_uhat(j) = 0.5*(U(j)+U(j+1)) - 0.5*((P(j+1)-P(j))/uhat(j)) + uhat(j)*(Area(j+1)-Area(j))/(Area(j)+Area(j+1));
                end
                uhat = new_uhat;
            else
                break
            end
        end
        
        for j=1:N
            if abs((P_old(j)-C(2*j-1))/C(2*j-1)*100) > err2 || abs((U_old(j)-C(2*j))/C(2*j)*100) > err2
                P_old(j) = C(2*j-1);
                U_old(j) = C(2*j);
            else
                break
            end
        end
    end

    % Error tracking
    Pstar = abs(((P(N) + 0.5*U(N)*U(N)) - P(1) - 0.5*U(1)*U(1)) / (0.5*U(1)*U(1)));
end