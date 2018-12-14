%% Scientific Computing  
% Group assignment for course Scientific Computing 
% By: Nerine Usman & Marianne Schaaphok
% Date: 25-11-2018

clear;
%% Parameters

dimension = 3; 
tol = 1e-4;
maxnorm = size(1,9); 
for p = 2
    n = 2^p; 
    h = 1/n; 
    
    T_1 = -1*diag(ones(n,1),-1) + 4*eye(n+1) -1*diag(ones(n,1),1);

    %% Matrix and vector construction in 2D
    if dimension == 2 
        % Construct matrix A without elimination (ook hier nog even naar de
        % randvoorwaarden kijken) 
        A(1:(n+1),01:(n+1)) = T_1; 
        A(1:n+1, n+2: 2*n+2) = -eye(n+1); 

        for i = 2:n 
            A((i-1)*(n+1)+1:i*(n+1),(i-2)*(n+1)+1:(i-1)*(n+1)) = -eye(n+1); 
            A((i-1)*(n+1)+1:i*(n+1),(i-1)*(n+1)+1:i*(n+1)) = T_1; 
            A((i-1)*(n+1)+1:i*(n+1),(i)*(n+1)+1:(i+1)*(n+1)) = -eye(n+1);
        end 

        A((n)*(n+1)+1:(n+1)*(n+1),(n-1)*(n+1)+1:(n)*(n+1)) = -eye(n+1); 
        A((n)*(n+1)+1:(n+1)*(n+1),(n)*(n+1)+1:(n+1)*(n+1)) = T_1; 
        A = 1/(h^2)*A;
        
        % Construct vector f  without elimination (moeten we even
        % controleren) nog waarden van de randvoorwaarden toevoegen aan f
        % voor symmetry van A (zie ook pg 27)
        f = ones((n+1)^2,1);
        for i = 1:(n+1)
           for j = 1:(n+1)
                x = (i-1)*h;
                y = (j-1)*h; 
                f((i-1)*(n+1)+j) = (x^2 + y^2)*sin(x*y);
           end 
        end
        
        % Construct exact solution (moet nog even goed)
        u_ex = zeros(size(f)); 

    %% Matrix and vector construction in 3D
    elseif dimension == 3 

        fprintf('Not yet implemented\n')

    else 
        fprintf('Please choose dimension 2 or 3\n')
    end 

    %% Direct Solver
    % Compute Cholesky decomposition
    R = chol(A); 

    % Firect solving algorithm for norm(r)<tolerance 
    r = f;
    u = zeros(size(f)); 
    i = 0;
    while norm(r)>tol
        y = R\r; 
        du = R'\y; 
        u = u+du; 
        r = f-A*u; 
        i = i+1; 
        if i>50
           fprintf('Could not converge within 50 iteration\n') 
           break
        end
    end
    error = norm(u - u_ex); 
    maxnorm(p-1) = error; 
end
